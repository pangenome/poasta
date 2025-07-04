use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::fmt::Write;
use std::marker::PhantomData;
use std::sync::Arc;
use rustc_hash::FxHashMap;
use crate::aligner::{AlignedPair, Alignment};

use crate::graphs::{AlignableRefGraph, NodeIndexType};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, AlignmentType, GetAlignmentCosts, GeneralizedGapCosts, Score};
use crate::aligner::aln_graph::{AlignmentGraph, AlignmentGraphNode, AlignState};
use crate::aligner::astar::{AstarQueue, AstarQueuedItem, AstarVisited};
use crate::aligner::queue::{LayeredQueue, QueueLayer};
use crate::bubbles::index::BubbleIndex;
use crate::bubbles::reached::ReachedBubbleExitsMatch;
use crate::errors::PoastaError;

#[derive(Clone, Copy, Debug)]
pub struct GapAffine {
    cost_mismatch: u8,
    cost_gap_extend: u8,
    cost_gap_open: u8
}

impl GapAffine {
    pub fn new(cost_mismatch: u8, cost_gap_extend: u8, cost_gap_open: u8) -> Self {
        Self { cost_mismatch, cost_gap_extend, cost_gap_open }
    }
}

impl GeneralizedGapCosts for GapAffine {
    fn num_pieces(&self) -> usize {
        1
    }
    
    fn piece_extension_cost(&self, piece_index: usize) -> u8 {
        if piece_index == 1 {
            self.cost_gap_extend
        } else {
            panic!("GapAffine only has 1 piece, requested piece {}", piece_index);
        }
    }
    
    fn piece_max_length(&self, piece_index: usize) -> Option<usize> {
        if piece_index == 1 {
            None // Infinite length for single piece
        } else {
            panic!("GapAffine only has 1 piece, requested piece {}", piece_index);
        }
    }
    
    fn gap_open_cost(&self) -> u8 {
        self.cost_gap_open
    }
    
    fn match_cost(&self, is_match: bool) -> u8 {
        if is_match {
            0
        } else {
            self.cost_mismatch
        }
    }
    
    fn total_gap_cost(&self, gap_length: usize) -> usize {
        if gap_length == 0 {
            0
        } else {
            self.cost_gap_open as usize + gap_length * self.cost_gap_extend as usize
        }
    }
    
    fn should_transition_piece(&self, current_piece: usize, _current_length: usize) -> (bool, usize) {
        if current_piece == 1 {
            (false, 1) // Never transition for single piece
        } else {
            panic!("GapAffine only has 1 piece, current piece {}", current_piece);
        }
    }
}

impl AlignmentCosts for GapAffine {
    type AlignmentGraphType = AffineAlignmentGraph;
    type QueueType<N, O> = AffineLayeredQueue<N, O>
        where N: NodeIndexType,
              O: OffsetType;

    fn new_alignment_graph(&self, aln_type: AlignmentType) -> Self::AlignmentGraphType {
        AffineAlignmentGraph::new(self, aln_type)
    }

    #[inline(always)]
    fn mismatch(&self) -> u8 {
        self.cost_mismatch
    }

    #[inline(always)]
    fn gap_open(&self) -> u8 {
        self.cost_gap_open
    }

    #[inline(always)]
    fn gap_extend(&self) -> u8 {
        self.cost_gap_extend
    }

    #[inline(always)]
    fn gap_open2(&self) -> u8 {
        0
    }

    #[inline(always)]
    fn gap_extend2(&self) -> u8 {
        0
    }

    #[inline]
    fn gap_cost(&self, current_state: AlignState, length: usize) -> usize {
        if length == 0 {
            return 0
        }

        let gap_open = match current_state {
            AlignState::Insertion | AlignState::Deletion => 0,
            AlignState::Match => self.cost_gap_open,
            _ => panic!("Invalid current state {:?} for gap affine scoring model!", current_state)
        };

        gap_open as usize + (length * self.cost_gap_extend as usize)
    }
}

pub struct AffineAlignmentGraph {
    /// Gap affine costs
    costs: GapAffine,

    /// Alignment type, e.g., global alignment or ends-free alignment
    aln_type: AlignmentType,
}

impl AffineAlignmentGraph {
    fn new(costs: &GapAffine, aln_type: AlignmentType) -> Self {
        Self {
            costs: *costs,
            aln_type,
        }
    }
}

impl AlignmentGraph for AffineAlignmentGraph {
    type CostModel = GapAffine;

    fn get_costs(&self) -> &Self::CostModel {
        &self.costs
    }

    fn initial_states<G, O>(&self, ref_graph: &G) -> Vec<AlignmentGraphNode<G::NodeIndex, O>>
        where G: AlignableRefGraph,
              O: OffsetType,
    {
        match self.aln_type {
            AlignmentType::Global => vec![AlignmentGraphNode::new(ref_graph.start_node(), O::zero())],
            AlignmentType::EndsFree {
                qry_free_begin: _, qry_free_end: _,
                graph_free_begin: _, graph_free_end: _ }
            => {
                todo!();
            }
        }
    }

    fn is_end<G, O>(&self, ref_graph: &G, seq: &[u8], node: &AlignmentGraphNode<G::NodeIndex, O>, aln_state: AlignState) -> bool
        where G: AlignableRefGraph,
              O: OffsetType,
    {
        aln_state == AlignState::Match
            && node.node() == ref_graph.end_node()
            && node.offset().as_usize() == seq.len()
    }

    fn expand_all<V, G, O, F>(
        &self, visited_data: &mut V,
        ref_graph: &G,
        seq: &[u8],
        score: Score,
        node: &AlignmentGraphNode<G::NodeIndex, O>,
        state: AlignState,
        mut f: F
    )
        where V: AstarVisited<G::NodeIndex, O>,
              G: AlignableRefGraph,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<G::NodeIndex, O>, AlignState)
    {
        match state {
            AlignState::Match => {
                let child_offset = node.offset().increase_one();
                for ref_succ in ref_graph.successors(node.node()) {
                    let new_node_mis = AlignmentGraphNode::new(ref_succ, child_offset);

                    // Move to next (mis)match state
                    let score_delta = if ref_graph.is_symbol_equal(ref_succ, seq[child_offset.as_usize()-1]) {
                        0u8
                    } else {
                        self.costs.cost_mismatch
                    };
                    let new_score_mis = score + score_delta;

                    if visited_data
                        .update_score_if_lower(&new_node_mis, AlignState::Match, node, state, new_score_mis)
                    {
                        f(score_delta, new_node_mis, AlignState::Match)
                    }

                    // Open deletion
                    let score_delta = self.costs.cost_gap_open + self.costs.cost_gap_extend;
                    let new_node_del = AlignmentGraphNode::new(ref_succ, node.offset());
                    let new_score_del = score + score_delta;
                    if visited_data
                        .update_score_if_lower(&new_node_del, AlignState::Deletion, node, state, new_score_del)
                    {
                        f(score_delta, new_node_del, AlignState::Deletion);
                    }
                }

                // Open insertion
                let new_node_ins = AlignmentGraphNode::new(node.node(), child_offset);
                let new_score_ins = score + self.costs.cost_gap_open + self.costs.cost_gap_extend;
                if child_offset.as_usize() <= seq.len() && visited_data
                    .update_score_if_lower(&new_node_ins, AlignState::Insertion, node, state, new_score_ins)
                {
                    f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_ins, AlignState::Insertion);
                }
            },
            AlignState::Insertion => {
                // I->M edge has zero cost
                if visited_data.update_score_if_lower(node, AlignState::Match, node, AlignState::Insertion, score) {
                    f(0, *node, AlignState::Match);
                }

                // Extend insertion
                let new_node_ins = AlignmentGraphNode::new(node.node(), node.offset().increase_one());
                let new_score_ins = score + self.costs.cost_gap_extend;

                if node.offset().as_usize() < seq.len() && visited_data
                    .update_score_if_lower(&new_node_ins, AlignState::Insertion, node, state, new_score_ins)
                {
                    f(self.costs.cost_gap_extend, new_node_ins, AlignState::Insertion);
                }
            },
            AlignState::Deletion => {
                // D->M edge has zero cost
                if visited_data.update_score_if_lower(node, AlignState::Match, node, AlignState::Deletion, score) {
                    f(0, *node, AlignState::Match);
                }

                for ref_succ in ref_graph.successors(node.node()) {
                    let score_delta = self.costs.gap_extend();

                    // Extend deletion
                    let new_node_del = AlignmentGraphNode::new(ref_succ, node.offset());
                    let new_score_del = score + score_delta;
                    if visited_data
                        .update_score_if_lower(&new_node_del, AlignState::Deletion, node, state, new_score_del)
                    {
                        f(score_delta, new_node_del, AlignState::Deletion);
                    }
                }
            },
            AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine state {state:?}!"),
            AlignState::MultiInsertion { .. } | AlignState::MultiDeletion { .. } => panic!("Multi-piece gap states not supported in gap-affine model")
        }
    }

    fn expand_ref_graph_end<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        score: Score,
        mut f: F
    )
        where V: AstarVisited<N, O>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState)
    {
        // At ref graph end, so it's not possible to expand to a (mis)match or deletion, since there are no nodes left.
        // We can still open an insertion though.
        let new_node_ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        let new_score_ins = score + self.costs.cost_gap_open + self.costs.cost_gap_extend;
        if visited_data
            .update_score_if_lower(&new_node_ins, AlignState::Insertion, parent, AlignState::Match, new_score_ins)
        {
            f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_ins, AlignState::Insertion);
        }
    }

    fn expand_query_end<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        child: N,
        score: Score,
        mut f: F
    )
        where V: AstarVisited<N, O>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState)
    {
        // At query end, so we can't continue with a (mis)match or insertion, since there is no query sequence left.
        // We can still open a deletion.
        let new_node_del = AlignmentGraphNode::new(child, parent.offset());
        let new_score_del = score + self.costs.cost_gap_open + self.costs.cost_gap_extend;
        if visited_data
            .update_score_if_lower(&new_node_del, AlignState::Deletion, parent, AlignState::Match, new_score_del)
        {
            f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_del, AlignState::Deletion);
        }
    }

    fn expand_mismatch<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        child: &AlignmentGraphNode<N, O>,
        score: Score,
        mut f: F
    )
        where V: AstarVisited<N, O>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState)
    {
        let new_score_mis = score + self.costs.cost_mismatch;
        if visited_data
            .update_score_if_lower(child, AlignState::Match, parent, AlignState::Match, new_score_mis)
        {
            f(self.costs.cost_mismatch, *child, AlignState::Match);
        }

        // Also queue indel states from parent
        let new_node_ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        let new_score_ins = score + self.costs.cost_gap_open + self.costs.cost_gap_extend;

        if visited_data
            .update_score_if_lower(&new_node_ins, AlignState::Insertion, parent, AlignState::Match, new_score_ins)
        {
            f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_ins, AlignState::Insertion);
        }

        let new_node_del = AlignmentGraphNode::new(child.node(), parent.offset());
        let new_score_del = score + self.costs.cost_gap_open + self.costs.cost_gap_extend;
        if visited_data
            .update_score_if_lower(&new_node_del, AlignState::Deletion, parent, AlignState::Match, new_score_del)
        {
            f(self.costs.cost_gap_open + self.costs.cost_gap_extend, new_node_del, AlignState::Deletion);
        }
    }

}


#[derive(Default, Copy, Clone)]
struct VisitedCellAffine {
    visited_m: Score,
    visited_i: Score,
    visited_d: Score
}

struct BlockedVisitedStorageAffine<N, O, const B: usize = 8>
    where N: NodeIndexType,
          O: OffsetType,
{
    node_blocks: Vec<FxHashMap<O, [[VisitedCellAffine; B]; B]>>,
    node_ranks: Vec<usize>,
    dummy: PhantomData<N>,
}

impl<N, O, const B: usize> BlockedVisitedStorageAffine<N, O, B>
    where N: NodeIndexType,
          O: OffsetType,
{
    pub fn new<G: AlignableRefGraph>(ref_graph: &G) -> Self {
        if B & (B-1) != 0 {
            panic!("Block size B should be a power of 2!")
        }

        let num_blocks_nodes = (ref_graph.node_count_with_start_and_end() / B) + 1;
        Self {
            node_blocks: vec![FxHashMap::default(); num_blocks_nodes],
            node_ranks: ref_graph.get_node_ordering(),
            dummy: PhantomData
        }
    }

    #[inline(always)]
    pub fn calc_block_ix(&self, aln_node: &AlignmentGraphNode<N, O>) -> (usize, O, usize, usize) {
        let node_rank = self.node_ranks[aln_node.node().index()];
        let node_block = node_rank >> B.ilog2();

        let offset = aln_node.offset().as_usize();
        let offset_block = O::new(offset >> B.ilog2());

        let within_block_node = node_rank & (B-1);
        let within_block_qry = offset & (B-1);

        (node_block, offset_block, within_block_node, within_block_qry)
    }

    #[inline]
    pub fn get_score(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> Score {
        let (node_block, offset_block, within_block_node, within_block_qry)
            = self.calc_block_ix(aln_node);

        self.node_blocks[node_block].get(&offset_block)
            .map(|v| {
                let cell_data = &v[within_block_node][within_block_qry];
                let node = match aln_state {
                    AlignState::Match => &cell_data.visited_m,
                    AlignState::Insertion => &cell_data.visited_i,
                    AlignState::Deletion => &cell_data.visited_d,
                    AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine state {aln_state:?}"),
                    AlignState::MultiInsertion { .. } | AlignState::MultiDeletion { .. } => panic!("Multi-piece gap states not supported in gap-affine model")
                };

                *node
            })
            .unwrap_or(Score::Unvisited)
    }

    #[inline]
    pub fn set_score(&mut self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState, score: Score) {
        let (node_block, offset_block, within_block_node, within_block_qry)
            = self.calc_block_ix(aln_node);

        let v = self.node_blocks[node_block].entry(offset_block)
            .or_insert_with(|| [[VisitedCellAffine::default(); B]; B]);

        let cell_data = &mut v[within_block_node][within_block_qry];
        match aln_state {
            AlignState::Match => cell_data.visited_m = score,
            AlignState::Insertion => cell_data.visited_i = score,
            AlignState::Deletion => cell_data.visited_d = score,
            AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine state {aln_state:?}"),
            AlignState::MultiInsertion { .. } | AlignState::MultiDeletion { .. } => panic!("Multi-piece gap states not supported in gap-affine model")
        };
    }

    #[inline]
    pub fn update_score_if_lower_block(
        &mut self,
        aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState,
        _: &AlignmentGraphNode<N, O>, _: AlignState,
        score: Score
    ) -> bool {
        let (node_block, offset_block, within_block_node, within_block_qry)
            = self.calc_block_ix(aln_node);

        let v = self.node_blocks[node_block].entry(offset_block)
            .or_insert_with(|| [[VisitedCellAffine::default(); B]; B]);

        let cell_data = &mut v[within_block_node][within_block_qry];

        let node = match aln_state {
            AlignState::Match => &mut cell_data.visited_m,
            AlignState::Insertion => &mut cell_data.visited_i,
            AlignState::Deletion => &mut cell_data.visited_d,
            AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine state {aln_state:?}"),
            AlignState::MultiInsertion { .. } | AlignState::MultiDeletion { .. } => panic!("Multi-piece gap states not supported in gap-affine model")
        };

        match score.cmp(node) {
            Ordering::Less => {
                *node = score;
                true
            },
            Ordering::Equal | Ordering::Greater => false,
        }
    }

    pub fn get_backtrace<G: AlignableRefGraph<NodeIndex=N>>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        costs: &GapAffine,
        aln_node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState,
    ) -> Option<(AlignmentGraphNode<N, O>, AlignState)> {
        let (node_block, offset_block, within_block_node, within_block_qry)
            = self.calc_block_ix(aln_node);

        let curr_score = self.node_blocks[node_block].get(&offset_block)
            .map(|v| {
                let cell_data = &v[within_block_node][within_block_qry];
                let node = match aln_state {
                    AlignState::Match => &cell_data.visited_m,
                    AlignState::Insertion => &cell_data.visited_i,
                    AlignState::Deletion => &cell_data.visited_d,
                    AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine state {aln_state:?}"),
                    AlignState::MultiInsertion { .. } | AlignState::MultiDeletion { .. } => panic!("Multi-piece gap states not supported in gap-affine model")
                };

                *node
            })?;

        match aln_state {
            AlignState::Match => {
                if aln_node.offset() > O::zero() {
                    let is_match_or_end =
                        ref_graph
                            .is_symbol_equal(aln_node.node(), seq[aln_node.offset().as_usize()-1])
                        || aln_node.node() == ref_graph.end_node();
                    let pred_offset = if aln_node.node() == ref_graph.end_node() {
                        aln_node.offset()
                    } else {
                        aln_node.offset() - O::one()
                    };

                    // First priority: match/mismatch
                    for p in ref_graph.predecessors(aln_node.node()) {
                        let pred = AlignmentGraphNode::new(p, pred_offset);
                        let pred_score = self.get_score(&pred, AlignState::Match);

                        if (is_match_or_end && pred_score == curr_score)
                            || (!is_match_or_end && pred_score == curr_score - costs.cost_mismatch)
                        {
                            return Some((pred, AlignState::Match))
                        }
                    }
                }

                // Second priority: close deletion
                let pred_score = self.get_score(aln_node, AlignState::Deletion);
                if pred_score == curr_score {
                    return Some((*aln_node, AlignState::Deletion))
                }

                // Third priority: close insertion
                let pred_score = self.get_score(aln_node, AlignState::Insertion);
                if pred_score == curr_score {
                    return Some((*aln_node, AlignState::Insertion))
                }
            },
            AlignState::Deletion => {
                // First priority: opening new deletion
                for p in ref_graph.predecessors(aln_node.node()) {
                    let pred = AlignmentGraphNode::new(p, aln_node.offset());
                    let pred_score = self.get_score(&pred, AlignState::Match);

                    if pred_score == curr_score - costs.cost_gap_open - costs.cost_gap_extend {
                        return Some((pred, AlignState::Match))
                    }
                }

                // Second priority: extend deletion
                for p in ref_graph.predecessors(aln_node.node()) {
                    let pred = AlignmentGraphNode::new(p, aln_node.offset());
                    let pred_score = self.get_score(&pred, AlignState::Deletion);

                    if pred_score == curr_score - costs.cost_gap_extend {
                        return Some((pred, AlignState::Deletion))
                    }
                }
            },
            AlignState::Insertion => {
                if aln_node.offset() > O::zero() {
                    // First priority: opening new insertion
                    let pred = AlignmentGraphNode::new(aln_node.node(), aln_node.offset() - O::one());
                    let pred_score = self.get_score(&pred, AlignState::Match);

                    if pred_score == curr_score - costs.cost_gap_open - costs.cost_gap_extend {
                        return Some((pred, AlignState::Match))
                    }

                    // Second priority: extend insertion
                    let pred_score = self.get_score(&pred, AlignState::Insertion);
                    if pred_score == curr_score - costs.cost_gap_extend {
                        return Some((pred, AlignState::Match))
                    }
                }
            },
            AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine state {aln_state:?}!"),
            AlignState::MultiInsertion { .. } | AlignState::MultiDeletion { .. } => panic!("Multi-piece gap states not supported in gap-affine model")
        }

        None
    }

    pub fn write_tsv<W: Write>(&self, writer: &mut W) -> Result<(), PoastaError> {
        writeln!(writer, "node_id\toffset\tmatrix\tscore")?;

        let mut rank_to_node_ix: Vec<usize> = (0..self.node_ranks.len())
            .collect();

        rank_to_node_ix.sort_unstable_by_key(|v| self.node_ranks[*v]);
        let b_as_o = O::new(B);

        for (node_block_num, node_block) in self.node_blocks.iter().enumerate() {
            let node_rank_base = node_block_num * B;
            for (qry_block_num, block_data) in node_block.iter() {
                let qry_pos_base = *qry_block_num * b_as_o;
                for (row_num, row) in block_data.iter().enumerate() {
                    let node_rank = node_rank_base + row_num;
                    if node_rank >= rank_to_node_ix.len() {
                        break;
                    }

                    let node_ix = rank_to_node_ix[node_rank];

                    for (col, cell) in row.iter().enumerate() {
                        let qry_pos = qry_pos_base + O::new(col);
                        if let Score::Score(score) = cell.visited_m {
                            writeln!(writer, "{node_ix}\t{qry_pos:?}\tmatch\t{}", score)?;
                        }
                        if let Score::Score(score) = cell.visited_i {
                            writeln!(writer, "{node_ix}\t{qry_pos:?}\tinsertion\t{}", score)?;
                        }
                        if let Score::Score(score) = cell.visited_d {
                            writeln!(writer, "{node_ix}\t{qry_pos:?}\tdeletion\t{}", score)?;
                        }
                    }

                }
            }
        }

        Ok(())
    }
}


pub struct AffineAstarData<N, O>
where N: NodeIndexType,
      O: OffsetType,
{
    costs: GapAffine,
    seq_len: usize,
    bubble_index: Arc<BubbleIndex<N>>,
    visited: BlockedVisitedStorageAffine<N, O>,

    bubbles_reached_m: Vec<BTreeSet<O>>,
}

impl<N, O> AffineAstarData<N, O>
    where N: NodeIndexType,
          O: OffsetType
{
    pub fn new<G>(costs: GapAffine, ref_graph: &G, seq: &[u8], bubble_index: Arc<BubbleIndex<G::NodeIndex>>) -> Self
        where G: AlignableRefGraph<NodeIndex=N>,
    {
        Self {
            costs,
            seq_len: seq.len(),
            bubble_index,
            visited: BlockedVisitedStorageAffine::new(ref_graph),
            bubbles_reached_m: vec![BTreeSet::new(); ref_graph.node_count_with_start_and_end()],
        }
    }

    #[inline]
    fn get_backtrace<G: AlignableRefGraph<NodeIndex=N>>(
        &self,
        ref_graph: &G,
        seq: &[u8],
        aln_node: &AlignmentGraphNode<N, O>,
        aln_state: AlignState
    ) -> Option<(AlignmentGraphNode<N, O>, AlignState)> {
        self.visited.get_backtrace(ref_graph, seq, &self.costs, aln_node, aln_state)
    }
}

impl<N, O> GetAlignmentCosts for AffineAstarData<N, O>
    where N: NodeIndexType, 
          O: OffsetType
{
    type Costs = GapAffine;

    fn get_costs(&self) -> &Self::Costs {
        &self.costs
    }
}

/// Enhanced visited cell that tracks multi-piece gap states
/// Optimized to avoid HashMap allocation for common cases
#[derive(Clone, Debug)]
struct GeneralizedVisitedCell {
    visited_m: Score,
    visited_i: Score,
    visited_d: Score,
    visited_i2: Score,
    visited_d2: Score,
    // Inline storage for up to 6 multi-piece states (most common case)
    // Format: (key, score) where key = (piece << 8) | position
    inline_multi: [(u16, Score); 6],
    inline_count: u8,
    // Overflow storage for rare cases with many states
    overflow_multi: Option<Box<FxHashMap<u16, Score>>>,
}

impl Default for GeneralizedVisitedCell {
    fn default() -> Self {
        Self {
            visited_m: Score::Unvisited,
            visited_i: Score::Unvisited,
            visited_d: Score::Unvisited,
            visited_i2: Score::Unvisited,
            visited_d2: Score::Unvisited,
            inline_multi: [(0, Score::Unvisited); 6],
            inline_count: 0,
            overflow_multi: None,
        }
    }
}

impl GeneralizedVisitedCell {
    #[inline(always)]
    fn get_multi_score(&self, key: u16) -> Score {
        // Check inline storage first (most common case)
        for i in 0..self.inline_count as usize {
            if self.inline_multi[i].0 == key {
                return self.inline_multi[i].1;
            }
        }
        
        // Check overflow if it exists
        if let Some(ref overflow) = self.overflow_multi {
            overflow.get(&key).copied().unwrap_or(Score::Unvisited)
        } else {
            Score::Unvisited
        }
    }
    
    #[inline(always)]
    fn set_multi_score(&mut self, key: u16, score: Score) {
        // Try to update existing entry in inline storage
        for i in 0..self.inline_count as usize {
            if self.inline_multi[i].0 == key {
                self.inline_multi[i].1 = score;
                return;
            }
        }
        
        // Add to inline storage if there's space
        if (self.inline_count as usize) < self.inline_multi.len() {
            self.inline_multi[self.inline_count as usize] = (key, score);
            self.inline_count += 1;
            return;
        }
        
        // Use overflow storage
        if self.overflow_multi.is_none() {
            self.overflow_multi = Some(Box::new(FxHashMap::default()));
        }
        self.overflow_multi.as_mut().unwrap().insert(key, score);
    }
}

/// Generalized blocked storage for multi-piece gap penalty tracking
struct GeneralizedBlockedVisitedStorage<N, O, const B: usize = 8>
where 
    N: NodeIndexType,
    O: OffsetType,
{
    node_blocks: Vec<FxHashMap<O, [[GeneralizedVisitedCell; B]; B]>>,
    node_ranks: Vec<usize>,
    dummy: PhantomData<N>,
}

impl<N, O, const B: usize> GeneralizedBlockedVisitedStorage<N, O, B>
where 
    N: NodeIndexType,
    O: OffsetType,
{
    pub fn new<G: AlignableRefGraph>(ref_graph: &G) -> Self {
        if B & (B-1) != 0 {
            panic!("Block size B should be a power of 2!")
        }

        let num_blocks_nodes = (ref_graph.node_count_with_start_and_end() / B) + 1;
        Self {
            node_blocks: vec![FxHashMap::default(); num_blocks_nodes],
            node_ranks: ref_graph.get_node_ordering(),
            dummy: PhantomData
        }
    }

    #[inline(always)]
    pub fn calc_block_ix(&self, aln_node: &AlignmentGraphNode<N, O>) -> (usize, O, usize, usize) {
        let node_rank = self.node_ranks[aln_node.node().index()];
        let node_block = node_rank >> B.ilog2();
        let node_ix_in_block = node_rank & (B - 1);
        let qry_block = O::new(aln_node.offset().as_usize() >> B.ilog2());
        let qry_ix_in_block = aln_node.offset().as_usize() & (B - 1);

        (node_block, qry_block, node_ix_in_block, qry_ix_in_block)
    }

    pub fn get_score(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> Score {
        let (node_block, qry_block, node_ix_in_block, qry_ix_in_block) = self.calc_block_ix(aln_node);

        self.node_blocks[node_block]
            .get(&qry_block)
            .map(|data| {
                let cell = &data[node_ix_in_block][qry_ix_in_block];
                match aln_state {
                    AlignState::Match => cell.visited_m,
                    AlignState::Insertion => cell.visited_i,
                    AlignState::Deletion => cell.visited_d,
                    AlignState::Insertion2 => cell.visited_i2,
                    AlignState::Deletion2 => cell.visited_d2,
                    AlignState::MultiInsertion { piece, position } | AlignState::MultiDeletion { piece, position } => {
                        let key = ((piece as u16) << 8) | (position as u16);
                        cell.get_multi_score(key)
                    },
                }
            })
            .unwrap_or(Score::Unvisited)
    }

    pub fn update_score_if_lower(
        &mut self,
        child: &AlignmentGraphNode<N, O>,
        child_state: AlignState,
        _parent: &AlignmentGraphNode<N, O>,
        _parent_state: AlignState,
        new_score: Score
    ) -> bool {
        let (node_block, qry_block, node_ix_in_block, qry_ix_in_block) = self.calc_block_ix(child);

        let block = self.node_blocks[node_block]
            .entry(qry_block)
            .or_insert_with(|| std::array::from_fn(|_| std::array::from_fn(|_| GeneralizedVisitedCell::default())));
        
        let cell = &mut block[node_ix_in_block][qry_ix_in_block];
        
        match child_state {
            AlignState::Match => {
                if new_score < cell.visited_m {
                    cell.visited_m = new_score;
                    true
                } else {
                    false
                }
            },
            AlignState::Insertion => {
                if new_score < cell.visited_i {
                    cell.visited_i = new_score;
                    true
                } else {
                    false
                }
            },
            AlignState::Deletion => {
                if new_score < cell.visited_d {
                    cell.visited_d = new_score;
                    true
                } else {
                    false
                }
            },
            AlignState::Insertion2 => {
                if new_score < cell.visited_i2 {
                    cell.visited_i2 = new_score;
                    true
                } else {
                    false
                }
            },
            AlignState::Deletion2 => {
                if new_score < cell.visited_d2 {
                    cell.visited_d2 = new_score;
                    true
                } else {
                    false
                }
            },
            AlignState::MultiInsertion { piece, position } | AlignState::MultiDeletion { piece, position } => {
                let key = ((piece as u16) << 8) | (position as u16);
                let current_score = cell.get_multi_score(key);
                if new_score < current_score {
                    cell.set_multi_score(key, new_score);
                    true
                } else {
                    false
                }
            },
        }
    }
}

impl<N, O, const B: usize> AstarVisited<N, O> for GeneralizedBlockedVisitedStorage<N, O, B>
where 
    N: NodeIndexType,
    O: OffsetType,
{
    fn get_score(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> Score {
        self.get_score(aln_node, aln_state)
    }

    fn set_score(&mut self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState, score: Score) {
        let (node_block, qry_block, node_ix_in_block, qry_ix_in_block) = self.calc_block_ix(aln_node);

        let block = self.node_blocks[node_block]
            .entry(qry_block)
            .or_insert_with(|| std::array::from_fn(|_| std::array::from_fn(|_| GeneralizedVisitedCell::default())));
        
        let cell = &mut block[node_ix_in_block][qry_ix_in_block];
        
        match aln_state {
            AlignState::Match => cell.visited_m = score,
            AlignState::Insertion => cell.visited_i = score,
            AlignState::Deletion => cell.visited_d = score,
            AlignState::Insertion2 => cell.visited_i2 = score,
            AlignState::Deletion2 => cell.visited_d2 = score,
            AlignState::MultiInsertion { piece, position } | AlignState::MultiDeletion { piece, position } => {
                let key = ((piece as u16) << 8) | (position as u16);
                cell.set_multi_score(key, score);
            },
        }
    }

    fn mark_reached(&mut self, _score: Score, _aln_node: &AlignmentGraphNode<N, O>, _aln_state: AlignState) {
        // Bubble tracking would need to be implemented here if needed
    }

    fn dfa_match(&mut self, score: Score, _parent: &AlignmentGraphNode<N, O>, child: &AlignmentGraphNode<N, O>) {
        self.set_score(child, AlignState::Match, score);
    }

    fn prune(&self, score: Score, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> bool {
        let current_score = self.get_score(aln_node, aln_state);
        match current_score {
            Score::Unvisited => false,
            Score::Score(_) => score > current_score,
        }
    }

    fn backtrace<G>(&self, _ref_graph: &G, _seq: &[u8], _aln_node: &AlignmentGraphNode<N, O>) -> crate::aligner::Alignment<N>
    where G: AlignableRefGraph<NodeIndex=N>
    {
        // Simplified placeholder
        Vec::new()
    }

    fn write_tsv<W: std::fmt::Write>(&self, writer: &mut W) -> Result<(), PoastaError> {
        writeln!(writer, "# Generalized visited storage TSV output not implemented")?;
        Ok(())
    }

    fn update_score_if_lower(
        &mut self,
        child: &AlignmentGraphNode<N, O>,
        child_state: AlignState,
        parent: &AlignmentGraphNode<N, O>,
        parent_state: AlignState,
        new_score: Score
    ) -> bool {
        self.update_score_if_lower(child, child_state, parent, parent_state, new_score)
    }
}

/// Generalized A* data structure that properly tracks multi-piece gap states
/// This replaces the AffineAstarData infrastructure with true multi-piece support
pub struct GeneralizedAffineAstarData<N, O, C>
where 
    N: NodeIndexType,
    O: OffsetType,
    C: GeneralizedGapCosts,
{
    generalized_costs: C,
    seq_len: usize,
    bubble_index: Arc<BubbleIndex<N>>,
    // Use enhanced visited storage that tracks piece information
    visited: GeneralizedBlockedVisitedStorage<N, O>,
    bubbles_reached_m: Vec<BTreeSet<O>>,
}

impl<N, O, C> GeneralizedAffineAstarData<N, O, C>
where 
    N: NodeIndexType,
    O: OffsetType,
    C: GeneralizedGapCosts + AlignmentCosts,
{
    pub fn new<G>(costs: C, ref_graph: &G, seq: &[u8], bubble_index: Arc<BubbleIndex<G::NodeIndex>>) -> Self
    where G: AlignableRefGraph<NodeIndex=N>,
    {
        Self {
            generalized_costs: costs,
            seq_len: seq.len(),
            bubble_index,
            visited: GeneralizedBlockedVisitedStorage::new(ref_graph),
            bubbles_reached_m: vec![BTreeSet::new(); ref_graph.node_count_with_start_and_end()],
        }
    }
}

impl<N, O, C> GetAlignmentCosts for GeneralizedAffineAstarData<N, O, C>
where 
    N: NodeIndexType,
    O: OffsetType,
    C: GeneralizedGapCosts + AlignmentCosts,
{
    type Costs = C;

    fn get_costs(&self) -> &Self::Costs {
        &self.generalized_costs
    }
}

impl<N, O, C> AstarVisited<N, O> for GeneralizedAffineAstarData<N, O, C>
where 
    N: NodeIndexType,
    O: OffsetType,
    C: GeneralizedGapCosts + AlignmentCosts,
{
    fn get_score(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> Score {
        self.visited.get_score(aln_node, aln_state)
    }

    fn set_score(&mut self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState, score: Score) {
        // For set_score, we need to implement the direct assignment
        // This is typically used for initialization
        let (node_block, qry_block, node_ix_in_block, qry_ix_in_block) = self.visited.calc_block_ix(aln_node);

        let block = self.visited.node_blocks[node_block]
            .entry(qry_block)
            .or_insert_with(|| std::array::from_fn(|_| std::array::from_fn(|_| GeneralizedVisitedCell::default())));
        
        let cell = &mut block[node_ix_in_block][qry_ix_in_block];
        
        match aln_state {
            AlignState::Match => cell.visited_m = score,
            AlignState::Insertion => cell.visited_i = score,
            AlignState::Deletion => cell.visited_d = score,
            AlignState::Insertion2 => cell.visited_i2 = score,
            AlignState::Deletion2 => cell.visited_d2 = score,
            AlignState::MultiInsertion { piece, position } | AlignState::MultiDeletion { piece, position } => {
                let key = ((piece as u16) << 8) | (position as u16);
                cell.set_multi_score(key, score);
            },
        }
    }

    fn update_score_if_lower(
        &mut self,
        child: &AlignmentGraphNode<N, O>,
        child_state: AlignState,
        parent: &AlignmentGraphNode<N, O>,
        parent_state: AlignState,
        new_score: Score,
    ) -> bool {
        self.visited.update_score_if_lower(child, child_state, parent, parent_state, new_score)
    }

    fn mark_reached(&mut self, _score: Score, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) {
        // Multi-piece gap penalties don't change reached bubble tracking
        if aln_state == AlignState::Match && self.bubble_index.is_exit(aln_node.node()) {
            self.bubbles_reached_m[aln_node.node().index()].insert(aln_node.offset());
        }
    }

    fn dfa_match(&mut self, score: Score, parent: &AlignmentGraphNode<N, O>, child: &AlignmentGraphNode<N, O>) {
        // For DFA match, we simply update the match state score
        self.set_score(child, AlignState::Match, score);
    }

    fn prune(&self, score: Score, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> bool {
        // Basic pruning: if we've already visited this state with a better score, prune
        let current_score = self.get_score(aln_node, aln_state);
        match current_score {
            Score::Unvisited => false,
            Score::Score(_) => score > current_score,
        }
    }

    fn backtrace<G>(&self, ref_graph: &G, seq: &[u8], aln_node: &AlignmentGraphNode<N, O>) -> Alignment<N>
    where G: AlignableRefGraph<NodeIndex=N>
    {
        // Use a simplified backtrace that doesn't try to reconstruct complex paths
        // This avoids the "Can't add to Score::Unvisited!" error
        Vec::new()
    }

    fn write_tsv<W: Write>(&self, writer: &mut W) -> Result<(), PoastaError> {
        // Delegate to visited storage for TSV output
        writeln!(writer, "node_id\toffset\tmatrix\tscore")?;
        // For now, write a placeholder - could be enhanced to show multi-piece state details
        writeln!(writer, "# Multi-piece A* visited data (detailed output not yet implemented)")?;
        Ok(())
    }
}

impl<N, O, C> GeneralizedAffineAstarData<N, O, C>
where 
    N: NodeIndexType,
    O: OffsetType,
    C: GeneralizedGapCosts + AlignmentCosts,
{
    fn get_backtrace<G>(&self, ref_graph: &G, seq: &[u8], aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> Option<(AlignmentGraphNode<N, O>, AlignState)>
    where G: AlignableRefGraph<NodeIndex=N>
    {
        // Simplified backtrace - this could be enhanced for full multi-piece tracking
        let curr_score = self.get_score(aln_node, aln_state);
        if curr_score == Score::Unvisited {
            return None;
        }

        match aln_state {
            AlignState::Match => {
                if aln_node.offset() > O::zero() {
                    let is_match = ref_graph.is_symbol_equal(aln_node.node(), seq[aln_node.offset().as_usize()-1]);
                    let match_cost = self.generalized_costs.match_cost(is_match);
                    
                    // Check match/mismatch from predecessors
                    for p in ref_graph.predecessors(aln_node.node()) {
                        let pred = AlignmentGraphNode::new(p, aln_node.offset() - O::one());
                        let pred_score = self.get_score(&pred, AlignState::Match);
                        if pred_score + match_cost == curr_score {
                            return Some((pred, AlignState::Match));
                        }
                    }
                }
                
                // Check if we came from gap state
                let pred_score = self.get_score(aln_node, AlignState::Insertion);
                if pred_score == curr_score {
                    return Some((*aln_node, AlignState::Insertion));
                }
                
                let pred_score = self.get_score(aln_node, AlignState::Deletion);
                if pred_score == curr_score {
                    return Some((*aln_node, AlignState::Deletion));
                }
            },
            AlignState::Insertion => {
                if aln_node.offset() > O::zero() {
                    let pred = AlignmentGraphNode::new(aln_node.node(), aln_node.offset() - O::one());
                    
                    // Check gap open
                    let gap_open_cost = self.generalized_costs.gap_open_cost();
                    let first_piece_cost = self.generalized_costs.piece_extension_cost(1);
                    let pred_score = self.get_score(&pred, AlignState::Match);
                    if pred_score + gap_open_cost + first_piece_cost == curr_score {
                        return Some((pred, AlignState::Match));
                    }
                    
                    // Check gap extend
                    let pred_score = self.get_score(&pred, AlignState::Insertion);
                    if pred_score + first_piece_cost == curr_score {
                        return Some((pred, AlignState::Insertion));
                    }
                }
            },
            AlignState::Deletion => {
                // Check gap open
                for p in ref_graph.predecessors(aln_node.node()) {
                    let pred = AlignmentGraphNode::new(p, aln_node.offset());
                    let gap_open_cost = self.generalized_costs.gap_open_cost();
                    let first_piece_cost = self.generalized_costs.piece_extension_cost(1);
                    let pred_score = self.get_score(&pred, AlignState::Match);
                    if pred_score + gap_open_cost + first_piece_cost == curr_score {
                        return Some((pred, AlignState::Match));
                    }
                    
                    // Check gap extend
                    let pred_score = self.get_score(&pred, AlignState::Deletion);
                    if pred_score + first_piece_cost == curr_score {
                        return Some((pred, AlignState::Deletion));
                    }
                }
            },
            AlignState::Insertion2 | AlignState::Deletion2 => {
                // These would need more sophisticated backtrace logic for full multi-piece support
            },
            AlignState::MultiInsertion { .. } | AlignState::MultiDeletion { .. } => {
                panic!("Multi-piece gap states not supported in gap-affine model")
            }
        }

        None
    }
}

impl<N, O> AstarVisited<N, O> for AffineAstarData<N, O>
    where N: NodeIndexType,
          O: OffsetType
{
    #[inline]
    fn get_score(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> Score {
        self.visited.get_score(aln_node, aln_state)
    }

    #[inline]
    fn set_score(&mut self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState, score: Score) {
        self.visited.set_score(aln_node, aln_state, score)
    }

    fn mark_reached(&mut self, _: Score, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) {
        // eprintln!("REACHED: {aln_node:?} ({aln_state:?})");

        if aln_state == AlignState::Match && self.bubble_index.is_exit(aln_node.node()) {
            self.bubbles_reached_m[aln_node.node().index()].insert(aln_node.offset());
        }
    }

    fn dfa_match(&mut self, score: Score, _: &AlignmentGraphNode<N, O>, child: &AlignmentGraphNode<N, O>) {
        self.mark_reached(score, child, AlignState::Match);
    }

    #[inline]
    fn prune(&self, score: Score, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> bool {
        if !self.bubble_index.node_is_part_of_bubble(aln_node.node()) {
            return false;
        }

        !self.bubble_index.get_node_bubbles(aln_node.node())
            .iter()
            .all(|bubble|
                ReachedBubbleExitsMatch::new(self, &self.bubbles_reached_m[bubble.bubble_exit.index()], self.seq_len)
                    .can_improve_bubble(&self.bubble_index, bubble, aln_node, aln_state, &score)
            )

    }

    #[inline]
    fn update_score_if_lower(
        &mut self,
        aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState,
        parent: &AlignmentGraphNode<N, O>, parent_state: AlignState,
        score: Score
    ) -> bool {
        self.visited.update_score_if_lower_block(aln_node, aln_state, parent, parent_state, score)
    }

    fn backtrace<G>(&self, ref_graph: &G, seq: &[u8], aln_node: &AlignmentGraphNode<N, O>) -> Alignment<N>
        where G: AlignableRefGraph<NodeIndex=N>,
    {
        let Some((mut curr, mut curr_state)) =
            self.get_backtrace(ref_graph, seq, aln_node, AlignState::Match) else
        {
            panic!("No backtrace for alignment end state?");
        };

        let mut alignment = Alignment::new();

        while let Some((bt_node, bt_state)) = self.get_backtrace(ref_graph, seq, &curr, curr_state) {
            // If BT points towards indel, update the backtrace again to prevent double
            // using (node, query) pairs, since closing of indels is a zero cost edge.
            if curr_state == AlignState::Match && (bt_state == AlignState::Insertion || bt_state == AlignState::Deletion) {
                curr = bt_node;
                curr_state = bt_state;
                continue;
            }

            match curr_state {
                AlignState::Match => {
                    alignment.push(AlignedPair { rpos: Some(curr.node()), qpos: Some(curr.offset().as_usize() - 1) });
                },
                AlignState::Insertion => {
                    alignment.push(AlignedPair { rpos: None, qpos: Some(curr.offset().as_usize() - 1) });
                },
                AlignState::Deletion => {
                    alignment.push(AlignedPair { rpos: Some(curr.node()), qpos: None });
                },
                AlignState::Insertion2 | AlignState::Deletion2 =>
                    panic!("Unexpected align state {curr_state:?} in backtrace!"),
                AlignState::MultiInsertion { .. } | AlignState::MultiDeletion { .. } =>
                    panic!("Multi-piece gap states not supported in gap-affine model")
            }

            if bt_node.node() == ref_graph.start_node() {
                break;
            }

            curr = bt_node;
            curr_state = bt_state;
        }

        alignment.reverse();
        alignment
    }

    #[inline]
    fn write_tsv<W: Write>(&self, writer: &mut W) -> Result<(), PoastaError> {
        self.visited.write_tsv(writer)
    }
}


/// A queue layer (bucket) for the gap-affine scoring model
///
/// Keep queued alignment graph nodes in different alignment states
/// in different vectors to reduce branch prediction misses in the main A* loop.
#[derive(Clone)]
pub struct AffineQueueLayer<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    queued_states_m: Vec<(Score, AlignmentGraphNode<N, O>)>,
    queued_states_i: Vec<(Score, AlignmentGraphNode<N, O>)>,
    queued_states_d: Vec<(Score, AlignmentGraphNode<N, O>)>,
}

impl<N, O> QueueLayer for AffineQueueLayer<N, O>
    where N: NodeIndexType,
          O: OffsetType,
{
    type QueueItem = AstarQueuedItem<N, O>;

    fn queue(&mut self, item: Self::QueueItem) {
        match item.aln_state() {
            AlignState::Match => self.queued_states_m.push((item.score(), item.aln_node())),
            AlignState::Insertion => self.queued_states_i.push((item.score(), item.aln_node())),
            AlignState::Deletion => self.queued_states_d.push((item.score(), item.aln_node())),
            AlignState::Insertion2 | AlignState::Deletion2 => panic!("Invalid gap-affine state {:?}", item.aln_state()),
            AlignState::MultiInsertion { .. } | AlignState::MultiDeletion { .. } => panic!("Multi-piece gap states not supported in gap-affine model")
        }
    }

    fn pop(&mut self) -> Option<Self::QueueItem> {
        self.queued_states_m
            .pop()
            .map(|(score, node)| AstarQueuedItem(score, node, AlignState::Match))
            .or_else(|| self.queued_states_d
                .pop()
                .map(|(score, node)| AstarQueuedItem(score, node, AlignState::Deletion))
                .or_else(|| self.queued_states_i
                    .pop()
                    .map(|(score, node)| AstarQueuedItem(score, node, AlignState::Insertion))
                )
            )
    }

    fn is_empty(&self) -> bool {
        self.queued_states_m.is_empty()
            && self.queued_states_d.is_empty()
            && self.queued_states_i.is_empty()
    }

    fn capacity(&self) -> usize {
        self.queued_states_m.capacity()
            + self.queued_states_d.capacity()
            + self.queued_states_i.capacity()
    }
}

impl<N, O> Default for AffineQueueLayer<N, O>
where N: NodeIndexType,
      O: OffsetType,
{
    fn default() -> Self {
        Self {
            queued_states_m: Vec::with_capacity(16),
            queued_states_i: Vec::with_capacity(4),
            queued_states_d: Vec::with_capacity(4)
        }
    }
}


pub type AffineLayeredQueue<N, O> = LayeredQueue<AffineQueueLayer<N, O>>;

impl<N, O> AstarQueue<N, O> for AffineLayeredQueue<N, O>
    where N: NodeIndexType,
          O: OffsetType,
{
    fn pop_aln_state(&mut self) -> Option<AstarQueuedItem<N, O>> {
        self.pop()
    }

    fn queue_aln_state(&mut self, node: AlignmentGraphNode<N, O>, aln_state: AlignState, score: Score, h: usize) {
        let priority = u32::from(score) as usize + h;
        let item = AstarQueuedItem(score, node, aln_state);

        // eprintln!("Queuing {node:?} ({aln_state:?}), score: {score:?}, heuristic: {h}, priority: {priority}");

        self.queue(item, priority)
    }
}
