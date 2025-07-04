use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::fmt::Write;
use std::marker::PhantomData;
use std::sync::Arc;
use rustc_hash::FxHashMap;

use crate::graphs::{AlignableRefGraph, NodeIndexType};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{GeneralizedGapCosts, GetAlignmentCosts, Score};
use crate::aligner::aln_graph::{AlignmentGraph, AlignmentGraphNode, GeneralizedAlignState, BasicAlignState, AlignState};
use crate::aligner::astar::{AstarQueue, AstarQueuedItem, AstarVisited};
use crate::aligner::queue::{LayeredQueue, QueueLayer};
use crate::bubbles::index::BubbleIndex;
use crate::bubbles::reached::ReachedBubbleExitsMatch;
use crate::errors::PoastaError;
use crate::aligner::{AlignedPair, Alignment};

/// Unified A* data structure for generalized gap cost models
/// Works with both traditional affine and multi-piece gap penalties
pub struct GeneralizedGapAstarData<N, O, C>
where 
    N: NodeIndexType,
    O: OffsetType,
    C: GeneralizedGapCosts,
{
    costs: C,
    seq_len: usize,
    bubble_index: Arc<BubbleIndex<N>>,
    visited: BlockedVisitedStorageGeneralized<N, O>,
    bubbles_reached_m: Vec<BTreeSet<O>>,
}

impl<N, O, C> GeneralizedGapAstarData<N, O, C>
where 
    N: NodeIndexType,
    O: OffsetType,
    C: GeneralizedGapCosts,
{
    pub fn new<G>(costs: C, ref_graph: &G, seq: &[u8], bubble_index: Arc<BubbleIndex<G::NodeIndex>>) -> Self
    where G: AlignableRefGraph<NodeIndex=N>,
    {
        Self {
            costs,
            seq_len: seq.len(),
            bubble_index,
            visited: BlockedVisitedStorageGeneralized::new(ref_graph),
            bubbles_reached_m: vec![BTreeSet::new(); ref_graph.node_count_with_start_and_end()],
        }
    }
}

impl<N, O, C> AstarVisited<N, O> for GeneralizedGapAstarData<N, O, C>
where 
    N: NodeIndexType,
    O: OffsetType,
    C: GeneralizedGapCosts,
{
    fn get_score(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> Score {
        // Convert legacy AlignState to GeneralizedAlignState
        let generalized_state = self.convert_legacy_state(aln_state);
        self.visited.get_score(aln_node, generalized_state)
    }

    fn set_score(&mut self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState, score: Score) {
        let generalized_state = self.convert_legacy_state(aln_state);
        self.visited.set_score(aln_node, generalized_state, score);
    }

    fn update_score_if_lower(
        &mut self,
        aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState,
        parent: &AlignmentGraphNode<N, O>, parent_state: AlignState,
        score: Score
    ) -> bool {
        let generalized_state = self.convert_legacy_state(aln_state);
        self.visited.update_score_if_lower(aln_node, generalized_state, score)
    }

    fn mark_reached(&mut self, _score: Score, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) {
        if aln_state == AlignState::Match && self.bubble_index.is_exit(aln_node.node()) {
            self.bubbles_reached_m[aln_node.node().index()].insert(aln_node.offset());
        }
    }

    fn dfa_match(&mut self, score: Score, _parent: &AlignmentGraphNode<N, O>, child: &AlignmentGraphNode<N, O>) {
        self.mark_reached(score, child, AlignState::Match);
    }

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

    fn backtrace<G>(&self, ref_graph: &G, seq: &[u8], aln_node: &AlignmentGraphNode<N, O>) -> Alignment<N>
    where G: AlignableRefGraph<NodeIndex=N>,
    {
        // For now, simplified backtrace - can be enhanced later
        Alignment::new()
    }

    fn write_tsv<W: Write>(&self, writer: &mut W) -> Result<(), PoastaError> {
        self.visited.write_tsv(writer)
    }
}

impl<N, O, C> GeneralizedGapAstarData<N, O, C>
where 
    N: NodeIndexType,
    O: OffsetType,
    C: GeneralizedGapCosts,
{
    /// Convert legacy AlignState to GeneralizedAlignState for backward compatibility
    fn convert_legacy_state(&self, state: AlignState) -> GeneralizedAlignState {
        match state {
            AlignState::Match => GeneralizedAlignState::match_state(),
            AlignState::Insertion => GeneralizedAlignState::new_gap_state(BasicAlignState::Insertion),
            AlignState::Deletion => GeneralizedAlignState::new_gap_state(BasicAlignState::Deletion),
            AlignState::Insertion2 => GeneralizedAlignState::gap_state(BasicAlignState::Insertion, 2, 1),
            AlignState::Deletion2 => GeneralizedAlignState::gap_state(BasicAlignState::Deletion, 2, 1),
        }
    }
}

/// Enhanced visited cell structure for generalized gap tracking
#[derive(Default, Copy, Clone)]
struct VisitedCellGeneralized {
    visited_m: Score,
    /// Gap states: indexed by [piece_index-1][piece_position-1]
    /// We support up to 4 pieces with up to 16 positions each for efficiency
    visited_i: [[Score; 16]; 4],
    visited_d: [[Score; 16]; 4],
}

/// Blocked storage for generalized gap alignment scores
struct BlockedVisitedStorageGeneralized<N, O, const B: usize = 8>
where 
    N: NodeIndexType,
    O: OffsetType,
{
    node_blocks: Vec<FxHashMap<O, [[VisitedCellGeneralized; B]; B]>>,
    node_ranks: Vec<usize>,
    dummy: PhantomData<N>,
}

impl<N, O, const B: usize> BlockedVisitedStorageGeneralized<N, O, B>
where 
    N: NodeIndexType,
    O: OffsetType,
{
    pub fn new<G: AlignableRefGraph>(ref_graph: &G) -> Self {
        if B & (B-1) != 0 {
            panic!("Block size B should be a power of 2!");
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
    pub fn get_score(&self, aln_node: &AlignmentGraphNode<N, O>, state: GeneralizedAlignState) -> Score {
        let (node_block, offset_block, within_block_node, within_block_qry)
            = self.calc_block_ix(aln_node);

        self.node_blocks[node_block].get(&offset_block)
            .map(|v| {
                let cell_data = &v[within_block_node][within_block_qry];
                match state.basic_state {
                    BasicAlignState::Match => cell_data.visited_m,
                    BasicAlignState::Insertion => {
                        if state.piece_index == 0 {
                            Score::Unvisited
                        } else {
                            let piece_idx = (state.piece_index - 1) as usize;
                            let pos_idx = (state.piece_position - 1) as usize;
                            if piece_idx < 4 && pos_idx < 16 {
                                cell_data.visited_i[piece_idx][pos_idx]
                            } else {
                                Score::Unvisited
                            }
                        }
                    }
                    BasicAlignState::Deletion => {
                        if state.piece_index == 0 {
                            Score::Unvisited
                        } else {
                            let piece_idx = (state.piece_index - 1) as usize;
                            let pos_idx = (state.piece_position - 1) as usize;
                            if piece_idx < 4 && pos_idx < 16 {
                                cell_data.visited_d[piece_idx][pos_idx]
                            } else {
                                Score::Unvisited
                            }
                        }
                    }
                }
            })
            .unwrap_or(Score::Unvisited)
    }

    #[inline]
    pub fn set_score(&mut self, aln_node: &AlignmentGraphNode<N, O>, state: GeneralizedAlignState, score: Score) {
        let (node_block, offset_block, within_block_node, within_block_qry)
            = self.calc_block_ix(aln_node);

        let v = self.node_blocks[node_block].entry(offset_block)
            .or_insert_with(|| [[VisitedCellGeneralized::default(); B]; B]);

        let cell_data = &mut v[within_block_node][within_block_qry];
        match state.basic_state {
            BasicAlignState::Match => cell_data.visited_m = score,
            BasicAlignState::Insertion => {
                if state.piece_index > 0 {
                    let piece_idx = (state.piece_index - 1) as usize;
                    let pos_idx = (state.piece_position - 1) as usize;
                    if piece_idx < 4 && pos_idx < 16 {
                        cell_data.visited_i[piece_idx][pos_idx] = score;
                    }
                }
            }
            BasicAlignState::Deletion => {
                if state.piece_index > 0 {
                    let piece_idx = (state.piece_index - 1) as usize;
                    let pos_idx = (state.piece_position - 1) as usize;
                    if piece_idx < 4 && pos_idx < 16 {
                        cell_data.visited_d[piece_idx][pos_idx] = score;
                    }
                }
            }
        }
    }

    #[inline]
    pub fn update_score_if_lower(
        &mut self,
        aln_node: &AlignmentGraphNode<N, O>, state: GeneralizedAlignState,
        score: Score
    ) -> bool {
        let current_score = self.get_score(aln_node, state);
        match score.cmp(&current_score) {
            Ordering::Less => {
                self.set_score(aln_node, state, score);
                true
            },
            Ordering::Equal | Ordering::Greater => false,
        }
    }

    pub fn write_tsv<W: Write>(&self, writer: &mut W) -> Result<(), PoastaError> {
        writeln!(writer, "node_id\toffset\tmatrix\tpiece\tposition\tscore")?;
        // Simplified implementation
        Ok(())
    }
}