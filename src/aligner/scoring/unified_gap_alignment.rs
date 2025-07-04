
use crate::graphs::{AlignableRefGraph, NodeIndexType};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, AlignmentType, GeneralizedGapCosts, Score};
use crate::aligner::aln_graph::{AlignmentGraph, AlignmentGraphNode, AlignState};
use crate::aligner::astar::AstarVisited;
use crate::aligner::scoring::gap_affine::{AffineLayeredQueue};

/// Unified alignment graph that implements generalized multi-piece gap penalties
/// This single implementation handles both traditional affine (1-piece) and multi-piece gaps
pub struct UnifiedGapAlignmentGraph<C>
where C: GeneralizedGapCosts,
{
    costs: UnifiedAlignmentCosts<C>,
    aln_type: AlignmentType,
}

impl<C> UnifiedGapAlignmentGraph<C>
where C: GeneralizedGapCosts,
{
    pub fn new(costs: C, aln_type: AlignmentType) -> Self {
        Self { 
            costs: UnifiedAlignmentCosts::new(costs),
            aln_type 
        }
    }
}

impl<C> AlignmentGraph for UnifiedGapAlignmentGraph<C>
where C: GeneralizedGapCosts,
{
    type CostModel = UnifiedAlignmentCosts<C>;

    fn get_costs(&self) -> &Self::CostModel {
        &self.costs
    }

    fn initial_states<G, O>(&self, ref_graph: &G) -> Vec<AlignmentGraphNode<G::NodeIndex, O>>
    where 
        G: AlignableRefGraph,
        O: OffsetType,
    {
        match self.aln_type {
            AlignmentType::Global => vec![AlignmentGraphNode::new(ref_graph.start_node(), O::zero())],
            AlignmentType::EndsFree { .. } => {
                todo!("Ends-free alignment not yet implemented for unified gap model");
            }
        }
    }

    fn is_end<G, O>(&self, ref_graph: &G, seq: &[u8], node: &AlignmentGraphNode<G::NodeIndex, O>, aln_state: AlignState) -> bool
    where 
        G: AlignableRefGraph,
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
    where 
        V: AstarVisited<G::NodeIndex, O>,
        G: AlignableRefGraph,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<G::NodeIndex, O>, AlignState)
    {
        match state {
            AlignState::Match => {
                self.expand_from_match(visited_data, ref_graph, seq, score, node, &mut f);
            },
            AlignState::Insertion => {
                self.expand_insertion_unified(visited_data, ref_graph, seq, score, node, 1, 1, &mut f);
            },
            AlignState::Deletion => {
                self.expand_deletion_unified(visited_data, ref_graph, seq, score, node, 1, 1, &mut f);
            },
            AlignState::Insertion2 => {
                // Legacy support: treat as piece 2, position 1
                self.expand_insertion_unified(visited_data, ref_graph, seq, score, node, 2, 1, &mut f);
            },
            AlignState::Deletion2 => {
                // Legacy support: treat as piece 2, position 1
                self.expand_deletion_unified(visited_data, ref_graph, seq, score, node, 2, 1, &mut f);
            },
            AlignState::MultiInsertion { piece, position } => {
                self.expand_insertion_unified(visited_data, ref_graph, seq, score, node, piece as usize, position as usize, &mut f);
            },
            AlignState::MultiDeletion { piece, position } => {
                self.expand_deletion_unified(visited_data, ref_graph, seq, score, node, piece as usize, position as usize, &mut f);
            },
        }
    }

    fn expand_ref_graph_end<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        score: Score,
        mut f: F
    )
    where 
        V: AstarVisited<N, O>,
        N: NodeIndexType,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState)
    {
        // At ref graph end, open insertion
        let new_node_ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        let gap_open_cost = self.costs.costs.gap_open_cost();
        let first_piece_cost = self.costs.costs.piece_extension_cost(1);
        let new_score_ins = score + gap_open_cost + first_piece_cost;
        let new_ins_state = AlignState::multi_insertion(1, 1);
        
        if visited_data
            .update_score_if_lower(&new_node_ins, new_ins_state, parent, AlignState::Match, new_score_ins)
        {
            f(gap_open_cost + first_piece_cost, new_node_ins, new_ins_state);
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
    where 
        V: AstarVisited<N, O>,
        N: NodeIndexType,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState)
    {
        // At query end, open deletion
        let new_node_del = AlignmentGraphNode::new(child, parent.offset());
        let gap_open_cost = self.costs.costs.gap_open_cost();
        let first_piece_cost = self.costs.costs.piece_extension_cost(1);
        let new_score_del = score + gap_open_cost + first_piece_cost;
        let new_del_state = AlignState::multi_deletion(1, 1);
        
        if visited_data
            .update_score_if_lower(&new_node_del, new_del_state, parent, AlignState::Match, new_score_del)
        {
            f(gap_open_cost + first_piece_cost, new_node_del, new_del_state);
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
    where 
        V: AstarVisited<N, O>,
        N: NodeIndexType,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState)
    {
        let mismatch_cost = self.costs.costs.match_cost(false);
        let new_score_mis = score + mismatch_cost;
        
        if visited_data
            .update_score_if_lower(child, AlignState::Match, parent, AlignState::Match, new_score_mis)
        {
            f(mismatch_cost, *child, AlignState::Match);
        }

        // Also queue indel states from parent
        let gap_open_cost = self.costs.costs.gap_open_cost();
        let first_piece_cost = self.costs.costs.piece_extension_cost(1);
        
        // Open insertion
        let new_node_ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        let new_score_ins = score + gap_open_cost + first_piece_cost;
        let new_ins_state = AlignState::multi_insertion(1, 1);

        if visited_data
            .update_score_if_lower(&new_node_ins, new_ins_state, parent, AlignState::Match, new_score_ins)
        {
            f(gap_open_cost + first_piece_cost, new_node_ins, new_ins_state);
        }

        // Open deletion
        let new_node_del = AlignmentGraphNode::new(child.node(), parent.offset());
        let new_score_del = score + gap_open_cost + first_piece_cost;
        let new_del_state = AlignState::multi_deletion(1, 1);
        
        if visited_data
            .update_score_if_lower(&new_node_del, new_del_state, parent, AlignState::Match, new_score_del)
        {
            f(gap_open_cost + first_piece_cost, new_node_del, new_del_state);
        }
    }
}

impl<C> UnifiedGapAlignmentGraph<C>
where C: GeneralizedGapCosts,
{
    /// Expand from match state - unified for all gap models
    fn expand_from_match<V, G, O, F>(
        &self,
        visited_data: &mut V,
        ref_graph: &G,
        seq: &[u8],
        score: Score,
        node: &AlignmentGraphNode<G::NodeIndex, O>,
        f: &mut F
    )
    where 
        V: AstarVisited<G::NodeIndex, O>,
        G: AlignableRefGraph,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<G::NodeIndex, O>, AlignState)
    {
        let child_offset = node.offset().increase_one();
        
        // Expand to match/mismatch states
        for ref_succ in ref_graph.successors(node.node()) {
            let new_node_mis = AlignmentGraphNode::new(ref_succ, child_offset);

            // Calculate match/mismatch cost
            let is_match = ref_graph.is_symbol_equal(ref_succ, seq[child_offset.as_usize()-1]);
            let score_delta = self.costs.costs.match_cost(is_match);
            let new_score_mis = score + score_delta;

            if visited_data
                .update_score_if_lower(&new_node_mis, AlignState::Match, node, AlignState::Match, new_score_mis)
            {
                f(score_delta, new_node_mis, AlignState::Match);
            }

            // Open deletion (start new gap in first piece, first position)
            let gap_open_cost = self.costs.costs.gap_open_cost();
            let first_piece_cost = self.costs.costs.piece_extension_cost(1);
            let score_delta = gap_open_cost + first_piece_cost;
            let new_node_del = AlignmentGraphNode::new(ref_succ, node.offset());
            let new_score_del = score + score_delta;
            let new_del_state = AlignState::multi_deletion(1, 1);
            
            if visited_data
                .update_score_if_lower(&new_node_del, new_del_state, node, AlignState::Match, new_score_del)
            {
                f(score_delta, new_node_del, new_del_state);
            }
        }

        // Open insertion (start new gap in first piece, first position)
        let gap_open_cost = self.costs.costs.gap_open_cost();
        let first_piece_cost = self.costs.costs.piece_extension_cost(1);
        let new_node_ins = AlignmentGraphNode::new(node.node(), child_offset);
        let new_score_ins = score + gap_open_cost + first_piece_cost;
        let new_ins_state = AlignState::multi_insertion(1, 1);
        
        if child_offset.as_usize() <= seq.len() && visited_data
            .update_score_if_lower(&new_node_ins, new_ins_state, node, AlignState::Match, new_score_ins)
        {
            f(gap_open_cost + first_piece_cost, new_node_ins, new_ins_state);
        }
    }

    /// Unified insertion expansion that handles multi-piece logic
    fn expand_insertion_unified<V, G, O, F>(
        &self,
        visited_data: &mut V,
        ref_graph: &G,
        seq: &[u8],
        score: Score,
        node: &AlignmentGraphNode<G::NodeIndex, O>,
        current_piece: usize,
        current_position: usize,
        f: &mut F
    )
    where 
        V: AstarVisited<G::NodeIndex, O>,
        G: AlignableRefGraph,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<G::NodeIndex, O>, AlignState)
    {
        // Create the current state
        let current_state = AlignState::multi_insertion(current_piece as u8, current_position as u8);
        
        // Insertion -> Match: zero cost transition
        if visited_data.update_score_if_lower(node, AlignState::Match, node, current_state, score) {
            f(0, *node, AlignState::Match);
        }

        // Extend insertion within current piece or transition to next piece
        let new_node_ins = AlignmentGraphNode::new(node.node(), node.offset().increase_one());
        
        if node.offset().as_usize() < seq.len() {
            let (should_transition, next_piece) = self.costs.costs.should_transition_piece(current_piece, current_position);
            
            let (extension_cost, next_state) = if should_transition {
                // Transition to next piece
                (self.costs.costs.piece_extension_cost(next_piece), 
                 AlignState::multi_insertion(next_piece as u8, 1))
            } else {
                // Continue in current piece
                (self.costs.costs.piece_extension_cost(current_piece), 
                 AlignState::multi_insertion(current_piece as u8, (current_position + 1) as u8))
            };
            
            let new_score_ins = score + extension_cost;

            if visited_data
                .update_score_if_lower(&new_node_ins, next_state, node, current_state, new_score_ins)
            {
                f(extension_cost, new_node_ins, next_state);
            }
        }
    }

    /// Unified deletion expansion that handles multi-piece logic
    fn expand_deletion_unified<V, G, O, F>(
        &self,
        visited_data: &mut V,
        ref_graph: &G,
        _seq: &[u8],
        score: Score,
        node: &AlignmentGraphNode<G::NodeIndex, O>,
        current_piece: usize,
        current_position: usize,
        f: &mut F
    )
    where 
        V: AstarVisited<G::NodeIndex, O>,
        G: AlignableRefGraph,
        O: OffsetType,
        F: FnMut(u8, AlignmentGraphNode<G::NodeIndex, O>, AlignState)
    {
        // Create the current state
        let current_state = AlignState::multi_deletion(current_piece as u8, current_position as u8);
        
        // Deletion -> Match: zero cost transition
        if visited_data.update_score_if_lower(node, AlignState::Match, node, current_state, score) {
            f(0, *node, AlignState::Match);
        }

        // Extend deletion within current piece or transition to next piece
        for ref_succ in ref_graph.successors(node.node()) {
            let (should_transition, next_piece) = self.costs.costs.should_transition_piece(current_piece, current_position);
            
            let (extension_cost, next_state) = if should_transition {
                // Transition to next piece
                (self.costs.costs.piece_extension_cost(next_piece), 
                 AlignState::multi_deletion(next_piece as u8, 1))
            } else {
                // Continue in current piece
                (self.costs.costs.piece_extension_cost(current_piece), 
                 AlignState::multi_deletion(current_piece as u8, (current_position + 1) as u8))
            };

            let new_node_del = AlignmentGraphNode::new(ref_succ, node.offset());
            let new_score_del = score + extension_cost;
            
            if visited_data
                .update_score_if_lower(&new_node_del, next_state, node, current_state, new_score_del)
            {
                f(extension_cost, new_node_del, next_state);
            }
        }
    }
}

/// Unified alignment costs implementation that works with any GeneralizedGapCosts
#[derive(Copy, Clone)]
pub struct UnifiedAlignmentCosts<C>
where C: GeneralizedGapCosts,
{
    costs: C,
}

impl<C> UnifiedAlignmentCosts<C>
where C: GeneralizedGapCosts,
{
    pub fn new(costs: C) -> Self {
        Self { costs }
    }
}

impl<C> AlignmentCosts for UnifiedAlignmentCosts<C>
where C: GeneralizedGapCosts,
{
    type AlignmentGraphType = UnifiedGapAlignmentGraph<C>;
    type QueueType<N, O> = AffineLayeredQueue<N, O>
        where N: NodeIndexType,
              O: OffsetType;

    fn new_alignment_graph(&self, aln_type: AlignmentType) -> Self::AlignmentGraphType {
        UnifiedGapAlignmentGraph::new(self.costs, aln_type)
    }

    fn mismatch(&self) -> u8 {
        self.costs.match_cost(false)
    }
    
    fn gap_open(&self) -> u8 {
        self.costs.gap_open_cost()
    }
    
    fn gap_extend(&self) -> u8 {
        self.costs.piece_extension_cost(1)
    }
    
    fn gap_open2(&self) -> u8 {
        0 // Not applicable for unified model
    }
    
    fn gap_extend2(&self) -> u8 {
        if self.costs.num_pieces() > 1 {
            self.costs.piece_extension_cost(2)
        } else {
            0
        }
    }

    fn gap_cost(&self, current_state: AlignState, length: usize) -> usize {
        match current_state {
            AlignState::Insertion | AlignState::Deletion => 0, // Already in gap
            AlignState::Match => self.costs.total_gap_cost(length),
            AlignState::Insertion2 | AlignState::Deletion2 => 0, // Already in gap
            AlignState::MultiInsertion { .. } | AlignState::MultiDeletion { .. } => 0, // Already in multi-piece gap
        }
    }
}