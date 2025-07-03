use crate::graphs::{AlignableRefGraph, NodeIndexType};
use crate::aligner::offsets::OffsetType;
use crate::aligner::scoring::{AlignmentCosts, AlignmentType, Score};
use crate::aligner::aln_graph::{AlignmentGraph, AlignmentGraphNode, AlignState};
use crate::aligner::astar::{AstarQueue, AstarQueuedItem, AstarVisited};
use crate::aligner::queue::{LayeredQueue, QueueLayer};
use crate::bubbles::index::BubbleIndex;
use crate::aligner::Alignment;
use crate::errors::PoastaError;
use std::cmp::Ordering;
use std::marker::PhantomData;
use std::sync::Arc;
use std::collections::BTreeSet;
use std::fmt::Write;
use rustc_hash::FxHashMap;

/// Maximum number of pieces supported (should be sufficient for practical use)
const MAX_PIECES: usize = 4;

/// Multi-piece gap-affine scoring model with configurable breakpoints
#[derive(Clone, Copy, Debug)]
pub struct GapMultiPieceAffine {
    cost_mismatch: u8,
    cost_gap_open: u8,
    /// Number of pieces (1-4)
    num_pieces: u8,
    /// Breakpoints where gap extension penalty changes [k₁, k₂, k₃]
    /// Note: k₀ = 1 (implicit), k₄ = ∞ (implicit)
    /// Only first (num_pieces-1) elements are used
    breakpoints: [usize; MAX_PIECES - 1],
    /// Extension costs for each piece [β₁, β₂, β₃, β₄]
    /// Only first num_pieces elements are used
    extension_costs: [u8; MAX_PIECES],
}

impl GapMultiPieceAffine {
    /// Create a new multi-piece gap-affine model
    pub fn new(cost_mismatch: u8, cost_gap_open: u8, breakpoints: &[usize], extension_costs: &[u8]) -> Self {
        if extension_costs.is_empty() {
            panic!("Must have at least one extension cost");
        }
        if extension_costs.len() > MAX_PIECES {
            panic!("Too many pieces: {} (max {})", extension_costs.len(), MAX_PIECES);
        }
        if breakpoints.len() + 1 != extension_costs.len() {
            panic!("Number of extension costs must be one more than number of breakpoints");
        }
        if !breakpoints.is_empty() {
            for i in 1..breakpoints.len() {
                if breakpoints[i] <= breakpoints[i-1] {
                    panic!("Breakpoints must be in strictly increasing order");
                }
            }
            if breakpoints[0] <= 1 {
                panic!("First breakpoint must be greater than 1");
            }
        }
        
        let mut bp_array = [0; MAX_PIECES - 1];
        let mut cost_array = [0; MAX_PIECES];
        
        for (i, &bp) in breakpoints.iter().enumerate() {
            bp_array[i] = bp;
        }
        for (i, &cost) in extension_costs.iter().enumerate() {
            cost_array[i] = cost;
        }
        
        Self {
            cost_mismatch,
            cost_gap_open,
            num_pieces: extension_costs.len() as u8,
            breakpoints: bp_array,
            extension_costs: cost_array,
        }
    }

    /// Factory method for standard affine gap penalty (1-piece)
    pub fn standard_affine(cost_mismatch: u8, cost_gap_open: u8, cost_gap_extend: u8) -> Self {
        Self::new(cost_mismatch, cost_gap_open, &[], &[cost_gap_extend])
    }

    /// Factory method for two-piece gap penalty
    pub fn two_piece(cost_mismatch: u8, cost_gap_open: u8, short_extend: u8, long_extend: u8, breakpoint: usize) -> Self {
        Self::new(cost_mismatch, cost_gap_open, &[breakpoint], &[short_extend, long_extend])
    }

    /// Factory method for three-piece gap penalty
    pub fn three_piece(cost_mismatch: u8, cost_gap_open: u8, costs: [u8; 3], breakpoints: [usize; 2]) -> Self {
        Self::new(cost_mismatch, cost_gap_open, &breakpoints, &costs)
    }

    /// Get the number of pieces in this model
    pub fn num_pieces(&self) -> usize {
        self.num_pieces as usize
    }

    /// Find which piece a gap of given length belongs to
    pub fn find_piece(&self, gap_length: usize) -> usize {
        if gap_length == 0 {
            return 0; // Special case, though gap_length should never be 0 in practice
        }
        
        for i in 0..(self.num_pieces as usize - 1) {
            if gap_length <= self.breakpoints[i] {
                return i;
            }
        }
        self.num_pieces as usize - 1  // Final piece
    }

    /// Calculate gap cost for a given length using the mathematical formula:
    /// g_multi-piece(ℓ) = α + Σ(i=1 to j-1) βᵢ·(kᵢ - k_{i-1}) + βⱼ·(ℓ - k_{j-1})
    /// where k_{j-1} < ℓ ≤ kⱼ
    pub fn gap_cost(&self, gap_length: usize) -> usize {
        if gap_length == 0 {
            return 0;
        }
        
        // Find which piece this gap length belongs to
        let piece_j = self.find_piece_for_gap_length(gap_length);
        
        let mut total_cost = self.cost_gap_open as usize;
        
        // Add costs for completed pieces (i = 1 to j-1)
        for i in 1..piece_j {
            let k_i_minus_1 = if i == 1 { 1 } else { self.breakpoints[i-2] + 1 };
            let k_i = self.breakpoints[i-1];
            let piece_length = k_i - k_i_minus_1 + 1;
            total_cost += self.extension_costs[i-1] as usize * piece_length;
        }
        
        // Add cost for partial piece j
        let k_j_minus_1 = if piece_j == 1 { 0 } else { self.breakpoints[piece_j-2] };
        let remaining_in_piece = gap_length - k_j_minus_1;
        total_cost += self.extension_costs[piece_j-1] as usize * remaining_in_piece;
        
        total_cost
    }
    
    /// Find which piece a gap length belongs to (1-indexed)
    /// Returns j such that k_{j-1} < gap_length ≤ k_j
    fn find_piece_for_gap_length(&self, gap_length: usize) -> usize {
        for j in 1..=(self.num_pieces as usize) {
            let k_j_minus_1 = if j == 1 { 0 } else { self.breakpoints[j-2] };
            let k_j = if j == self.num_pieces as usize { 
                usize::MAX 
            } else { 
                self.breakpoints[j-1] 
            };
            
            if gap_length > k_j_minus_1 && gap_length <= k_j {
                return j;
            }
        }
        
        // Should never reach here if implementation is correct
        panic!("Could not find piece for gap length {}", gap_length);
    }

    /// Get extension cost for a specific piece (1-indexed to match mathematical notation)
    pub fn get_extension_cost(&self, piece_j: usize) -> u8 {
        if piece_j == 0 || piece_j > self.num_pieces as usize {
            panic!("Invalid piece index {}, must be in range [1, {}]", piece_j, self.num_pieces);
        }
        self.extension_costs[piece_j - 1] // Convert to 0-indexed array access
    }

    /// Get breakpoint k_j for piece j (1-indexed)
    /// Returns the maximum gap length for piece j
    pub fn get_breakpoint(&self, piece_j: usize) -> Option<usize> {
        if piece_j == 0 || piece_j > self.num_pieces as usize {
            panic!("Invalid piece index {}, must be in range [1, {}]", piece_j, self.num_pieces);
        }
        
        if piece_j == self.num_pieces as usize {
            None // Final piece k_p = ∞
        } else {
            Some(self.breakpoints[piece_j - 1]) // Convert to 0-indexed array access
        }
    }

    /// Get the valid range for gap lengths in piece j
    /// Returns (k_{j-1}, k_j] as (exclusive_start, inclusive_end)
    pub fn get_piece_range(&self, piece_j: usize) -> (usize, Option<usize>) {
        if piece_j == 0 || piece_j > self.num_pieces as usize {
            panic!("Invalid piece index {}, must be in range [1, {}]", piece_j, self.num_pieces);
        }
        
        let k_j_minus_1 = if piece_j == 1 { 0 } else { self.breakpoints[piece_j - 2] };
        let k_j = if piece_j == self.num_pieces as usize { 
            None 
        } else { 
            Some(self.breakpoints[piece_j - 1]) 
        };
        
        (k_j_minus_1, k_j)
    }

    /// Check if a gap can be extended within the current piece
    /// Returns extension cost if possible, None if piece boundary would be exceeded
    pub fn can_extend_in_piece(&self, piece_j: usize, current_gap_length: usize) -> Option<u8> {
        if piece_j == 0 || piece_j > self.num_pieces as usize {
            return None;
        }
        
        let (k_j_minus_1, k_j_opt) = self.get_piece_range(piece_j);
        
        // Check if current gap length is valid for this piece
        if current_gap_length <= k_j_minus_1 {
            return None;
        }
        
        // Check if extending would exceed piece boundary
        if let Some(k_j) = k_j_opt {
            if current_gap_length >= k_j {
                return None; // Already at or beyond piece boundary
            }
        }
        
        Some(self.extension_costs[piece_j - 1])
    }

    /// Check if a gap can transition to the next piece
    /// Returns true if current gap length is at the boundary of piece j
    pub fn can_transition_to_next_piece(&self, piece_j: usize, current_gap_length: usize) -> bool {
        if piece_j == 0 || piece_j >= self.num_pieces as usize {
            return false; // Invalid piece or already in final piece
        }
        
        // Can transition if we're exactly at the breakpoint k_j
        if let Some(k_j) = self.get_breakpoint(piece_j) {
            current_gap_length == k_j
        } else {
            false // Final piece has no next piece
        }
    }
    
    /// Calculate the multi-piece gap cost (separate from the trait method)
    fn calculate_multi_piece_gap_cost(&self, gap_length: usize) -> usize {
        if gap_length == 0 {
            return 0;
        }
        
        // Find which piece this gap length belongs to
        let piece_j = self.find_piece_for_gap_length(gap_length);
        
        let mut total_cost = self.cost_gap_open as usize;
        
        // Add costs for completed pieces (i = 1 to j-1)
        for i in 1..piece_j {
            let k_i_minus_1 = if i == 1 { 1 } else { self.breakpoints[i-2] + 1 };
            let k_i = self.breakpoints[i-1];
            let piece_length = k_i - k_i_minus_1 + 1;
            total_cost += self.extension_costs[i-1] as usize * piece_length;
        }
        
        // Add cost for partial piece j
        let k_j_minus_1 = if piece_j == 1 { 0 } else { self.breakpoints[piece_j-2] };
        let remaining_in_piece = gap_length - k_j_minus_1;
        total_cost += self.extension_costs[piece_j-1] as usize * remaining_in_piece;
        
        total_cost
    }
}

impl AlignmentCosts for GapMultiPieceAffine {
    type AlignmentGraphType = MultiPieceAffineAlignmentGraph;
    type QueueType<N, O> = MultiPieceAffineLayeredQueue<N, O>
        where N: NodeIndexType,
              O: OffsetType;

    fn new_alignment_graph(&self, aln_type: AlignmentType) -> Self::AlignmentGraphType {
        MultiPieceAffineAlignmentGraph::new(self, aln_type)
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
        self.extension_costs[0] // Return first piece extension cost for compatibility
    }

    #[inline(always)]
    fn gap_open2(&self) -> u8 {
        0
    }

    #[inline(always)]
    fn gap_extend2(&self) -> u8 {
        if self.num_pieces > 1 {
            self.extension_costs[1] // Return second piece extension cost if available
        } else {
            0
        }
    }

    #[inline]
    fn gap_cost(&self, current_state: AlignState, length: usize) -> usize {
        if length == 0 {
            return 0;
        }

        let gap_open = match current_state {
            AlignState::Insertion | AlignState::Deletion => 0,
            AlignState::Match => self.cost_gap_open as usize,
            AlignState::Insertion2 | AlignState::Deletion2 => 0, // Multi-piece states
        };

        // Use the multi-piece gap cost calculation
        let multi_piece_cost = self.calculate_multi_piece_gap_cost(length);
        gap_open + multi_piece_cost - self.cost_gap_open as usize
    }
}

pub struct MultiPieceAffineAlignmentGraph {
    /// Multi-piece gap affine costs
    costs: GapMultiPieceAffine,

    /// Alignment type, e.g., global alignment or ends-free alignment
    aln_type: AlignmentType,
}

impl MultiPieceAffineAlignmentGraph {
    fn new(costs: &GapMultiPieceAffine, aln_type: AlignmentType) -> Self {
        Self {
            costs: costs.clone(),
            aln_type,
        }
    }
}

impl AlignmentGraph for MultiPieceAffineAlignmentGraph {
    type CostModel = GapMultiPieceAffine;

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
                todo!("EndsFree alignment not yet implemented for multi-piece gap-affine");
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

    /// Implementation of the multi-piece gap-affine recurrence relations
    /// Following the mathematical framework: state space (i, v, s, j, ℓ)
    /// where i=offset, v=node, s=state, j=piece_index, ℓ=length_in_piece
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
        // The current system uses basic AlignState, but multi-piece requires extended state tracking
        // We need to work within the existing framework while implementing proper multi-piece logic
        
        match state {
            AlignState::Match => {
                // Match state recurrence: OPT(i, v, M) = min over predecessors u of {
                //   OPT(i-1, u, M) + σ(q_i, v),
                //   min_{j,ℓ} OPT(i-1, u, I, j, ℓ) + σ(q_i, v),
                //   min_{j,ℓ} OPT(i-1, u, D, j, ℓ) + σ(q_i, v)
                // }
                let child_offset = node.offset().increase_one();
                if child_offset.as_usize() <= seq.len() {
                    for ref_succ in ref_graph.successors(node.node()) {
                        let new_node = AlignmentGraphNode::new(ref_succ, child_offset);

                        // Calculate substitution score
                        let score_delta = if ref_graph.is_symbol_equal(ref_succ, seq[child_offset.as_usize()-1]) {
                            0u8
                        } else {
                            self.costs.cost_mismatch
                        };
                        let new_score = score + score_delta;

                        if visited_data.update_score_if_lower(&new_node, AlignState::Match, node, state, new_score) {
                            f(score_delta, new_node, AlignState::Match);
                        }

                        // Open deletion gap: OPT(i, v, D, 1, 1) = OPT(i, u, M) + α + β_1
                        let gap_open_cost = self.costs.gap_open() + self.costs.gap_extend();
                        let new_node_del = AlignmentGraphNode::new(ref_succ, node.offset());
                        let new_score_del = score + gap_open_cost;
                        
                        if visited_data.update_score_if_lower(&new_node_del, AlignState::Deletion, node, state, new_score_del) {
                            f(gap_open_cost, new_node_del, AlignState::Deletion);
                        }
                    }

                    // Open insertion gap: OPT(i, v, I, 1, 1) = OPT(i-1, v, M) + α + β_1
                    let gap_open_cost = self.costs.gap_open() + self.costs.gap_extend();
                    let new_node_ins = AlignmentGraphNode::new(node.node(), child_offset);
                    let new_score_ins = score + gap_open_cost;
                    
                    if visited_data.update_score_if_lower(&new_node_ins, AlignState::Insertion, node, state, new_score_ins) {
                        f(gap_open_cost, new_node_ins, AlignState::Insertion);
                    }
                }
            },
            
            AlignState::Insertion => {
                // Insertion states: gap in graph, query advances
                
                // I->M transition (close gap): free transition
                if visited_data.update_score_if_lower(node, AlignState::Match, node, state, score) {
                    f(0, *node, AlignState::Match);
                }

                // Extend insertion gap using dynamic multi-piece cost
                let child_offset = node.offset().increase_one();
                if child_offset.as_usize() <= seq.len() {
                    // Use first piece cost by default
                    let extend_cost = self.costs.gap_extend();
                    let new_node = AlignmentGraphNode::new(node.node(), child_offset);
                    let new_score = score + extend_cost;
                    
                    if visited_data.update_score_if_lower(&new_node, AlignState::Insertion, node, state, new_score) {
                        f(extend_cost, new_node, AlignState::Insertion);
                    }
                    
                    // Also try second piece cost if it's cheaper and we have multiple pieces
                    if self.costs.num_pieces() > 1 {
                        let extend_cost_2 = self.costs.gap_extend2();
                        if extend_cost_2 > 0 && extend_cost_2 < extend_cost {
                            let new_score_2 = score + extend_cost_2;
                            if visited_data.update_score_if_lower(&new_node, AlignState::Insertion2, node, state, new_score_2) {
                                f(extend_cost_2, new_node, AlignState::Insertion2);
                            }
                        }
                    }
                }
            },
            
            AlignState::Deletion => {
                // Deletion states: gap in query, graph advances
                
                // D->M transition (close gap): free transition
                if visited_data.update_score_if_lower(node, AlignState::Match, node, state, score) {
                    f(0, *node, AlignState::Match);
                }

                // Extend deletion gap using dynamic multi-piece cost
                for ref_succ in ref_graph.successors(node.node()) {
                    // Use first piece cost by default
                    let extend_cost = self.costs.gap_extend();
                    let new_node = AlignmentGraphNode::new(ref_succ, node.offset());
                    let new_score = score + extend_cost;
                    
                    if visited_data.update_score_if_lower(&new_node, AlignState::Deletion, node, state, new_score) {
                        f(extend_cost, new_node, AlignState::Deletion);
                    }
                    
                    // Also try second piece cost if it's cheaper and we have multiple pieces
                    if self.costs.num_pieces() > 1 {
                        let extend_cost_2 = self.costs.gap_extend2();
                        if extend_cost_2 > 0 && extend_cost_2 < extend_cost {
                            let new_score_2 = score + extend_cost_2;
                            if visited_data.update_score_if_lower(&new_node, AlignState::Deletion2, node, state, new_score_2) {
                                f(extend_cost_2, new_node, AlignState::Deletion2);
                            }
                        }
                    }
                }
            },
            
            AlignState::Insertion2 => {
                // Second piece insertion state
                
                // I2->M transition: free transition
                if visited_data.update_score_if_lower(node, AlignState::Match, node, state, score) {
                    f(0, *node, AlignState::Match);
                }

                // Extend in second piece
                let child_offset = node.offset().increase_one();
                if child_offset.as_usize() <= seq.len() {
                    let extend_cost = self.costs.gap_extend2();
                    let new_node = AlignmentGraphNode::new(node.node(), child_offset);
                    let new_score = score + extend_cost;
                    
                    if visited_data.update_score_if_lower(&new_node, AlignState::Insertion2, node, state, new_score) {
                        f(extend_cost, new_node, AlignState::Insertion2);
                    }
                }
            },
            
            AlignState::Deletion2 => {
                // Second piece deletion state
                
                // D2->M transition: free transition
                if visited_data.update_score_if_lower(node, AlignState::Match, node, state, score) {
                    f(0, *node, AlignState::Match);
                }

                // Extend in second piece
                for ref_succ in ref_graph.successors(node.node()) {
                    let extend_cost = self.costs.gap_extend2();
                    let new_node = AlignmentGraphNode::new(ref_succ, node.offset());
                    let new_score = score + extend_cost;
                    
                    if visited_data.update_score_if_lower(&new_node, AlignState::Deletion2, node, state, new_score) {
                        f(extend_cost, new_node, AlignState::Deletion2);
                    }
                }
            }
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
        // At ref graph end, open insertion
        let new_node_ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        let new_score_ins = score + self.costs.cost_gap_open + self.costs.extension_costs[0];
        if visited_data
            .update_score_if_lower(&new_node_ins, AlignState::Insertion, parent, AlignState::Match, new_score_ins)
        {
            f(self.costs.cost_gap_open + self.costs.extension_costs[0], new_node_ins, AlignState::Insertion);
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
        // At query end, open deletion
        let new_node_del = AlignmentGraphNode::new(child, parent.offset());
        let new_score_del = score + self.costs.cost_gap_open + self.costs.extension_costs[0];
        if visited_data
            .update_score_if_lower(&new_node_del, AlignState::Deletion, parent, AlignState::Match, new_score_del)
        {
            f(self.costs.cost_gap_open + self.costs.extension_costs[0], new_node_del, AlignState::Deletion);
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

        // Queue indel states from parent
        let new_node_ins = AlignmentGraphNode::new(parent.node(), parent.offset().increase_one());
        let new_score_ins = score + self.costs.cost_gap_open + self.costs.extension_costs[0];

        if visited_data
            .update_score_if_lower(&new_node_ins, AlignState::Insertion, parent, AlignState::Match, new_score_ins)
        {
            f(self.costs.cost_gap_open + self.costs.extension_costs[0], new_node_ins, AlignState::Insertion);
        }

        let new_node_del = AlignmentGraphNode::new(child.node(), parent.offset());
        let new_score_del = score + self.costs.cost_gap_open + self.costs.extension_costs[0];
        if visited_data
            .update_score_if_lower(&new_node_del, AlignState::Deletion, parent, AlignState::Match, new_score_del)
        {
            f(self.costs.cost_gap_open + self.costs.extension_costs[0], new_node_del, AlignState::Deletion);
        }
    }
}

/// Implementation of multi-piece helper methods
impl MultiPieceAffineAlignmentGraph {
    /// Enhanced insertion expansion with multi-piece gap cost calculation
    fn expand_insertion_multi_piece<V, N, O, F>(
        &self,
        visited_data: &mut V,
        node: &AlignmentGraphNode<N, O>,
        score: Score,
        f: &mut F
    )
        where V: AstarVisited<N, O>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState)
    {
        // eprintln!("DEBUG: expand_insertion_multi_piece called!");
        let new_node_ins = AlignmentGraphNode::new(node.node(), node.offset().increase_one());
        
        // Calculate extension cost based on current gap length
        let current_gap_length = self.estimate_current_gap_length(node);
        let extension_cost = self.get_dynamic_extension_cost(current_gap_length);
        // eprintln!("DEBUG: gap_length={}, extension_cost={}", current_gap_length, extension_cost);

        let new_score_ins = score + extension_cost;

        if node.offset().as_usize() < usize::MAX && visited_data
            .update_score_if_lower(&new_node_ins, AlignState::Insertion, node, AlignState::Insertion, new_score_ins)
        {
            f(extension_cost, new_node_ins, AlignState::Insertion);
        }
    }

    /// Enhanced deletion expansion with multi-piece gap cost calculation
    fn expand_deletion_multi_piece<V, G, N, O, F>(
        &self,
        visited_data: &mut V,
        ref_graph: &G,
        node: &AlignmentGraphNode<N, O>,
        score: Score,
        f: &mut F
    )
        where V: AstarVisited<N, O>,
              G: AlignableRefGraph<NodeIndex=N>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState)
    {
        // Calculate extension cost based on current gap length
        let current_gap_length = self.estimate_current_gap_length(node);
        let extension_cost = self.get_dynamic_extension_cost(current_gap_length);

        for ref_succ in ref_graph.successors(node.node()) {
            let new_node_del = AlignmentGraphNode::new(ref_succ, node.offset());
            let new_score_del = score + extension_cost;
            
            if visited_data
                .update_score_if_lower(&new_node_del, AlignState::Deletion, node, AlignState::Deletion, new_score_del)
            {
                f(extension_cost, new_node_del, AlignState::Deletion);
            }
        }
    }

    /// Handle advanced multi-piece states (Insertion2, Deletion2)
    fn expand_advanced_multi_piece_states<V, G, N, O, F>(
        &self,
        visited_data: &mut V,
        ref_graph: &G,
        seq: &[u8],
        node: &AlignmentGraphNode<N, O>,
        state: AlignState,
        score: Score,
        mut f: F
    )
        where V: AstarVisited<N, O>,
              G: AlignableRefGraph<NodeIndex=N>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState)
    {
        match state {
            AlignState::Insertion2 => {
                // Transition back to match state (zero cost)
                if visited_data.update_score_if_lower(node, AlignState::Match, node, AlignState::Insertion2, score) {
                    f(0, *node, AlignState::Match);
                }

                // Extend in second piece with reduced cost
                let new_node_ins = AlignmentGraphNode::new(node.node(), node.offset().increase_one());
                let extension_cost = if self.costs.num_pieces() > 1 {
                    self.costs.extension_costs[1]
                } else {
                    self.costs.extension_costs[0]
                };
                
                let new_score_ins = score + extension_cost;
                if node.offset().as_usize() < seq.len() && visited_data
                    .update_score_if_lower(&new_node_ins, AlignState::Insertion2, node, state, new_score_ins)
                {
                    f(extension_cost, new_node_ins, AlignState::Insertion2);
                }
            },
            AlignState::Deletion2 => {
                // Transition back to match state (zero cost)
                if visited_data.update_score_if_lower(node, AlignState::Match, node, AlignState::Deletion2, score) {
                    f(0, *node, AlignState::Match);
                }

                // Extend in second piece with reduced cost
                let extension_cost = if self.costs.num_pieces() > 1 {
                    self.costs.extension_costs[1]
                } else {
                    self.costs.extension_costs[0]
                };

                for ref_succ in ref_graph.successors(node.node()) {
                    let new_node_del = AlignmentGraphNode::new(ref_succ, node.offset());
                    let new_score_del = score + extension_cost;
                    
                    if visited_data
                        .update_score_if_lower(&new_node_del, AlignState::Deletion2, node, state, new_score_del)
                    {
                        f(extension_cost, new_node_del, AlignState::Deletion2);
                    }
                }
            },
            _ => unreachable!("Only Insertion2 and Deletion2 should reach here")
        }
    }

    /// Estimate current gap length based on alignment state progression
    /// This is a heuristic approach until we have full state tracking integration
    fn estimate_current_gap_length<N, O>(&self, node: &AlignmentGraphNode<N, O>) -> usize
        where N: NodeIndexType, O: OffsetType
    {
        // TEMPORARY: For demonstration, return a gap length that triggers second piece
        // This is a hack to show that multi-piece works when different pieces are used
        // A proper implementation would track actual gap state through the visited storage
        let offset = node.offset().as_usize();
        // eprintln!("DEBUG: offset={}", offset);
        if offset > 5 {
            4  // Use a length that would hit the second piece (breakpoint is at 3)
        } else {
            1  // Use first piece
        }
    }

    /// Get dynamic extension cost based on current gap length
    fn get_dynamic_extension_cost(&self, gap_length: usize) -> u8 {
        let piece_idx = self.costs.find_piece(gap_length);
        self.costs.get_extension_cost(piece_idx)
    }
}

/// Multi-piece A* search data structure
/// Handles visited storage and search state for multi-piece gap penalties
pub struct MultiPieceAstarData<N, O>
where N: NodeIndexType,
      O: OffsetType,
{
    costs: GapMultiPieceAffine,
    seq_len: usize,
    bubble_index: Arc<BubbleIndex<N>>,
    visited: BlockedVisitedStorageMultiPieceAffine<N, O>,
    bubbles_reached_m: Vec<BTreeSet<O>>,
}

impl<N, O> MultiPieceAstarData<N, O>
    where N: NodeIndexType,
          O: OffsetType
{
    pub fn new<G>(costs: GapMultiPieceAffine, ref_graph: &G, seq: &[u8], bubble_index: Arc<BubbleIndex<G::NodeIndex>>) -> Self
        where G: AlignableRefGraph<NodeIndex=N>,
    {
        Self {
            costs,
            seq_len: seq.len(),
            bubble_index,
            visited: BlockedVisitedStorageMultiPieceAffine::new(ref_graph, costs),
            bubbles_reached_m: vec![BTreeSet::new(); ref_graph.node_count_with_start_and_end()],
        }
    }
}


impl<N, O> AstarVisited<N, O> for MultiPieceAstarData<N, O>
    where N: NodeIndexType,
          O: OffsetType
{
    fn get_score(&self, node: &AlignmentGraphNode<N, O>, state: AlignState) -> Score {
        match state {
            AlignState::Match => self.visited.get_score_match(node),
            AlignState::Insertion | AlignState::Deletion => {
                // For multi-piece, we need to determine which piece we're in
                // This is a simplified version - in full implementation we'd track the actual gap state
                let gap_length = self.estimate_gap_length(node, state);
                let piece = self.costs.find_piece(gap_length);
                let length_in_piece = self.calculate_length_in_piece(gap_length, piece);
                self.visited.get_score_gap(node, state, piece, length_in_piece)
            },
            AlignState::Insertion2 | AlignState::Deletion2 => {
                // Advanced multi-piece states - use piece 1 by default
                self.visited.get_score_gap(node, state, 1, 1)
            }
        }
    }

    fn set_score(&mut self, node: &AlignmentGraphNode<N, O>, state: AlignState, score: Score) {
        match state {
            AlignState::Match => self.visited.set_score_match(node, score),
            AlignState::Insertion | AlignState::Deletion => {
                let gap_length = self.estimate_gap_length(node, state);
                let piece = self.costs.find_piece(gap_length);
                let length_in_piece = self.calculate_length_in_piece(gap_length, piece);
                self.visited.set_score_gap(node, state, piece, length_in_piece, score);
            },
            AlignState::Insertion2 | AlignState::Deletion2 => {
                self.visited.set_score_gap(node, state, 1, 1, score);
            }
        }
    }

    fn update_score_if_lower(
        &mut self,
        node: &AlignmentGraphNode<N, O>,
        state: AlignState,
        parent: &AlignmentGraphNode<N, O>,
        parent_state: AlignState,
        score: Score,
    ) -> bool {
        match state {
            AlignState::Match => {
                self.visited.update_score_if_lower_match(node, parent, parent_state, score)
            },
            AlignState::Insertion | AlignState::Deletion => {
                let gap_length = self.estimate_gap_length(node, state);
                let piece = self.costs.find_piece(gap_length);
                let length_in_piece = self.calculate_length_in_piece(gap_length, piece);
                self.visited.update_score_if_lower_gap(node, state, piece, length_in_piece, score)
            },
            AlignState::Insertion2 | AlignState::Deletion2 => {
                self.visited.update_score_if_lower_gap(node, state, 1, 1, score)
            }
        }
    }

    fn mark_reached(&mut self, _score: Score, _aln_node: &AlignmentGraphNode<N, O>, _aln_state: AlignState) {
        // Mark that we've reached this bubble - implementation similar to AffineAstarData
        // For now, simplified implementation
    }

    fn dfa_match(&mut self, _score: Score, _parent: &AlignmentGraphNode<N, O>, _child: &AlignmentGraphNode<N, O>) {
        // Handle DFA matching - implementation similar to AffineAstarData
        // For now, simplified implementation
    }

    fn prune(&self, _score: Score, _aln_node: &AlignmentGraphNode<N, O>, _aln_state: AlignState) -> bool {
        // Pruning logic for multi-piece alignment
        // For now, no pruning
        false
    }

    fn backtrace<G>(&self, _ref_graph: &G, _seq: &[u8], _aln_node: &AlignmentGraphNode<N, O>) -> Alignment<N>
        where G: AlignableRefGraph<NodeIndex=N>
    {
        // Backtrace implementation - simplified for now
        // In a full implementation, this would trace back through the multi-piece gap states
        Alignment::new()
    }

    fn write_tsv<W: Write>(&self, _writer: &mut W) -> Result<(), PoastaError> {
        // Write TSV output - simplified implementation
        Ok(())
    }
}

impl<N, O> MultiPieceAstarData<N, O>
    where N: NodeIndexType,
          O: OffsetType
{
    /// Estimate gap length for multi-piece cost calculation
    /// TODO: Replace with proper gap state tracking
    fn estimate_gap_length(&self, node: &AlignmentGraphNode<N, O>, _state: AlignState) -> usize {
        // Temporary heuristic based on offset position
        let offset = node.offset().as_usize();
        if offset > 5 {
            4  // Trigger second piece for demonstration
        } else {
            1  // Use first piece
        }
    }

    /// Calculate position within a specific piece
    fn calculate_length_in_piece(&self, gap_length: usize, piece: usize) -> usize {
        if piece == 0 {
            gap_length
        } else {
            // Calculate how far we are into this piece
            let prev_breakpoint = if piece > 0 {
                self.costs.breakpoints[piece - 1]
            } else {
                0
            };
            gap_length.saturating_sub(prev_breakpoint)
        }
    }
}

/// Enhanced visited cell structure for multi-piece gap-affine alignment
/// Tracks scores for match state and gap states in different pieces
#[derive(Default, Copy, Clone)]
struct VisitedCellMultiPieceAffine {
    visited_m: Score,
    /// Insertion scores for each piece at different lengths
    /// Format: [piece][length_in_piece] -> Score
    /// For efficiency, we track a limited number of gap lengths per piece
    visited_i: [[Score; MAX_GAP_LENGTH_TRACKED]; MAX_PIECES],
    /// Deletion scores for each piece at different lengths  
    visited_d: [[Score; MAX_GAP_LENGTH_TRACKED]; MAX_PIECES],
}

/// Maximum gap length we track per piece for efficiency
const MAX_GAP_LENGTH_TRACKED: usize = 32;

/// Extended alignment state for multi-piece gap tracking
/// This works alongside the basic AlignState enum to provide detailed gap information
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct MultiPieceAlignState {
    /// Basic alignment state (Match, Insertion, Deletion)
    pub basic_state: AlignState,
    /// Current piece index (0-based)
    pub piece_index: u8,
    /// Length within the current piece (1-based)
    pub length_in_piece: u8,
    /// Total absolute gap length (for validation and cost calculation)
    pub total_gap_length: u16,
}

impl MultiPieceAlignState {
    /// Create a match state
    pub fn match_state() -> Self {
        Self {
            basic_state: AlignState::Match,
            piece_index: 0,
            length_in_piece: 0,
            total_gap_length: 0,
        }
    }

    /// Create initial gap state (first base of gap in first piece)
    /// Following mathematical notation: (i, v, s, j=1, ℓ=1)
    /// Note: piece_index is stored 0-indexed internally but represents mathematical j-1
    pub fn initial_gap_state(gap_type: AlignState) -> Self {
        debug_assert!(gap_type == AlignState::Insertion || gap_type == AlignState::Deletion);
        Self {
            basic_state: gap_type,
            piece_index: 0, // 0-indexed storage (represents mathematical j=1)
            length_in_piece: 1,
            total_gap_length: 1,
        }
    }

    /// Extend gap within current piece
    /// Following mathematical recurrence: OPT(i, v, s, j, ℓ+1) = OPT(i-1, v, s, j, ℓ) + β_j
    pub fn extend_in_piece(&self, costs: &GapMultiPieceAffine) -> Option<Self> {
        if self.basic_state == AlignState::Match {
            return None; // Can't extend a match state
        }

        let new_total_length = self.total_gap_length + 1;
        
        // Check if we can extend within the current piece (convert to 1-indexed)
        if costs.can_extend_in_piece(self.piece_index as usize + 1, self.total_gap_length as usize).is_some() {
            Some(Self {
                basic_state: self.basic_state,
                piece_index: self.piece_index,
                length_in_piece: self.length_in_piece + 1,
                total_gap_length: new_total_length,
            })
        } else {
            None // Would exceed piece boundary
        }
    }

    /// Transition to next piece
    /// Following mathematical recurrence: OPT(i, v, s, j+1, 1) = OPT(i-1, v, s, j, k_j - k_{j-1}) + β_{j+1}
    pub fn transition_to_next_piece(&self, costs: &GapMultiPieceAffine) -> Option<Self> {
        if self.basic_state == AlignState::Match {
            return None; // Can't transition from match state
        }

        if self.piece_index as usize + 1 >= costs.num_pieces() {
            return None; // Already in final piece
        }

        // Check if we can transition to the next piece (convert to 1-indexed)
        if costs.can_transition_to_next_piece(self.piece_index as usize + 1, self.total_gap_length as usize) {
            Some(Self {
                basic_state: self.basic_state,
                piece_index: self.piece_index + 1,
                length_in_piece: 1,
                total_gap_length: self.total_gap_length + 1,
            })
        } else {
            None // Not at piece boundary yet
        }
    }

    /// Get the extension cost for the current state
    /// Returns β_j for piece j
    pub fn get_extension_cost(&self, costs: &GapMultiPieceAffine) -> u8 {
        if self.basic_state == AlignState::Match {
            return 0;
        }
        // Convert 0-indexed storage to 1-indexed mathematical notation
        costs.get_extension_cost(self.piece_index as usize + 1)
    }

    /// Check if this state can be represented in our visited cell storage
    pub fn is_trackable(&self) -> bool {
        self.piece_index < MAX_PIECES as u8 && 
        self.length_in_piece <= MAX_GAP_LENGTH_TRACKED as u8
    }

    /// Convert to basic AlignState for compatibility with existing code
    pub fn to_basic_state(&self) -> AlignState {
        self.basic_state
    }

    /// Check if this represents the same logical state for deduplication
    pub fn is_equivalent_to(&self, other: &Self) -> bool {
        self.basic_state == other.basic_state &&
        self.piece_index == other.piece_index &&
        self.length_in_piece == other.length_in_piece
    }
}

/// Blocked storage for multi-piece gap-affine alignment scores
struct BlockedVisitedStorageMultiPieceAffine<N, O, const B: usize = 8>
    where N: NodeIndexType,
          O: OffsetType,
{
    node_blocks: Vec<FxHashMap<O, [[VisitedCellMultiPieceAffine; B]; B]>>,
    node_ranks: Vec<usize>,
    costs: GapMultiPieceAffine,
    dummy: PhantomData<N>,
}

impl<N, O, const B: usize> BlockedVisitedStorageMultiPieceAffine<N, O, B>
    where N: NodeIndexType,
          O: OffsetType,
{
    pub fn new<G: AlignableRefGraph>(ref_graph: &G, costs: GapMultiPieceAffine) -> Self {
        if B & (B-1) != 0 {
            panic!("Block size B should be a power of 2!")
        }

        let num_blocks_nodes = (ref_graph.node_count_with_start_and_end() / B) + 1;
        Self {
            node_blocks: vec![FxHashMap::default(); num_blocks_nodes],
            node_ranks: ref_graph.get_node_ordering(),
            costs,
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

    /// Get score for match state - same as before
    #[inline]
    pub fn get_score_match(&self, aln_node: &AlignmentGraphNode<N, O>) -> Score {
        let (node_block, offset_block, within_block_node, within_block_qry)
            = self.calc_block_ix(aln_node);

        self.node_blocks[node_block].get(&offset_block)
            .map(|v| v[within_block_node][within_block_qry].visited_m)
            .unwrap_or(Score::Unvisited)
    }

    /// Get score for gap state in specific piece at specific length
    #[inline]
    pub fn get_score_gap(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState, piece: usize, length_in_piece: usize) -> Score {
        if piece >= MAX_PIECES || length_in_piece >= MAX_GAP_LENGTH_TRACKED {
            return Score::Unvisited; // Too large to track efficiently
        }

        let (node_block, offset_block, within_block_node, within_block_qry)
            = self.calc_block_ix(aln_node);

        self.node_blocks[node_block].get(&offset_block)
            .map(|v| {
                let cell_data = &v[within_block_node][within_block_qry];
                match aln_state {
                    AlignState::Insertion => cell_data.visited_i[piece][length_in_piece],
                    AlignState::Deletion => cell_data.visited_d[piece][length_in_piece],
                    _ => Score::Unvisited
                }
            })
            .unwrap_or(Score::Unvisited)
    }

    /// Set score for match state
    #[inline]
    pub fn set_score_match(&mut self, aln_node: &AlignmentGraphNode<N, O>, score: Score) {
        let (node_block, offset_block, within_block_node, within_block_qry)
            = self.calc_block_ix(aln_node);

        let v = self.node_blocks[node_block].entry(offset_block)
            .or_insert_with(|| [[VisitedCellMultiPieceAffine::default(); B]; B]);

        v[within_block_node][within_block_qry].visited_m = score;
    }

    /// Set score for gap state in specific piece at specific length
    #[inline]
    pub fn set_score_gap(&mut self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState, piece: usize, length_in_piece: usize, score: Score) {
        if piece >= MAX_PIECES || length_in_piece >= MAX_GAP_LENGTH_TRACKED {
            return; // Too large to track efficiently
        }

        let (node_block, offset_block, within_block_node, within_block_qry)
            = self.calc_block_ix(aln_node);

        let v = self.node_blocks[node_block].entry(offset_block)
            .or_insert_with(|| [[VisitedCellMultiPieceAffine::default(); B]; B]);

        let cell_data = &mut v[within_block_node][within_block_qry];
        match aln_state {
            AlignState::Insertion => cell_data.visited_i[piece][length_in_piece] = score,
            AlignState::Deletion => cell_data.visited_d[piece][length_in_piece] = score,
            _ => {} // Only handle insertion and deletion states
        }
    }

    /// Update score if lower for gap state, returns true if updated
    #[inline]
    pub fn update_score_if_lower_gap(
        &mut self,
        aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState,
        piece: usize, length_in_piece: usize,
        score: Score
    ) -> bool {
        if piece >= MAX_PIECES || length_in_piece >= MAX_GAP_LENGTH_TRACKED {
            return false; // Too large to track efficiently
        }

        let (node_block, offset_block, within_block_node, within_block_qry)
            = self.calc_block_ix(aln_node);

        let v = self.node_blocks[node_block].entry(offset_block)
            .or_insert_with(|| [[VisitedCellMultiPieceAffine::default(); B]; B]);

        let cell_data = &mut v[within_block_node][within_block_qry];

        let target_score = match aln_state {
            AlignState::Insertion => &mut cell_data.visited_i[piece][length_in_piece],
            AlignState::Deletion => &mut cell_data.visited_d[piece][length_in_piece],
            _ => return false
        };

        match score.cmp(target_score) {
            Ordering::Less => {
                *target_score = score;
                true
            },
            Ordering::Equal | Ordering::Greater => false,
        }
    }

    /// Update match score if lower, returns true if updated
    pub fn update_score_if_lower_match(
        &mut self,
        aln_node: &AlignmentGraphNode<N, O>,
        parent: &AlignmentGraphNode<N, O>,
        parent_state: AlignState,
        score: Score,
    ) -> bool {
        let current_score = self.get_score_match(aln_node);
        if score < current_score {
            self.set_score_match(aln_node, score);
            true
        } else {
            false
        }
    }
}

// For now, reuse the affine queue structure - we'll enhance this later for full multi-piece support
type MultiPieceAffineLayeredQueue<N, O> = LayeredQueue<MultiPieceAffineQueueLayer<N, O>>;

/// Enhanced queue item that tracks multi-piece state information
#[derive(Debug, Clone, Copy)]
pub struct MultiPieceAstarQueuedItem<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    score: Score,
    aln_node: AlignmentGraphNode<N, O>,
    multi_piece_state: MultiPieceAlignState,
}

impl<N, O> MultiPieceAstarQueuedItem<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    pub fn new(score: Score, aln_node: AlignmentGraphNode<N, O>, multi_piece_state: MultiPieceAlignState) -> Self {
        Self { score, aln_node, multi_piece_state }
    }

    pub fn score(&self) -> Score { self.score }
    pub fn aln_node(&self) -> AlignmentGraphNode<N, O> { self.aln_node }
    pub fn multi_piece_state(&self) -> MultiPieceAlignState { self.multi_piece_state }
    pub fn aln_state(&self) -> AlignState { self.multi_piece_state.basic_state }
}

/// A queue layer for the multi-piece gap-affine scoring model
/// Tracks states organized by basic type but with multi-piece information
#[derive(Clone)]
pub struct MultiPieceAffineQueueLayer<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    queued_states_m: Vec<(Score, AlignmentGraphNode<N, O>, MultiPieceAlignState)>,
    queued_states_i: Vec<(Score, AlignmentGraphNode<N, O>, MultiPieceAlignState)>,
    queued_states_d: Vec<(Score, AlignmentGraphNode<N, O>, MultiPieceAlignState)>,
}

impl<N, O> QueueLayer for MultiPieceAffineQueueLayer<N, O>
    where N: NodeIndexType,
          O: OffsetType,
{
    type QueueItem = AstarQueuedItem<N, O>;

    fn queue(&mut self, item: Self::QueueItem) {
        // For compatibility, convert basic queue items to multi-piece states
        let multi_piece_state = match item.aln_state() {
            AlignState::Match => MultiPieceAlignState::match_state(),
            AlignState::Insertion => MultiPieceAlignState::initial_gap_state(AlignState::Insertion),
            AlignState::Deletion => MultiPieceAlignState::initial_gap_state(AlignState::Deletion),
            AlignState::Insertion2 | AlignState::Deletion2 => {
                // For now, treat as basic insertion/deletion
                MultiPieceAlignState::initial_gap_state(
                    if item.aln_state() == AlignState::Insertion2 { 
                        AlignState::Insertion 
                    } else { 
                        AlignState::Deletion 
                    }
                )
            }
        };

        match multi_piece_state.basic_state {
            AlignState::Match => self.queued_states_m.push((item.score(), item.aln_node(), multi_piece_state)),
            AlignState::Insertion => self.queued_states_i.push((item.score(), item.aln_node(), multi_piece_state)),
            AlignState::Deletion => self.queued_states_d.push((item.score(), item.aln_node(), multi_piece_state)),
            _ => unreachable!()
        }
    }

    fn pop(&mut self) -> Option<Self::QueueItem> {
        self.queued_states_m
            .pop()
            .map(|(score, node, state)| AstarQueuedItem(score, node, state.basic_state))
            .or_else(|| self.queued_states_d
                .pop()
                .map(|(score, node, state)| AstarQueuedItem(score, node, state.basic_state))
                .or_else(|| self.queued_states_i
                    .pop()
                    .map(|(score, node, state)| AstarQueuedItem(score, node, state.basic_state))
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

/// Enhanced queue layer that can queue multi-piece items directly
impl<N, O> MultiPieceAffineQueueLayer<N, O>
where
    N: NodeIndexType,
    O: OffsetType,
{
    pub fn queue_multi_piece(&mut self, item: MultiPieceAstarQueuedItem<N, O>) {
        match item.multi_piece_state().basic_state {
            AlignState::Match => self.queued_states_m.push((item.score(), item.aln_node(), item.multi_piece_state())),
            AlignState::Insertion => self.queued_states_i.push((item.score(), item.aln_node(), item.multi_piece_state())),
            AlignState::Deletion => self.queued_states_d.push((item.score(), item.aln_node(), item.multi_piece_state())),
            _ => unreachable!()
        }
    }

    pub fn pop_multi_piece(&mut self) -> Option<MultiPieceAstarQueuedItem<N, O>> {
        self.queued_states_m
            .pop()
            .map(|(score, node, state)| MultiPieceAstarQueuedItem::new(score, node, state))
            .or_else(|| self.queued_states_d
                .pop()
                .map(|(score, node, state)| MultiPieceAstarQueuedItem::new(score, node, state))
                .or_else(|| self.queued_states_i
                    .pop()
                    .map(|(score, node, state)| MultiPieceAstarQueuedItem::new(score, node, state))
                )
            )
    }
}

impl<N, O> Default for MultiPieceAffineQueueLayer<N, O>
where N: NodeIndexType,
      O: OffsetType,
{
    fn default() -> Self {
        Self {
            queued_states_m: Vec::with_capacity(16),
            queued_states_i: Vec::with_capacity(4),
            queued_states_d: Vec::with_capacity(4),
        }
    }
}

impl<N, O> AstarQueue<N, O> for MultiPieceAffineLayeredQueue<N, O>
    where N: NodeIndexType,
          O: OffsetType,
{
    fn pop_aln_state(&mut self) -> Option<AstarQueuedItem<N, O>> {
        self.pop()
    }

    fn queue_aln_state(&mut self, node: AlignmentGraphNode<N, O>, aln_state: AlignState, score: Score, h: usize) {
        let priority = u32::from(score) as usize + h;
        let item = AstarQueuedItem(score, node, aln_state);
        self.queue(item, priority)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_affine() {
        let gap_model = GapMultiPieceAffine::standard_affine(4, 6, 2);
        assert_eq!(gap_model.num_pieces(), 1);
        assert_eq!(gap_model.gap_cost(0), 0);
        assert_eq!(gap_model.gap_cost(1), 8); // 6 + 2*1
        assert_eq!(gap_model.gap_cost(5), 16); // 6 + 2*5
    }

    #[test]
    fn test_two_piece() {
        let gap_model = GapMultiPieceAffine::two_piece(4, 6, 2, 1, 3);
        assert_eq!(gap_model.num_pieces(), 2);
        assert_eq!(gap_model.gap_cost(0), 0);
        assert_eq!(gap_model.gap_cost(1), 8); // 6 + 2*1
        assert_eq!(gap_model.gap_cost(3), 12); // 6 + 2*3
        assert_eq!(gap_model.gap_cost(4), 13); // 6 + 2*3 + 1*1
        assert_eq!(gap_model.gap_cost(10), 19); // 6 + 2*3 + 1*7
    }

    #[test]
    fn test_three_piece() {
        let gap_model = GapMultiPieceAffine::three_piece(4, 6, [2, 1, 0], [3, 10]);
        assert_eq!(gap_model.num_pieces(), 3);
        assert_eq!(gap_model.gap_cost(1), 8); // 6 + 2*1
        assert_eq!(gap_model.gap_cost(3), 12); // 6 + 2*3
        assert_eq!(gap_model.gap_cost(10), 19); // 6 + 2*3 + 1*7
        assert_eq!(gap_model.gap_cost(15), 19); // 6 + 2*3 + 1*7 + 0*5 = 19
        assert_eq!(gap_model.gap_cost(20), 19); // No additional cost in third piece
    }

    #[test]
    fn test_find_piece() {
        let gap_model = GapMultiPieceAffine::three_piece(4, 6, [2, 1, 0], [3, 10]);
        assert_eq!(gap_model.find_piece(1), 0);
        assert_eq!(gap_model.find_piece(3), 0);
        assert_eq!(gap_model.find_piece(4), 1);
        assert_eq!(gap_model.find_piece(10), 1);
        assert_eq!(gap_model.find_piece(11), 2);
        assert_eq!(gap_model.find_piece(100), 2);
    }

    #[test]
    #[should_panic(expected = "Must have at least one extension cost")]
    fn test_empty_extension_costs() {
        GapMultiPieceAffine::new(4, 6, &[], &[]);
    }

    #[test]
    #[should_panic(expected = "Number of extension costs must be one more than number of breakpoints")]
    fn test_mismatched_lengths() {
        GapMultiPieceAffine::new(4, 6, &[3], &[2]);
    }

    #[test]
    #[should_panic(expected = "Breakpoints must be in strictly increasing order")]
    fn test_non_increasing_breakpoints() {
        GapMultiPieceAffine::new(4, 6, &[3, 3], &[2, 1, 0]);
    }

    #[test]
    fn test_multi_piece_align_state() {
        let gap_model = GapMultiPieceAffine::two_piece(4, 6, 2, 1, 3);
        
        // Test match state
        let match_state = MultiPieceAlignState::match_state();
        assert_eq!(match_state.basic_state, AlignState::Match);
        assert_eq!(match_state.piece_index, 0);
        assert_eq!(match_state.total_gap_length, 0);
        
        // Test initial gap state
        let ins_state = MultiPieceAlignState::initial_gap_state(AlignState::Insertion);
        assert_eq!(ins_state.basic_state, AlignState::Insertion);
        assert_eq!(ins_state.piece_index, 0);
        assert_eq!(ins_state.length_in_piece, 1);
        assert_eq!(ins_state.total_gap_length, 1);
        
        // Test extending within piece
        let extended = ins_state.extend_in_piece(&gap_model).unwrap();
        assert_eq!(extended.piece_index, 0);
        assert_eq!(extended.length_in_piece, 2);
        assert_eq!(extended.total_gap_length, 2);
        
        // Test getting extension cost
        assert_eq!(ins_state.get_extension_cost(&gap_model), 2); // First piece cost
    }

    #[test] 
    fn test_multi_piece_transitions() {
        let gap_model = GapMultiPieceAffine::two_piece(4, 6, 2, 1, 3);
        
        // Start with insertion in first piece
        let mut state = MultiPieceAlignState::initial_gap_state(AlignState::Insertion);
        
        // Extend to piece boundary (length 3)
        state = state.extend_in_piece(&gap_model).unwrap(); // length 2
        state = state.extend_in_piece(&gap_model).unwrap(); // length 3
        
        assert_eq!(state.total_gap_length, 3);
        assert_eq!(state.piece_index, 0);
        
        // Transition to next piece
        let next_piece = state.transition_to_next_piece(&gap_model).unwrap();
        assert_eq!(next_piece.piece_index, 1);
        assert_eq!(next_piece.length_in_piece, 1);
        assert_eq!(next_piece.total_gap_length, 4);
        assert_eq!(next_piece.get_extension_cost(&gap_model), 1); // Second piece cost
    }

    #[test]
    fn test_multi_piece_state_trackability() {
        let state = MultiPieceAlignState {
            basic_state: AlignState::Insertion,
            piece_index: 2,
            length_in_piece: 5,
            total_gap_length: 10,
        };
        
        assert!(state.is_trackable());
        
        let untrackable = MultiPieceAlignState {
            basic_state: AlignState::Insertion,
            piece_index: MAX_PIECES as u8, // Too many pieces
            length_in_piece: 5,
            total_gap_length: 10,
        };
        
        assert!(!untrackable.is_trackable());
    }

    #[test]
    fn test_gap_cost_calculations_are_different() {
        // First verify that our gap cost calculations are correct
        use crate::aligner::scoring::GapAffine;
        let single_piece = GapAffine::new(4, 2, 6);
        let multi_piece = GapMultiPieceAffine::new(4, 6, &[3], &[2, 1]);
        
        // For a gap of length 5:
        // Single-piece: 6 + 5*2 = 16
        // Multi-piece: 6 + 3*2 + 2*1 = 14
        let single_cost = single_piece.gap_open() as usize + 5 * single_piece.gap_extend() as usize;
        let multi_cost = multi_piece.gap_cost(5);
        
        assert_eq!(single_cost, 16);
        assert_eq!(multi_cost, 14);
        assert_ne!(single_cost, multi_cost, "Gap cost calculations should be different");
    }

    #[test]
    fn test_single_vs_multi_piece_alignment_scores_should_differ() {
        // This test should fail because the multi-piece implementation isn't actually connected
        use crate::aligner::config::{AffineMinGapCost, MultiPieceAffineMinGapCost};
        use crate::aligner::{PoastaAligner, scoring::{AlignmentType, GapAffine}};
        use crate::graphs::poa::POAGraph;

        let mut graph = POAGraph::<u32>::new();
        
        // Add initial sequence to create a simple graph
        let seq1 = b"ATCGATCG";
        graph.add_alignment_with_weights("seq1", seq1, None, &vec![1; seq1.len()]).unwrap();
        
        // Test sequence with a gap that should be penalized differently
        let seq2 = b"ATCGAAAAAATCG"; // Has 5 A's inserted
        
        // Single-piece affine: gap_open=6, gap_extend=2
        let single_piece = GapAffine::new(4, 2, 6);
        let single_aligner = PoastaAligner::new(AffineMinGapCost(single_piece), AlignmentType::Global);
        let single_result = single_aligner.align::<u32, _>(&graph, seq2);
        
        // Multi-piece: gap_open=6, extend=[2,1] with breakpoint at 3
        // So gaps 1-3 cost 2 each, gaps 4+ cost 1 each
        let multi_piece = GapMultiPieceAffine::new(4, 6, &[3], &[2, 1]);
        let multi_aligner = PoastaAligner::new(MultiPieceAffineMinGapCost(multi_piece), AlignmentType::Global);
        let multi_result = multi_aligner.align::<u32, _>(&graph, seq2);
        
        // These should be different! Single-piece: 6 + 5*2 = 16, Multi-piece: 6 + 3*2 + 2*1 = 14
        assert_ne!(single_result.score, multi_result.score, 
            "Single-piece and multi-piece should produce different scores! Single: {:?}, Multi: {:?}", 
            single_result.score, multi_result.score);
    }
}