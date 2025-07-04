use super::offsets::OffsetType;
use crate::graphs::{AlignableRefGraph, NodeIndexType};
use std::fmt::Debug;
use crate::aligner::astar::AstarVisited;
use crate::aligner::scoring::{AlignmentCosts, Score};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AlignState {
    Match,
    Deletion,
    Insertion,
    Deletion2, // For two-piece gap model (backward compatibility)
    Insertion2,
    // Multi-piece gap states with position tracking
    MultiInsertion { piece: u8, position: u8 },
    MultiDeletion { piece: u8, position: u8 },
}

impl AlignState {
    /// Check if this is a gap state (insertion or deletion)
    pub fn is_gap(&self) -> bool {
        matches!(self, 
            AlignState::Insertion | AlignState::Deletion | 
            AlignState::Insertion2 | AlignState::Deletion2 |
            AlignState::MultiInsertion { .. } | AlignState::MultiDeletion { .. }
        )
    }
    
    /// Check if this is an insertion state
    pub fn is_insertion(&self) -> bool {
        matches!(self, 
            AlignState::Insertion | AlignState::Insertion2 | AlignState::MultiInsertion { .. }
        )
    }
    
    /// Check if this is a deletion state  
    pub fn is_deletion(&self) -> bool {
        matches!(self, 
            AlignState::Deletion | AlignState::Deletion2 | AlignState::MultiDeletion { .. }
        )
    }
    
    /// Get piece and position for multi-piece states
    pub fn get_piece_info(&self) -> Option<(u8, u8)> {
        match self {
            AlignState::MultiInsertion { piece, position } | 
            AlignState::MultiDeletion { piece, position } => Some((*piece, *position)),
            _ => None,
        }
    }
    
    /// Create a new multi-piece insertion state
    pub fn multi_insertion(piece: u8, position: u8) -> Self {
        AlignState::MultiInsertion { piece, position }
    }
    
    /// Create a new multi-piece deletion state
    pub fn multi_deletion(piece: u8, position: u8) -> Self {
        AlignState::MultiDeletion { piece, position }
    }
    
    /// Convert to next position in same piece, or transition to next piece
    pub fn advance_gap_position(&self, should_transition: bool, next_piece: u8) -> Option<Self> {
        match self {
            AlignState::MultiInsertion { piece, position } => {
                if should_transition {
                    Some(AlignState::MultiInsertion { piece: next_piece, position: 1 })
                } else {
                    Some(AlignState::MultiInsertion { piece: *piece, position: position + 1 })
                }
            },
            AlignState::MultiDeletion { piece, position } => {
                if should_transition {
                    Some(AlignState::MultiDeletion { piece: next_piece, position: 1 })
                } else {
                    Some(AlignState::MultiDeletion { piece: *piece, position: position + 1 })
                }
            },
            _ => None,
        }
    }
}

/// Unified alignment state for generalized multi-piece gap penalties
/// Combines basic alignment operation with piece tracking information
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GeneralizedAlignState {
    /// Basic alignment operation
    pub basic_state: BasicAlignState,
    /// Current piece index (1-based, 0 means not in gap)
    pub piece_index: u8,
    /// Position within current piece (1-based)
    pub piece_position: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BasicAlignState {
    Match,
    Insertion,
    Deletion,
}

impl GeneralizedAlignState {
    /// Create a match state
    pub const fn match_state() -> Self {
        Self {
            basic_state: BasicAlignState::Match,
            piece_index: 0,
            piece_position: 0,
        }
    }
    
    /// Create a new gap state (first position of first piece)
    pub const fn new_gap_state(gap_type: BasicAlignState) -> Self {
        Self {
            basic_state: gap_type,
            piece_index: 1,
            piece_position: 1,
        }
    }
    
    /// Create a gap state at specific piece and position
    pub const fn gap_state(gap_type: BasicAlignState, piece_index: u8, piece_position: u8) -> Self {
        Self {
            basic_state: gap_type,
            piece_index,
            piece_position,
        }
    }
    
    /// Check if this is a match state
    pub const fn is_match(&self) -> bool {
        matches!(self.basic_state, BasicAlignState::Match)
    }
    
    /// Check if this is a gap state (insertion or deletion)
    pub const fn is_gap(&self) -> bool {
        !self.is_match()
    }
    
    /// Convert to legacy AlignState for backward compatibility
    pub fn to_legacy_align_state(&self) -> AlignState {
        match self.basic_state {
            BasicAlignState::Match => AlignState::Match,
            BasicAlignState::Insertion => {
                if self.piece_index <= 1 {
                    AlignState::Insertion
                } else {
                    AlignState::Insertion2
                }
            }
            BasicAlignState::Deletion => {
                if self.piece_index <= 1 {
                    AlignState::Deletion
                } else {
                    AlignState::Deletion2
                }
            }
        }
    }
}


/// A struct that represents a node in the alignment graph
///
/// It holds cursors to a node in the reference graph to which we are aligning,
/// and a query offset.
///
/// This struct does not store the alignment state (e.g., whether it's in Match,
/// insertion or deletion state). It is up to the specific alignment graph
/// implementations to retain that context.
///
/// ### See also
///
/// - [`AlignmentGraph`]
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct AlignmentGraphNode<N, O>
where
    N: Copy + Clone + Debug + Eq,
    O: Copy + Clone + Debug + Eq,
{
    node: N,
    offset: O,
}


impl<N, O> AlignmentGraphNode<N, O>
    where N: NodeIndexType,
          O: OffsetType,
{

    pub fn new(node: N, offset: O) -> Self {
        Self {
            node,
            offset,
        }
    }

    #[inline(always)]
    pub fn node(&self) -> N {
        self.node
    }

    #[inline(always)]
    pub fn offset(&self) -> O {
        self.offset
    }
}


pub trait AlignmentGraph {
    type CostModel: AlignmentCosts;

    fn get_costs(&self) -> &Self::CostModel;

    fn initial_states<G, O>(&self, ref_graph: &G) -> Vec<AlignmentGraphNode<G::NodeIndex, O>>
        where G: AlignableRefGraph,
              O: OffsetType;

    fn is_end<G, O>(&self, ref_graph: &G, seq: &[u8], node: &AlignmentGraphNode<G::NodeIndex, O>, aln_state: AlignState) -> bool
        where G: AlignableRefGraph,
              O: OffsetType;

    fn expand_all<V, G, O, F>(
        &self,
        visited_data: &mut V,
        ref_graph: &G,
        seq: &[u8],
        score: Score,
        node: &AlignmentGraphNode<G::NodeIndex, O>,
        state: AlignState,
        f: F
    )
        where V: AstarVisited<G::NodeIndex, O>,
              G: AlignableRefGraph,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<G::NodeIndex, O>, AlignState);

    fn expand_ref_graph_end<V, N, O, F>(&self, visited_data: &mut V, parent: &AlignmentGraphNode<N, O>, score: Score, f: F)
        where V: AstarVisited<N, O>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState);

    fn expand_query_end<V, N, O, F>(&self, visited_data: &mut V, parent: &AlignmentGraphNode<N, O>, child: N, score: Score, f: F)
        where V: AstarVisited<N, O>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState);

    fn expand_mismatch<V, N, O, F>(
        &self,
        visited_data: &mut V,
        parent: &AlignmentGraphNode<N, O>,
        child: &AlignmentGraphNode<N, O>,
        score: Score,
        f: F
    )
        where V: AstarVisited<N, O>,
              N: NodeIndexType,
              O: OffsetType,
              F: FnMut(u8, AlignmentGraphNode<N, O>, AlignState);
}
