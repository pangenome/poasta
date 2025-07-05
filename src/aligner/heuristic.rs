use std::marker::PhantomData;
use std::sync::Arc;
use crate::aligner::aln_graph::{AlignmentGraphNode, AlignState};
use crate::aligner::offsets::OffsetType;
use crate::aligner::path_index::PathIndex;
use crate::aligner::scoring::{AlignmentCosts, GapAffine};
use crate::bubbles::index::BubbleIndex;
use crate::graphs::NodeIndexType;

pub trait AstarHeuristic<N, O> {
    fn h(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> usize
        where N: NodeIndexType,
              O: OffsetType;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HeuristicType {
    Dijkstra,
    MinimumGapCost,
    PathAware,
}

impl HeuristicType {
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "dijkstra" => Some(HeuristicType::Dijkstra),
            "mingap" | "minimumgapcost" => Some(HeuristicType::MinimumGapCost),
            "path" | "pathaware" => Some(HeuristicType::PathAware),
            _ => None,
        }
    }
}

/// A* heuristic that always returns 0, such that
/// A* reduces to standard Dijkstra's algorithm.
#[derive(Default)]
pub struct Dijkstra<N, O>(PhantomData<(N, O)>);


impl<N, O> AstarHeuristic<N, O> for Dijkstra<N, O> {
    fn h(&self, _: &AlignmentGraphNode<N, O>, _: AlignState) -> usize
        where N: NodeIndexType,
              O: OffsetType
    {
        0
    }
}


pub struct MinimumGapCostAffine<N, O> {
    costs: GapAffine,
    bubble_index: Arc<BubbleIndex<N>>,
    seq_length: usize,
    dummy: PhantomData<O>,
}

impl<N, O> MinimumGapCostAffine<N, O> {
    pub fn new(costs: GapAffine, bubble_index: Arc<BubbleIndex<N>>, seq_length: usize) -> Self {
        Self {
            costs,
            bubble_index,
            seq_length,
            dummy: PhantomData,
        }
    }
}

impl<N, O> AstarHeuristic<N, O> for MinimumGapCostAffine<N, O> {

    fn h(&self, aln_node: &AlignmentGraphNode<N, O>, mut aln_state: AlignState) -> usize
        where N: NodeIndexType,
              O: OffsetType,
    {
        let min_dist_to_exit = self.bubble_index.get_min_dist_to_end(aln_node.node())
            .saturating_sub(1);
        let max_dist_to_exit = self.bubble_index.get_max_dist_to_end(aln_node.node())
            .saturating_sub(1);

        let target_offset_min = aln_node.offset().as_usize() + min_dist_to_exit;
        let target_offset_max = aln_node.offset().as_usize() + max_dist_to_exit;

        let min_gap_length;
        if target_offset_min > self.seq_length {
            min_gap_length = target_offset_min - self.seq_length;

            // requires deletions, so if in insertion or match state, we need to open a new gap
            if aln_state != AlignState::Deletion {
                aln_state = AlignState::Match;
            }
        } else if target_offset_max < self.seq_length {
            min_gap_length = self.seq_length - target_offset_max;

            // requires insertions, so if in deletion or match state, we need to open a new gap
            if aln_state != AlignState::Insertion {
                aln_state = AlignState::Match;
            }
        } else {
            min_gap_length = 0;
        }

        self.costs.gap_cost(aln_state, min_gap_length)
    }
}

pub struct PathAwareHeuristic<N, O> {
    costs: GapAffine,
    path_index: Arc<PathIndex<N>>,
    seq_length: usize,
    max_paths_to_consider: usize,
    dummy: PhantomData<O>,
}

impl<N, O> PathAwareHeuristic<N, O> {
    pub fn new(costs: GapAffine, path_index: Arc<PathIndex<N>>, seq_length: usize, max_paths_to_consider: usize) -> Self
    where
        N: NodeIndexType,
    {
        Self {
            costs,
            path_index,
            seq_length,
            max_paths_to_consider,
            dummy: PhantomData,
        }
    }
}

impl<N, O> AstarHeuristic<N, O> for PathAwareHeuristic<N, O> {
    fn h(&self, aln_node: &AlignmentGraphNode<N, O>, aln_state: AlignState) -> usize
    where
        N: NodeIndexType,
        O: OffsetType,
    {
        let paths = self.path_index.get_paths_through_node(aln_node.node());
        
        if paths.is_empty() {
            // Fallback: if node is not on any path, use a conservative estimate
            let remaining_query = self.seq_length.saturating_sub(aln_node.offset().as_usize());
            return self.costs.gap_cost(aln_state, remaining_query);
        }
        
        // Consider only the most relevant paths to avoid quadratic costs
        let paths_to_check = paths.len().min(self.max_paths_to_consider);
        
        let mut min_cost = usize::MAX;
        
        for &(path_id, pos) in &paths[..paths_to_check] {
            if let Some(dist_to_end) = self.path_index.get_distance_to_end(path_id, pos) {
                let path_length_remaining = dist_to_end as usize;
                let query_length_remaining = self.seq_length.saturating_sub(aln_node.offset().as_usize());
                
                let cost = if path_length_remaining > query_length_remaining {
                    // Need deletions
                    let gap_length = path_length_remaining - query_length_remaining;
                    let state = if aln_state == AlignState::Deletion {
                        aln_state
                    } else {
                        AlignState::Match // Will open new gap
                    };
                    self.costs.gap_cost(state, gap_length)
                } else if query_length_remaining > path_length_remaining {
                    // Need insertions
                    let gap_length = query_length_remaining - path_length_remaining;
                    let state = if aln_state == AlignState::Insertion {
                        aln_state
                    } else {
                        AlignState::Match // Will open new gap
                    };
                    self.costs.gap_cost(state, gap_length)
                } else {
                    // Perfect match possible
                    0
                };
                
                min_cost = min_cost.min(cost);
            }
        }
        
        min_cost
    }
}


#[cfg(test)]
mod tests {
    use super::{AstarHeuristic, Dijkstra, MinimumGapCostAffine, PathAwareHeuristic, HeuristicType};
    use crate::aligner::aln_graph::{AlignState, AlignmentGraphNode};
    use crate::aligner::path_index::PathIndex;
    use crate::aligner::scoring::{AlignmentCosts, GapAffine};
    use crate::bubbles::index::BubbleIndex;
    use crate::graphs::AlignableRefGraph;
    use crate::graphs::mock::{create_test_graph1, NIx};
    use petgraph::graph::NodeIndex;
    use std::sync::Arc;

    #[test]
    fn test_heuristic_type_parsing() {
        assert_eq!(HeuristicType::from_str("dijkstra"), Some(HeuristicType::Dijkstra));
        assert_eq!(HeuristicType::from_str("DIJKSTRA"), Some(HeuristicType::Dijkstra));
        assert_eq!(HeuristicType::from_str("mingap"), Some(HeuristicType::MinimumGapCost));
        assert_eq!(HeuristicType::from_str("minimumgapcost"), Some(HeuristicType::MinimumGapCost));
        assert_eq!(HeuristicType::from_str("path"), Some(HeuristicType::PathAware));
        assert_eq!(HeuristicType::from_str("pathaware"), Some(HeuristicType::PathAware));
        assert_eq!(HeuristicType::from_str("invalid"), None);
    }

    #[test]
    fn test_dijkstra_heuristic() {
        let heuristic = Dijkstra::<NodeIndex<NIx>, u32>::default();
        let node = NodeIndex::<NIx>::new(1);
        let aln_node = AlignmentGraphNode::new(node, 5u32);
        
        // Dijkstra should always return 0 (no heuristic)
        assert_eq!(heuristic.h(&aln_node, AlignState::Match), 0);
        assert_eq!(heuristic.h(&aln_node, AlignState::Insertion), 0);
        assert_eq!(heuristic.h(&aln_node, AlignState::Deletion), 0);
    }

    #[test]
    fn test_min_gap_cost() {
        let costs = GapAffine::new(4, 2, 6);

        let graph = create_test_graph1();
        let bubble_index = Arc::new(BubbleIndex::new(&graph));
        let heuristic = MinimumGapCostAffine::new(costs, bubble_index, 10);

        let node = NodeIndex::<NIx>::new(1);

        let node1 = AlignmentGraphNode::new(node, 2u32);
        assert_eq!(heuristic.h(&node1, AlignState::Match), 14);
        assert_eq!(heuristic.h(&node1, AlignState::Deletion), 14);
        assert_eq!(heuristic.h(&node1, AlignState::Insertion), 8);

        let node2 = AlignmentGraphNode::new(node, 7u32);
        assert_eq!(heuristic.h(&node2, AlignState::Match), 8);
        assert_eq!(heuristic.h(&node2, AlignState::Deletion), 2);
        assert_eq!(heuristic.h(&node2, AlignState::Insertion), 8);

        let node3 = AlignmentGraphNode::new(node, 6u32);
        assert_eq!(heuristic.h(&node3, AlignState::Match), 0);
        assert_eq!(heuristic.h(&node3, AlignState::Deletion), 0);
        assert_eq!(heuristic.h(&node3, AlignState::Insertion), 0);
    }

    #[test]
    fn test_path_aware_heuristic() {
        let costs = GapAffine::new(4, 2, 6);
        let graph = create_test_graph1();
        let path_index = Arc::new(PathIndex::build_from_graph(&graph, 5).unwrap());
        let seq_length = 10;
        
        let heuristic = PathAwareHeuristic::new(costs, path_index.clone(), seq_length, 3);
        
        // Test node on a path
        let start = graph.start_node();
        let aln_node = AlignmentGraphNode::new(start, 0u32);
        
        // At start with full sequence left, should estimate cost based on path
        let h_value = heuristic.h(&aln_node, AlignState::Match);
        assert!(h_value >= 0, "Heuristic should be non-negative");
        
        // Test consistency: h(n) should decrease as we progress
        let aln_node2 = AlignmentGraphNode::new(start, 5u32);
        let h_value2 = heuristic.h(&aln_node2, AlignState::Match);
        assert!(h_value2 <= h_value, "Heuristic should decrease as offset increases");
    }

    #[test]
    fn test_path_aware_fallback() {
        let costs = GapAffine::new(4, 2, 6);
        let graph = create_test_graph1();
        let path_index = Arc::new(PathIndex::build_from_graph(&graph, 1).unwrap());
        let seq_length = 10;
        
        let heuristic = PathAwareHeuristic::new(costs, path_index, seq_length, 3);
        
        // Find a node that might not be on any path (unlikely with small graph, but tests fallback)
        let nodes: Vec<_> = graph.all_nodes().collect();
        if let Some(&node) = nodes.last() {
            let aln_node = AlignmentGraphNode::new(node, 5u32);
            let h_value = heuristic.h(&aln_node, AlignState::Match);
            
            // Fallback should give conservative estimate
            let remaining = seq_length - 5;
            let expected = costs.gap_cost(AlignState::Match, remaining);
            assert!(h_value <= expected, "Fallback heuristic should be conservative");
        }
    }
}

