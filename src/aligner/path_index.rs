use crate::errors::PoastaError;
use crate::graphs::AlignableRefGraph;
use crate::graphs::NodeIndexType;
use smallvec::SmallVec;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct PathCoordinate<N> {
    pub path_id: u32,
    pub path_position: u32,
    pub node: N,
    pub alt_paths: SmallVec<[(u32, u32); 2]>,
}

#[derive(Debug, Clone)]
pub struct Path<N> {
    pub id: u32,
    pub name: String,
    pub nodes: Vec<N>,
    pub length: usize,
}

#[derive(Debug, Clone)]
pub struct PathDistanceInfo {
    pub path_id: u32,
    pub forward_distances: Vec<u32>,
    pub backward_distances: Vec<u32>,
}

#[derive(Debug)]
pub struct PathIndex<N> {
    paths: Vec<Path<N>>,
    node_to_paths: HashMap<N, SmallVec<[(u32, u32); 4]>>,
    path_distances: Vec<PathDistanceInfo>,
    max_paths_per_node: usize,
}

impl<N> PathIndex<N> {
    pub fn new(max_paths_per_node: usize) -> Self {
        Self {
            paths: Vec::new(),
            node_to_paths: HashMap::new(),
            path_distances: Vec::new(),
            max_paths_per_node,
        }
    }

    pub fn build_from_graph<G>(graph: &G, max_paths_per_node: usize) -> Result<Self, PoastaError>
    where
        G: AlignableRefGraph<NodeIndex = N>,
        N: NodeIndexType,
    {
        let mut index = Self::new(max_paths_per_node);
        
        // Find major paths through the graph using a greedy approach
        index.extract_major_paths(graph)?;
        
        // Compute distances for each path
        index.compute_path_distances();
        
        Ok(index)
    }

    fn extract_major_paths<G>(&mut self, graph: &G) -> Result<(), PoastaError>
    where
        G: AlignableRefGraph<NodeIndex = N>,
        N: NodeIndexType,
    {
        let mut visited_edges: HashSet<(N, N)> = HashSet::new();
        let mut path_id = 0;
        
        // Start from nodes with no incoming edges (or few incoming edges)
        let mut start_nodes: Vec<N> = Vec::new();
        
        // Find start node
        let start = graph.start_node();
        start_nodes.push(start);
        
        // Also consider nodes with high out-degree as potential path starts
        for node in graph.all_nodes() {
            let in_degree = graph.in_degree(node);
            let out_degree = graph.out_degree(node);
            
            if in_degree == 0 || (out_degree > 2 && in_degree == 1) {
                start_nodes.push(node);
            }
        }
        
        // Extract paths using a greedy approach
        for start_node in start_nodes {
            if visited_edges.iter().any(|(from, _)| *from == start_node) {
                continue;
            }
            
            let path = self.extract_path_from(graph, start_node, &mut visited_edges, path_id)?;
            if path.nodes.len() > 1 {
                self.add_path(path);
                path_id += 1;
            }
        }
        
        // If we have few paths, try to find more from unvisited high-degree nodes
        if self.paths.len() < 10 {
            self.extract_secondary_paths(graph, &mut visited_edges, &mut path_id)?;
        }
        
        Ok(())
    }

    fn extract_path_from<G>(
        &self,
        graph: &G,
        start: N,
        visited_edges: &mut HashSet<(N, N)>,
        path_id: u32,
    ) -> Result<Path<N>, PoastaError>
    where
        G: AlignableRefGraph<NodeIndex = N>,
        N: NodeIndexType,
    {
        let mut nodes = vec![start];
        let mut current = start;
        let mut length = 0;
        
        // Follow the path greedily, preferring unvisited edges
        while current != graph.end_node() {
            let neighbors: Vec<N> = graph.successors(current).collect();
            
            if neighbors.is_empty() {
                break;
            }
            
            // Prefer unvisited edges, then edges with highest out-degree
            let next = neighbors
                .iter()
                .filter(|&&n| !visited_edges.contains(&(current, n)))
                .max_by_key(|&&n| graph.out_degree(n))
                .or_else(|| neighbors.first())
                .copied();
            
            if let Some(next_node) = next {
                visited_edges.insert((current, next_node));
                nodes.push(next_node);
                length += 1;
                current = next_node;
            } else {
                break;
            }
        }
        
        Ok(Path {
            id: path_id,
            name: format!("path_{}", path_id),
            nodes,
            length,
        })
    }

    fn extract_secondary_paths<G>(
        &mut self,
        graph: &G,
        visited_edges: &mut HashSet<(N, N)>,
        path_id: &mut u32,
    ) -> Result<(), PoastaError>
    where
        G: AlignableRefGraph<NodeIndex = N>,
        N: NodeIndexType,
    {
        // Find nodes with unvisited outgoing edges
        let mut candidate_starts: Vec<(N, usize)> = Vec::new();
        
        for node in graph.all_nodes() {
            let unvisited_out = graph
                .successors(node)
                .filter(|&n| !visited_edges.contains(&(node, n)))
                .count();
            
            if unvisited_out > 0 {
                candidate_starts.push((node, unvisited_out));
            }
        }
        
        // Sort by number of unvisited edges
        candidate_starts.sort_by_key(|(_, count)| std::cmp::Reverse(*count));
        
        for (start, _) in candidate_starts.into_iter().take(20) {
            let path = self.extract_path_from(graph, start, visited_edges, *path_id)?;
            if path.nodes.len() > 3 {
                self.add_path(path);
                *path_id += 1;
            }
        }
        
        Ok(())
    }

    fn add_path(&mut self, path: Path<N>) 
    where
        N: NodeIndexType,
    {
        let path_id = path.id;
        
        // Update node_to_paths mapping
        for (pos, &node) in path.nodes.iter().enumerate() {
            let entry = self.node_to_paths.entry(node).or_insert_with(SmallVec::new);
            
            // Only keep the most relevant paths per node
            if entry.len() < self.max_paths_per_node {
                entry.push((path_id, pos as u32));
            }
        }
        
        self.paths.push(path);
    }

    fn compute_path_distances(&mut self) {
        for path in &self.paths {
            let mut forward_distances = vec![0u32; path.nodes.len()];
            let mut backward_distances = vec![0u32; path.nodes.len()];
            
            // Forward distances
            for i in 1..path.nodes.len() {
                forward_distances[i] = forward_distances[i - 1] + 1;
            }
            
            // Backward distances
            for i in (0..path.nodes.len() - 1).rev() {
                backward_distances[i] = backward_distances[i + 1] + 1;
            }
            
            self.path_distances.push(PathDistanceInfo {
                path_id: path.id,
                forward_distances,
                backward_distances,
            });
        }
    }

    pub fn get_paths_through_node(&self, node: N) -> &[(u32, u32)] 
    where
        N: NodeIndexType,
    {
        self.node_to_paths
            .get(&node)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    pub fn get_path(&self, path_id: u32) -> Option<&Path<N>> {
        self.paths.iter().find(|p| p.id == path_id)
    }

    pub fn get_path_length(&self, path_id: u32) -> usize {
        self.paths
            .iter()
            .find(|p| p.id == path_id)
            .map(|p| p.length)
            .unwrap_or(0)
    }

    pub fn get_distance_to_end(&self, path_id: u32, position: u32) -> Option<u32> {
        self.path_distances
            .iter()
            .find(|d| d.path_id == path_id)
            .and_then(|d| d.backward_distances.get(position as usize))
            .copied()
    }

    pub fn get_distance_from_start(&self, path_id: u32, position: u32) -> Option<u32> {
        self.path_distances
            .iter()
            .find(|d| d.path_id == path_id)
            .and_then(|d| d.forward_distances.get(position as usize))
            .copied()
    }

    pub fn num_paths(&self) -> usize {
        self.paths.len()
    }

    pub fn paths(&self) -> &[Path<N>] {
        &self.paths
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphs::mock::{create_test_graph1, create_test_graph2};

    #[test]
    fn test_path_index_creation() {
        let graph = create_test_graph1();
        
        let index = PathIndex::build_from_graph(&graph, 4).unwrap();
        
        assert!(index.num_paths() > 0);
        
        // Check that start node is on at least one path
        let start = graph.start_node();
        let paths = index.get_paths_through_node(start);
        assert!(!paths.is_empty(), "Start node should be on at least one path");
    }

    #[test]
    fn test_path_distances() {
        let graph = create_test_graph1();
        let index = PathIndex::build_from_graph(&graph, 4).unwrap();
        
        // Get first path
        if let Some(path) = index.paths().first() {
            // Check distances are computed correctly
            let start_dist = index.get_distance_from_start(path.id, 0);
            assert_eq!(start_dist, Some(0));
            
            let last_pos = path.nodes.len() - 1;
            let end_dist = index.get_distance_to_end(path.id, last_pos as u32);
            assert_eq!(end_dist, Some(0));
        }
    }

    #[test]
    fn test_max_paths_per_node() {
        let graph = create_test_graph2(); // Larger graph
        let max_paths = 3;
        let index = PathIndex::build_from_graph(&graph, max_paths).unwrap();
        
        // Check that no node has more paths than the limit
        for node in graph.all_nodes() {
            let paths = index.get_paths_through_node(node);
            assert!(paths.len() <= max_paths, 
                "Node should not have more than {} paths", max_paths);
        }
    }
}