use poasta::aligner::config::{AffineDijkstra, AffineMinGapCost, AffinePathAware};
use poasta::aligner::scoring::{AlignmentType, GapAffine};
use poasta::aligner::PoastaAligner;
use poasta::graphs::poa::POAGraph;

#[test]
fn test_heuristic_consistency() {
    // Create a simple graph with a reference sequence
    let mut graph = POAGraph::<u32>::new();
    let reference = b"ACGTACGTACGTACGTACGT";
    graph.add_alignment_with_weights("reference", reference, None, &vec![1; reference.len()]).unwrap();
    
    // Add some variations to make the graph more complex
    let variant1 = b"ACGTACGTAAGTACGTACGT"; // substitution
    graph.add_alignment_with_weights("variant1", variant1, None, &vec![1; variant1.len()]).unwrap();
    
    let variant2 = b"ACGTACGTACGTACGT"; // deletion
    graph.add_alignment_with_weights("variant2", variant2, None, &vec![1; variant2.len()]).unwrap();
    
    // Query sequence with multiple possible alignments
    let query = b"ACGTACGTAAGTACGTACGT";
    
    let costs = GapAffine::new(4, 2, 6);
    
    // Test with Dijkstra (no heuristic)
    let dijkstra_aligner = PoastaAligner::new(
        AffineDijkstra(costs),
        AlignmentType::Global
    );
    let dijkstra_result = dijkstra_aligner.align::<u32, _>(&graph, query);
    
    // Test with MinimumGapCost
    let mingap_aligner = PoastaAligner::new(
        AffineMinGapCost(costs),
        AlignmentType::Global
    );
    let mingap_result = mingap_aligner.align::<u32, _>(&graph, query);
    
    // Test with PathAware
    let path_aligner = PoastaAligner::new(
        AffinePathAware(costs),
        AlignmentType::Global
    );
    let path_result = path_aligner.align::<u32, _>(&graph, query);
    
    // All heuristics should produce the same optimal alignment score
    assert_eq!(dijkstra_result.score, mingap_result.score, 
        "Dijkstra and MinGap should produce same score");
    assert_eq!(dijkstra_result.score, path_result.score,
        "Dijkstra and PathAware should produce same score");
    
    // Heuristics should reduce the number of states visited
    assert!(mingap_result.num_visited <= dijkstra_result.num_visited,
        "MinGap should visit fewer or equal states than Dijkstra");
    assert!(path_result.num_visited <= dijkstra_result.num_visited,
        "PathAware should visit fewer or equal states than Dijkstra");
    
    println!("Dijkstra visited: {}, MinGap visited: {}, PathAware visited: {}",
        dijkstra_result.num_visited, mingap_result.num_visited, path_result.num_visited);
}

#[test]
fn test_path_aware_complex_graph() {
    // Create a more complex graph with multiple long paths
    let mut graph = POAGraph::<u32>::new();
    
    // Add a long reference
    let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    graph.add_alignment_with_weights("reference", reference, None, &vec![1; reference.len()]).unwrap();
    
    // Add multiple variants with different patterns
    let deletion_variant = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    graph.add_alignment_with_weights("deletion", deletion_variant, None, &vec![1; deletion_variant.len()]).unwrap();
    
    let insertion_variant = b"ACGTACGTACGTACGTAAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    graph.add_alignment_with_weights("insertion", insertion_variant, None, &vec![1; insertion_variant.len()]).unwrap();
    
    let complex_variant = b"ACGTACGTAAGTACCTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    graph.add_alignment_with_weights("complex", complex_variant, None, &vec![1; complex_variant.len()]).unwrap();
    
    // Query that combines features from different paths
    let query = b"ACGTACGTAAGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    
    let costs = GapAffine::new(4, 2, 6);
    
    // Compare PathAware with MinGap on complex graph
    let mingap_aligner = PoastaAligner::new(
        AffineMinGapCost(costs),
        AlignmentType::Global
    );
    let mingap_result = mingap_aligner.align::<u32, _>(&graph, query);
    
    let path_aligner = PoastaAligner::new(
        AffinePathAware(costs),
        AlignmentType::Global
    );
    let path_result = path_aligner.align::<u32, _>(&graph, query);
    
    // Should produce same alignment
    assert_eq!(mingap_result.score, path_result.score,
        "MinGap and PathAware should produce same score on complex graph");
    
    // PathAware might be more efficient on graphs with clear path structure
    println!("Complex graph - MinGap visited: {}, PathAware visited: {}",
        mingap_result.num_visited, path_result.num_visited);
}

#[test]
fn test_heuristic_global_alignment() {
    let mut graph = POAGraph::<u32>::new();
    let reference = b"ACGTACGTACGTACGT";
    graph.add_alignment_with_weights("reference", reference, None, &vec![1; reference.len()]).unwrap();
    
    // Add variant to make graph more complex
    let variant = b"ACGTAAGTACGTACGT";
    graph.add_alignment_with_weights("variant", variant, None, &vec![1; variant.len()]).unwrap();
    
    // Query that requires alignment choices
    let query = b"ACGTACGTAAGTACGT";
    
    let costs = GapAffine::new(4, 2, 6);
    
    // Test all three heuristics with global alignment
    let dijkstra_aligner = PoastaAligner::new(AffineDijkstra(costs), AlignmentType::Global);
    let mingap_aligner = PoastaAligner::new(AffineMinGapCost(costs), AlignmentType::Global);
    let path_aligner = PoastaAligner::new(AffinePathAware(costs), AlignmentType::Global);
    
    let dijkstra_result = dijkstra_aligner.align::<u32, _>(&graph, query);
    let mingap_result = mingap_aligner.align::<u32, _>(&graph, query);
    let path_result = path_aligner.align::<u32, _>(&graph, query);
    
    // For global alignment, all heuristics MUST produce the same optimal score
    assert_eq!(dijkstra_result.score, mingap_result.score,
        "Dijkstra and MinGap must produce same score in global alignment");
    assert_eq!(dijkstra_result.score, path_result.score,
        "Dijkstra and PathAware must produce same score in global alignment");
    
    // Heuristics should reduce search space
    println!("Global alignment - Dijkstra: {} states, MinGap: {} states, PathAware: {} states",
        dijkstra_result.num_visited, mingap_result.num_visited, path_result.num_visited);
    
    assert!(mingap_result.num_visited <= dijkstra_result.num_visited,
        "MinGap should not visit more states than Dijkstra");
    
    // PathAware might visit more states if paths don't align well with optimal solution
    // but it should still find the optimal score
    println!("PathAware visited {}% of Dijkstra's states", 
        (path_result.num_visited * 100) / dijkstra_result.num_visited);
}