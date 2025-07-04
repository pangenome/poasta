use poasta::graphs::AlignableRefGraph;
use poasta::io::{load_graph_from_fasta_msa, save_graph, load_graph};
use poasta::graphs::poa::POAGraphWithIx;

fn graph_counts(graph: &POAGraphWithIx) -> (usize, usize) {
    match graph {
        POAGraphWithIx::U8(g) => (g.node_count(), g.edge_count()),
        POAGraphWithIx::U16(g) => (g.node_count(), g.edge_count()),
        POAGraphWithIx::U32(g) => (g.node_count(), g.edge_count()),
        POAGraphWithIx::USIZE(g) => (g.node_count(), g.edge_count()),
    }
}

#[test]
fn test_save_and_reload_consistency() {
    let path = "tests/small_test.input.fa";
    let graph = load_graph_from_fasta_msa(path).expect("load fasta msa");

    let counts_before = graph_counts(&graph);

    let mut data = Vec::new();
    save_graph(&graph, &mut data).expect("save graph");

    let loaded = load_graph(data.as_slice()).expect("load graph");
    let counts_after = graph_counts(&loaded);

    assert_eq!(counts_before, counts_after);
}

