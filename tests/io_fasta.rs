use poasta::graphs::poa::POAGraph;
use poasta::io::fasta::poa_graph_to_fasta;

#[test]
fn poa_graph_to_fasta_basic() {
    let mut graph: POAGraph<u16> = POAGraph::new();
    let seq1 = b"AC";
    graph
        .add_alignment_with_weights("seq1", seq1, None, &vec![1; seq1.len()])
        .unwrap();
    let seq2 = b"ACGT";
    graph
        .add_alignment_with_weights("seq2", seq2, None, &vec![1; seq2.len()])
        .unwrap();

    let mut out = Vec::new();
    poa_graph_to_fasta(&graph, &mut out).unwrap();
    let s = String::from_utf8(out).unwrap();
    let expected = ">seq1\n---AC\n>seq2\nACGT--\n";
    assert_eq!(s, expected);
}

#[test]
fn poa_graph_to_fasta_empty() {
    let mut graph: POAGraph<u16> = POAGraph::new();
    graph
        .add_alignment_with_weights("empty", b"", None, &[])
        .unwrap();
    let mut out = Vec::new();
    poa_graph_to_fasta(&graph, &mut out).unwrap();
    let s = String::from_utf8(out).unwrap();
    let expected = ">empty\n";
    assert_eq!(s, expected);
}
