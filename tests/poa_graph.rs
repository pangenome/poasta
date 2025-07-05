use petgraph::graph::DiGraph;
use petgraph::algo::toposort;
use poasta::graphs::poa::POAGraph;
use poasta::graphs::AlignableRefGraph;
use poasta::aligner::alignment::{Alignment, AlignedPair};

fn build_mock_graph1() -> DiGraph<i64, ()> {
    let mut nmap = Vec::new();
    let mut g = DiGraph::<i64, ()>::default();

    for i in 1..=9 {
        nmap.push(g.add_node(i));
    }

    let edges = [
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (3, 7),
        (7, 8),
        (8, 9),
    ];

    for (s, t) in edges.iter() {
        g.add_edge(nmap[*s - 1], nmap[*t - 1], ());
    }

    let end_node = g.add_node(10);
    for n in g.node_indices().collect::<Vec<_>>() {
        if n != end_node && g.neighbors(n).count() == 0 {
            g.add_edge(n, end_node, ());
        }
    }

    g
}

#[test]
fn test_new_graph() {
    let g: POAGraph<u32> = POAGraph::new();
    assert_eq!(g.node_count(), 0);
    assert_eq!(g.node_count_with_start_and_end(), 2);
    assert_ne!(g.start_node(), g.end_node());
}

#[test]
fn test_add_alignment_and_get_node_ranks() {
    let mut g: POAGraph<u32> = POAGraph::new();

    let seq1 = b"123456";
    let weights1 = vec![1usize; seq1.len()];
    g.add_alignment_with_weights("first", seq1, None, &weights1).unwrap();

    let node_three = g
        .all_nodes()
        .find(|n| g.get_symbol(*n) == b'3')
        .unwrap();

    let seq2 = b"3789";
    let weights2 = vec![1usize; seq2.len()];
    let aln: Alignment<_> = vec![
        AlignedPair::new(Some(node_three), Some(0)),
        AlignedPair::new(None, Some(1)),
        AlignedPair::new(None, Some(2)),
        AlignedPair::new(None, Some(3)),
    ];
    g.add_alignment_with_weights("branch", seq2, Some(&aln), &weights2).unwrap();

    assert_eq!(g.node_count(), 9);
    assert_eq!(g.node_count_with_start_and_end(), 11);

    let mock = build_mock_graph1();
    let topo = toposort(&mock, None).unwrap();
    let mut expected = std::collections::HashMap::new();
    for (rank, node) in topo.iter().enumerate() {
        let val = mock[*node];
        if val <= 9 {
            expected.insert((val as u8 + b'0') as u8, rank + 1);
        }
    }

    let ranks = g.get_node_ranks();
    let mut actual = std::collections::HashMap::new();
    for n in g.all_nodes() {
        let symbol = g.get_symbol(n);
        if symbol == b'#' || symbol == b'$' {
            continue;
        }
        actual.insert(symbol, ranks[n.index()]);
    }

    assert_eq!(actual, expected);
}
