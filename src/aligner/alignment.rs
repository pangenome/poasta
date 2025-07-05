use crate::graphs::{AlignableRefGraph, NodeIndexType};

#[derive(Clone, Debug)]
pub struct AlignedPair<N>
where
    N: NodeIndexType
{
    /// Represents the node rank in the graph
    pub rpos: Option<N>,

    /// Query sequence position
    pub qpos: Option<usize>
}

impl<N> AlignedPair<N>
where
    N: NodeIndexType
{
    pub fn new(rpos: Option<N>, qpos: Option<usize>) -> Self {
        Self { rpos, qpos }
    }

    pub fn is_aligned(&self) -> bool {
        matches!((self.rpos, self.qpos), (Some(_), Some(_)))
    }

    pub fn is_indel(&self) -> bool {
        !self.is_aligned()
    }
    
    pub fn is_deletion(&self) -> bool {
        self.rpos.is_none() && self.qpos.is_some()
    }
    
    pub fn is_insertion(&self) -> bool {
        self.rpos.is_some() && self.qpos.is_none()
    }
}

pub type Alignment<N> = Vec<AlignedPair<N>>;

pub fn print_alignment<G, N>(graph: &G, sequence: &[u8], aln: &Alignment<N>) -> String
where
    G: AlignableRefGraph<NodeIndex=N>,
    N: NodeIndexType
{
    let mut graph_chars = Vec::with_capacity(aln.len());
    let mut aln_chars = Vec::with_capacity(aln.len());
    let mut query_chars = Vec::with_capacity(aln.len());

    for pair in aln {
        if pair.is_aligned() {
            let node = graph.get_symbol_char(pair.rpos.unwrap());
            let qry = char::from(sequence[pair.qpos.unwrap()]);

            graph_chars.push(node);
            aln_chars.push(if node == qry { '|' } else { '·' });
            query_chars.push(qry);
        } else if let Some(nix) = pair.rpos {
            let node = graph.get_symbol_char(nix);
            graph_chars.push(node);
            aln_chars.push(' ');
            query_chars.push('-');
        } else if let Some(qpos) = pair.qpos {
            let qry = char::from(sequence[qpos]);
            graph_chars.push('-');
            aln_chars.push(' ');
            query_chars.push(qry);
        }
    }

    format!(
        "{}\n{}\n{}",
        graph_chars.into_iter().collect::<String>(),
        aln_chars.into_iter().collect::<String>(),
        query_chars.into_iter().collect::<String>(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphs::mock::{create_test_graph1, NIx};
    use petgraph::graph::NodeIndex;

    type N = NodeIndex<NIx>;

    #[test]
    fn test_aligned_pair_is_aligned() {
        let pair = AlignedPair::<N>::new(Some(N::new(1)), Some(2));
        assert!(pair.is_aligned());
        assert!(!pair.is_insertion());
        assert!(!pair.is_deletion());
    }

    #[test]
    fn test_aligned_pair_is_insertion() {
        let pair = AlignedPair::<N>::new(Some(N::new(1)), None);
        assert!(pair.is_insertion());
        assert!(!pair.is_aligned());
        assert!(!pair.is_deletion());
    }

    #[test]
    fn test_aligned_pair_is_deletion() {
        let pair = AlignedPair::<N>::new(None, Some(0));
        assert!(pair.is_deletion());
        assert!(!pair.is_aligned());
        assert!(!pair.is_insertion());
    }

    #[test]
    fn test_print_alignment() {
        let graph = create_test_graph1();
        let seq = b"AC";
        let aln: Alignment<N> = vec![
            AlignedPair::new(Some(N::new(0)), Some(0)),
            AlignedPair::new(Some(N::new(1)), None),
            AlignedPair::new(None, Some(1)),
        ];
        let output = print_alignment(&graph, seq, &aln);
        assert_eq!(output, "---\n·  \nA-C");
    }
}
