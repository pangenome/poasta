use std::fmt::{self, Display, Write};

use itertools::Itertools;

use crate::aligner::Alignment;
use crate::graphs::poa::{POAGraph, POANodeIndex};
use crate::graphs::AlignableRefGraph;

use super::gfa::{Field, FieldValue};
use super::graph::GraphSegments;

/// Resolve which GFA segment a POA node belongs to.
#[derive(Copy, Clone)]
pub struct NodeSegmentResolver<'a, Ix>
where
    Ix: petgraph::graph::IndexType + serde::de::DeserializeOwned,
{
    graph: &'a POAGraph<Ix>,
    segments: &'a GraphSegments<Ix>,
}

impl<'a, Ix> NodeSegmentResolver<'a, Ix>
where
    Ix: petgraph::graph::IndexType + serde::de::DeserializeOwned,
{
    /// Create a new resolver for the given graph.
    pub fn new(graph: &'a POAGraph<Ix>, segments: &'a GraphSegments<Ix>) -> Self {
        Self { graph, segments }
    }

    /// Return the segment index and position of `node` within that segment.
    pub fn resolve(&self, node: POANodeIndex<Ix>) -> Option<(usize, usize)> {
        for (segment_ix, (&start, &end)) in
            self.segments.start_nodes.iter().zip(&self.segments.end_nodes).enumerate()
        {
            let mut curr = start;
            let mut pos = 0usize;

            loop {
                if curr == node {
                    return Some((segment_ix, pos));
                }

                if curr == end {
                    break;
                }

                curr = self.graph.successors(curr).next()?;
                pos += 1;
            }
        }

        None
    }
}


pub struct GAFRecord {
    pub query_name: String,
    pub query_length: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: char,
    pub graph_path: String,
    pub path_length: usize,
    pub path_aln_start: usize,
    pub path_aln_end: usize,
    pub num_matches: usize,
    pub aln_block_len: usize,
    pub mapping_quality: usize,
    pub additional_fields: Vec<Field>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graphs::poa::POAGraph;
    use crate::graphs::AlignableRefGraph;

    #[test]
    fn resolve_returns_segment_index_and_position() {
        let mut graph = POAGraph::<u32>::new();

        let weights = vec![1; 2];
        let (s1_start, s1_end) = graph
            .add_nodes_for_sequence(b"AC", &weights, 0, 2)
            .unwrap();
        let (s2_start, s2_end) = graph
            .add_nodes_for_sequence(b"GT", &weights, 0, 2)
            .unwrap();

        graph.add_edge(s1_end, s2_start, 0, 1);

        let mut segments = GraphSegments::default();
        segments.names.push("s1".to_string());
        segments.start_nodes.push(s1_start);
        segments.end_nodes.push(s1_end);
        segments.segment_lengths.push(2);

        segments.names.push("s2".to_string());
        segments.start_nodes.push(s2_start);
        segments.end_nodes.push(s2_end);
        segments.segment_lengths.push(2);

        let resolver = NodeSegmentResolver::new(&graph, &segments);

        let s1_second = graph.successors(s1_start).next().unwrap();
        assert_eq!(resolver.resolve(s1_start), Some((0, 0)));
        assert_eq!(resolver.resolve(s1_second), Some((0, 1)));

        let s2_second = graph.successors(s2_start).next().unwrap();
        assert_eq!(resolver.resolve(s2_start), Some((1, 0)));
        assert_eq!(resolver.resolve(s2_second), Some((1, 1)));

        assert_eq!(resolver.resolve(graph.start_node()), None);
    }
}

impl Display for GAFRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let fields_str = self.additional_fields.iter()
            .fold(String::new(), |mut output, f| {
                let _ = write!(output, "\t{}", f);
                
                output
            });
        
        write!(
            f, 
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.query_length,
            self.query_start,
            self.query_end,
            self.strand,
            self.graph_path,
            self.path_length,
            self.path_aln_start,
            self.path_aln_end,
            self.num_matches,
            self.aln_block_len,
            self.mapping_quality,
            fields_str.trim()
        )?;
        
        Ok(())
    }
}
            


pub fn alignment_to_gaf<Ix>(
    graph: &POAGraph<Ix>,
    graph_segments: &GraphSegments<Ix>,
    seq_name: &str,
    sequence: &[u8],
    alignment: &Alignment<POANodeIndex<Ix>>,
    resolver: &NodeSegmentResolver<Ix>,
) -> Option<GAFRecord>
where 
    Ix: petgraph::graph::IndexType + serde::de::DeserializeOwned,
{
    if alignment.is_empty() {
        return None;
    }
    
    let mut query_start = 0;
    let mut path_aln_start = 0;
    let mut path_segments = Vec::new();
    let mut cigar_ops = Vec::new();
    
    let mut at_aln_start = true;
    
    let mut last_match_segment_ix = 0;
    let mut last_match_segment_pos = 0;
    let mut num_matches = 0;
    for aln_pair in alignment.iter() {
        if at_aln_start {
            if aln_pair.is_insertion() {
                query_start += 1;
            } else if aln_pair.is_aligned() {
                let (segment_ix, segment_pos) = resolver
                    .resolve(aln_pair.rpos.unwrap())
                    .expect("node not found in any segment");
                path_aln_start = segment_pos;
                path_segments.push(segment_ix);
                cigar_ops.push(if graph.is_symbol_equal(aln_pair.rpos.unwrap(), sequence[aln_pair.qpos.unwrap()]) {
                    num_matches += 1;
                    '='
                } else {
                    'X'
                });
                    
                at_aln_start = false;
                last_match_segment_ix = path_segments.len() - 1;
                last_match_segment_pos = segment_pos;
            }
        } else {
            match (aln_pair.rpos, aln_pair.qpos) {
                (Some(node), Some(qpos)) => {
                    let (segment_ix, segment_pos) = resolver
                        .resolve(node)
                        .expect("node not found in any segment");

                    if path_segments.last().copied() != Some(segment_ix) {
                        path_segments.push(segment_ix);
                    }
                    
                    cigar_ops.push(if graph.is_symbol_equal(node, sequence[qpos]) {
                        num_matches += 1;
                        '='
                    } else {
                        'X'
                    });
                    
                    last_match_segment_ix = path_segments.len() - 1;
                    last_match_segment_pos = segment_pos;
                },

                (Some(node), None) => {
                    let (segment_ix, _) = resolver
                        .resolve(node)
                        .expect("node not found in any segment");

                    if path_segments.last().copied() != Some(segment_ix) {
                        path_segments.push(segment_ix);
                    }
                    
                    cigar_ops.push('D');
                },
                
                (None, Some(_)) => {
                    cigar_ops.push('I');
                },
                _ => unreachable!(),
            }
        }
    }
    
    let graph_path = path_segments[..=last_match_segment_ix].iter()
        .fold(String::new(), |mut output, segment_ix| {
            let _ = write!(output, ">{}", graph_segments.names[*segment_ix]);
            output
        });
    
    let path_length = path_segments[..=last_match_segment_ix].iter()
        .map(|segment_ix| graph_segments.segment_lengths[*segment_ix])
        .sum();
    
    let path_aln_end = path_length 
        - graph_segments.segment_lengths[path_segments[last_match_segment_ix]] 
        + last_match_segment_pos;
    
    let query_end = alignment.iter().rev()
        .find(|aln_pair| aln_pair.is_aligned())
        .unwrap()
        .qpos.unwrap();
    
    let mut cigar_rle = cigar_ops.iter()
        .chunk_by(|op| **op)
        .into_iter()
        .map(|(op, group)| (op, group.count()))
        .collect::<Vec<_>>();
    
    match cigar_rle.last() {
        Some(('I', _)) => {
            let removed = cigar_rle.pop().unwrap();
            eprintln!("Removed insertion from end {}", removed.1);
        },
        Some(('D', _)) => {
            let removed = cigar_rle.pop().unwrap();
            eprintln!("Removed insertion from end {}", removed.1);
        }
        _ => (),
    }
    
    let aln_block_len = cigar_rle.iter()
        .map(|(_, count)| *count)
        .sum();
    
    let cigar_string = cigar_rle.iter()
        .fold(String::new(), |mut output, (op, count)| {
            let _ = write!(output, "{}{}", *count, *op);
            output
        });
    
    Some(GAFRecord {
        query_name: seq_name.to_string(),
        query_length: sequence.len(),
        query_start,
        query_end,
        strand: '+',
        graph_path,
        path_length,
        path_aln_start,
        path_aln_end,
        num_matches,
        aln_block_len,
        mapping_quality: 60,
        additional_fields: vec![
            Field { tag: "cg".to_string(), value: FieldValue::String(cigar_string) },
        ],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::graph::{load_graph_from_gfa, POAGraphFromGFA};
    use crate::graphs::AlignableRefGraph;
    use std::path::Path;

    #[test]
    fn resolves_nodes_to_segments() {
        let gfa_path = Path::new("tests/test.gfa");
        let POAGraphFromGFA { graph, graph_segments } =
            load_graph_from_gfa::<u32>(gfa_path).unwrap();

        let resolver = NodeSegmentResolver::new(&graph, &graph_segments);

        // segment 0 (s1)
        let s1_start = graph_segments.start_nodes[0];
        let s1_second = graph.successors(s1_start).find(|&n| n != graph.end_node()).unwrap();
        let s1_end = graph_segments.end_nodes[0];
        assert_eq!(resolver.resolve(s1_start), Some((0, 0)));
        assert_eq!(resolver.resolve(s1_second), Some((0, 1)));
        assert_eq!(
            resolver.resolve(s1_end),
            Some((0, graph_segments.segment_lengths[0] - 1))
        );

        // segment 1 (s2)
        let s2_start = graph_segments.start_nodes[1];
        let s2_end = graph_segments.end_nodes[1];
        let s2_second = graph.successors(s2_start).find(|&n| n != graph.end_node()).unwrap();
        assert_eq!(resolver.resolve(s2_start), Some((1, 0)));
        assert_eq!(resolver.resolve(s2_second), Some((1, 1)));
        assert_eq!(
            resolver.resolve(s2_end),
            Some((1, graph_segments.segment_lengths[1] - 1))
        );

        // segment 3 (s4)
        let s4_start = graph_segments.start_nodes[3];
        let s4_end = graph_segments.end_nodes[3];
        assert_eq!(resolver.resolve(s4_start), Some((3, 0)));
        assert_eq!(
            resolver.resolve(s4_end),
            Some((3, graph_segments.segment_lengths[3] - 1))
        );

        // nodes outside of segments
        assert_eq!(resolver.resolve(graph.start_node()), None);
        assert_eq!(resolver.resolve(graph.end_node()), None);
    }
}
