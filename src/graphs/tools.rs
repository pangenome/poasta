use std::cell::RefCell;
use rustc_hash::FxHashSet;
use crate::graphs::AlignableRefGraph;

pub fn rev_postorder_nodes<G: AlignableRefGraph>(graph: &G) -> Vec<G::NodeIndex> {
    let mut ordered = Vec::with_capacity(graph.node_count_with_start_and_end());

    let mut stack = vec![
        (graph.start_node(), RefCell::new(graph.successors(graph.start_node())))
    ];

    let next_valid_child = |succ_iter: &RefCell<G::SuccessorIterator<'_>>, visited: &FxHashSet<G::NodeIndex>| {
        while let Some(child) = succ_iter.borrow_mut().next() {
            if !visited.contains(&child) {
                return Some(child)
            }
        }

        None
    };

    let mut visited = FxHashSet::default();
    while !stack.is_empty() {
        let (_, succ_iter) = stack.last().unwrap();

        if let Some(child) = next_valid_child(succ_iter, &visited) {
            visited.insert(child);
            stack.push((child, RefCell::new(graph.successors(child))))
        } else {
            let last = stack.pop().unwrap();
            ordered.push(last.0);
        }
    }

    ordered.reverse();
    ordered
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::mock::MockGraph;

    #[test]
    fn single_node_graph() {
        let mut g = MockGraph::default();
        let n = g.add_node(1);

        let order = rev_postorder_nodes(&g);
        assert_eq!(order, vec![n]);
    }

    #[test]
    fn small_branching_graph() {
        let mut g = MockGraph::default();
        let n1 = g.add_node(1);
        let n2 = g.add_node(2);
        let n3 = g.add_node(3);
        let n4 = g.add_node(4);

        g.add_edge(n1, n2, ());
        g.add_edge(n1, n3, ());
        g.add_edge(n2, n4, ());
        g.add_edge(n3, n4, ());

        let order = rev_postorder_nodes(&g);

        assert_eq!(order, vec![n1, n2, n3, n4]);
    }
}
