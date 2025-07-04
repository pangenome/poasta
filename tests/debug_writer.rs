use poasta::debug::{DebugOutputWriter, messages::DebugOutputMessage};
use tempfile::tempdir;

#[test]
fn debug_writer_creates_files() {
    let dir = tempdir().unwrap();
    let writer = DebugOutputWriter::init(dir.path());

    writer.log(DebugOutputMessage::IntermediateGraph { graph_dot: "digraph {}".to_string() });
    writer.log(DebugOutputMessage::AstarData { visited_tsv: "a\tb\n".to_string() });
    writer.log(DebugOutputMessage::Terminate);

    writer.join().unwrap();

    let dot_path = dir.path().join("graph_for_none.dot");
    assert!(dot_path.exists(), "missing dot file");

    let tsv_path = dir.path().join("astar_iterations/none.iter0.tsv");
    assert!(tsv_path.exists(), "missing tsv file");
}
