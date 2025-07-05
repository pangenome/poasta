use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;
use tempfile::tempdir;

#[test]
fn poasta_requires_subcommand() {
    Command::cargo_bin("poasta")
        .unwrap()
        .assert()
        .failure()
        .stderr(predicate::str::contains("No subcommand"));
}

#[test]
fn poasta_align_fasta_stdout() {
    Command::cargo_bin("poasta")
        .unwrap()
        .args(["align", "-O", "fasta", "tests/small_test.fa"])
        .assert()
        .success()
        .stdout(predicate::str::contains(">seq1"));
}

#[test]
fn poasta_align_and_view_gfa() {
    let dir = tempdir().unwrap();
    let out = dir.path().join("graph.poasta");

    Command::cargo_bin("poasta")
        .unwrap()
        .args(["align", "-o", out.to_str().unwrap(), "tests/small_test.fa"])
        .assert()
        .success();

    Command::cargo_bin("poasta")
        .unwrap()
        .args(["view", out.to_str().unwrap(), "-O", "gfa"])
        .assert()
        .success()
        .stdout(predicate::str::contains("S\t"));
}
