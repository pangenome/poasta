use assert_cmd::prelude::*;
use predicates::prelude::*;
use std::process::Command;

#[test]
fn lasagna_requires_subcommand() {
    Command::cargo_bin("lasagna")
        .unwrap()
        .assert()
        .failure()
        .stderr(predicate::str::contains("No subcommand"));
}

#[test]
fn lasagna_align_gaf_output() {
    Command::cargo_bin("lasagna")
        .unwrap()
        .args(["align", "tests/test.gfa", "tests/small_test.fa"])
        .assert()
        .success()
        .stdout(predicate::str::contains("seq1"));
}
