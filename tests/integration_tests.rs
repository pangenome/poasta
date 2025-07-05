use std::process::Command;
use std::fs;
use std::path::Path;

/// Test that the ends-free CLI command works
#[test]
fn test_ends_free_cli_basic() {
    // Create test data
    let test_file = "tests/test_ends_free_cli.fa";
    let test_content = ">ref\nATCGATCG\n>query\nCGATC\n";
    fs::write(test_file, test_content).expect("Failed to write test file");

    // Run poasta with ends-free alignment
    let output = Command::new("cargo")
        .args(&["run", "--release", "--bin", "poasta", "--", 
                "align", test_file, "-m", "ends-free", "-g", "8,24", "-e", "2,1", "-O", "fasta"])
        .output()
        .expect("Failed to execute command");

    // Should succeed
    assert!(output.status.success(), 
            "Command failed: {}", String::from_utf8_lossy(&output.stderr));
    
    // Should produce valid output
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains(">ref"), "Output should contain reference sequence");
    assert!(stdout.contains(">query"), "Output should contain query sequence");

    // Clean up
    let _ = fs::remove_file(test_file);
}

/// Test standard gap affine vs two-piece gap affine with ends-free
#[test]
fn test_gap_models_comparison() {
    let test_file = "tests/test_gap_models.fa";
    let test_content = ">ref\nAAAATTTTGGGGCCCC\n>query\nTTTTGGGG\n";
    fs::write(test_file, test_content).expect("Failed to write test file");

    // Test standard gap affine
    let output_standard = Command::new("cargo")
        .args(&["run", "--release", "--bin", "poasta", "--", 
                "align", test_file, "-m", "ends-free", "-g", "8", "-e", "2", "-O", "fasta"])
        .output()
        .expect("Failed to execute standard affine command");

    assert!(output_standard.status.success(), 
            "Standard affine command failed: {}", String::from_utf8_lossy(&output_standard.stderr));

    // Test two-piece gap affine
    let output_2piece = Command::new("cargo")
        .args(&["run", "--release", "--bin", "poasta", "--", 
                "align", test_file, "-m", "ends-free", "-g", "8,24", "-e", "2,1", "-O", "fasta"])
        .output()
        .expect("Failed to execute two-piece affine command");

    assert!(output_2piece.status.success(), 
            "Two-piece affine command failed: {}", String::from_utf8_lossy(&output_2piece.stderr));

    // Both should produce valid alignments
    let stdout_standard = String::from_utf8_lossy(&output_standard.stdout);
    let stdout_2piece = String::from_utf8_lossy(&output_2piece.stdout);
    
    assert!(stdout_standard.contains(">ref"), "Standard output should contain reference");
    assert!(stdout_2piece.contains(">ref"), "Two-piece output should contain reference");

    // Clean up
    let _ = fs::remove_file(test_file);
}

/// Test ends-free vs global alignment behavior
#[test]
fn test_ends_free_vs_global() {
    let test_file = "tests/test_alignment_modes.fa";
    let test_content = ">ref\nATCGATCGATCG\n>query\nGATC\n";
    fs::write(test_file, test_content).expect("Failed to write test file");

    // Test global alignment (might fail or give poor score)
    let _output_global = Command::new("cargo")
        .args(&["run", "--release", "--bin", "poasta", "--", 
                "align", test_file, "-m", "global", "-g", "8,24", "-e", "2,1", "-O", "fasta"])
        .output()
        .expect("Failed to execute global command");

    // Test ends-free alignment (should work well)
    let output_ends_free = Command::new("cargo")
        .args(&["run", "--release", "--bin", "poasta", "--", 
                "align", test_file, "-m", "ends-free", "-g", "8,24", "-e", "2,1", "-O", "fasta"])
        .output()
        .expect("Failed to execute ends-free command");

    // Ends-free should definitely succeed
    assert!(output_ends_free.status.success(), 
            "Ends-free command failed: {}", String::from_utf8_lossy(&output_ends_free.stderr));

    let stdout_ends_free = String::from_utf8_lossy(&output_ends_free.stdout);
    assert!(stdout_ends_free.contains(">ref"), "Ends-free output should contain reference");
    assert!(stdout_ends_free.contains(">query"), "Ends-free output should contain query");

    // Clean up
    let _ = fs::remove_file(test_file);
}

/// Test edge cases
#[test]
fn test_edge_cases() {
    // Test with empty query
    let test_file_empty = "tests/test_empty_query.fa";
    let test_content_empty = ">ref\nATCG\n>query\n\n";
    fs::write(test_file_empty, test_content_empty).expect("Failed to write test file");

    let output_empty = Command::new("cargo")
        .args(&["run", "--release", "--bin", "poasta", "--", 
                "align", test_file_empty, "-m", "ends-free", "-g", "8,24", "-e", "2,1", "-O", "fasta"])
        .output()
        .expect("Failed to execute empty query command");

    assert!(output_empty.status.success(), 
            "Empty query command failed: {}", String::from_utf8_lossy(&output_empty.stderr));

    // Test with single nucleotide
    let test_file_single = "tests/test_single_nuc.fa";
    let test_content_single = ">ref\nATCGATCG\n>query\nC\n";
    fs::write(test_file_single, test_content_single).expect("Failed to write test file");

    let output_single = Command::new("cargo")
        .args(&["run", "--release", "--bin", "poasta", "--", 
                "align", test_file_single, "-m", "ends-free", "-g", "8,24", "-e", "2,1", "-O", "fasta"])
        .output()
        .expect("Failed to execute single nucleotide command");

    assert!(output_single.status.success(), 
            "Single nucleotide command failed: {}", String::from_utf8_lossy(&output_single.stderr));

    // Clean up
    let _ = fs::remove_file(test_file_empty);
    let _ = fs::remove_file(test_file_single);
}

/// Test multiple sequences
#[test]
fn test_multiple_sequences() {
    let test_file = "tests/test_multiple_seqs.fa";
    let test_content = ">ref\nATCGATCGATCG\n>query1\nTCGA\n>query2\nGATC\n>query3\nCGAT\n";
    fs::write(test_file, test_content).expect("Failed to write test file");

    let output = Command::new("cargo")
        .args(&["run", "--release", "--bin", "poasta", "--", 
                "align", test_file, "-m", "ends-free", "-g", "8,24", "-e", "2,1", "-O", "fasta"])
        .output()
        .expect("Failed to execute multiple sequences command");

    assert!(output.status.success(), 
            "Multiple sequences command failed: {}", String::from_utf8_lossy(&output.stderr));

    let stdout = String::from_utf8_lossy(&output.stdout);
    
    // Should contain all sequences
    assert!(stdout.contains(">ref"), "Output should contain reference");
    assert!(stdout.contains(">query1"), "Output should contain query1");
    assert!(stdout.contains(">query2"), "Output should contain query2");
    assert!(stdout.contains(">query3"), "Output should contain query3");
    
    // Count the number of sequences (each sequence header starts with >)
    let seq_count = stdout.matches('>').count();
    assert_eq!(seq_count, 4, "Should have 4 sequences in output, found {}", seq_count);

    // Clean up
    let _ = fs::remove_file(test_file);
}

/// Test parameter validation
#[test]
fn test_parameter_validation() {
    let test_file = "tests/test_params.fa";
    let test_content = ">ref\nATCG\n>query\nTCG\n";
    fs::write(test_file, test_content).expect("Failed to write test file");

    // Test invalid two-piece parameters (should fall back to standard affine)
    let output_invalid = Command::new("cargo")
        .args(&["run", "--release", "--bin", "poasta", "--", 
                "align", test_file, "-m", "ends-free", "-g", "8,24", "-e", "2,2", "-O", "fasta"])
        .output()
        .expect("Failed to execute invalid params command");

    // Should still work (with warning) by falling back to standard affine
    assert!(output_invalid.status.success(), 
            "Invalid params command failed: {}", String::from_utf8_lossy(&output_invalid.stderr));

    let stderr = String::from_utf8_lossy(&output_invalid.stderr);
    assert!(stderr.contains("Warning") || stderr.contains("Using standard affine"), 
            "Should contain warning about invalid parameters");

    // Clean up
    let _ = fs::remove_file(test_file);
}

/// Test with the original failing data if available
#[test]
fn test_original_data_if_available() {
    if Path::new("partition714.del.fa").exists() {
        let output = Command::new("cargo")
            .args(&["run", "--release", "--bin", "poasta", "--", 
                    "align", "partition714.del.fa", "-m", "ends-free", "-n", "5", 
                    "-g", "8,24", "-e", "2,1", "-O", "fasta"])
            .output()
            .expect("Failed to execute original data command");

        assert!(output.status.success(), 
                "Original data command failed: {}", String::from_utf8_lossy(&output.stderr));

        let stdout = String::from_utf8_lossy(&output.stdout);
        
        // Should produce multiple aligned sequences
        let seq_count = stdout.matches('>').count();
        assert!(seq_count > 1, "Should produce multiple sequences, found {}", seq_count);
    }
}

/// Test performance with timing
#[test]
fn test_performance_basic() {
    let test_file = "tests/test_performance.fa";
    let test_content = ">ref\nATCGATCGATCGATCGATCGATCGATCGATCG\n>query\nCGATCGATCGATCG\n";
    fs::write(test_file, test_content).expect("Failed to write test file");

    let start = std::time::Instant::now();
    
    let output = Command::new("cargo")
        .args(&["run", "--release", "--bin", "poasta", "--", 
                "align", test_file, "-m", "ends-free", "-g", "8,24", "-e", "2,1", "-O", "fasta"])
        .output()
        .expect("Failed to execute performance test command");

    let duration = start.elapsed();

    assert!(output.status.success(), 
            "Performance test command failed: {}", String::from_utf8_lossy(&output.stderr));

    // Should complete reasonably quickly (within 10 seconds for this small test)
    assert!(duration.as_secs() < 10, 
            "Performance test took too long: {:?}", duration);

    // Clean up
    let _ = fs::remove_file(test_file);
}