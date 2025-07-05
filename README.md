<h1 align="center">🍝 POASTA</h1>
<h2 align="center">Fast, gap-affine sequence-to-graph and partial order aligner</h2>

<p>&nbsp;</p>

[![CI](https://github.com/broadinstitute/poasta/actions/workflows/ci.yml/badge.svg)](https://github.com/broadinstitute/poasta/actions/workflows/ci.yml) [![codecov](https://codecov.io/gh/broadinstitute/poasta/branch/main/graph/badge.svg)](https://codecov.io/gh/broadinstitute/poasta)


POASTA is a fast and optimal partial order aligner that supports gap-affine alignment penalties. Inspired by 
a recent [algorithm for pairwise alignment](https://github.com/smarco/WFA2-lib), it can exploit exact matches
between the query and the graph, greatly speeding up the alignment process.

## Installation

### Pre-built binaries

Pre-built executables for Linux and macOS are available from the
[GitHub releases page](https://github.com/broadinstitute/poasta/releases).
Download the archive that matches your platform and unpack it:

```bash
tar -xzf poasta-vX.Y.Z-x86_64-unknown-linux-gnu.tar.gz
sudo mv poasta /usr/local/bin/
```
Replace `X.Y.Z` with the release number you wish to install. The binary can
then be used directly from the command line.

### Conda

POASTA is packaged for bioconda. You can install it using `conda` or
`mamba` as follows:

```bash
conda install -c conda-forge -c bioconda poasta
```
This will fetch a pre-built binary suitable for your environment and add it
to your active conda environment's path.

### Building POASTA from source

#### Rust compiler

POASTA is written in Rust, and to build and install it, you'll need a recent version of the Rust compiler. The 
minimum supported Rust version is 1.70.

1. Download and install `rustup`: https://rustup.rs/
2. Run `rustup update`

#### Building POASTA

1. Clone the repository. 
    
   ```bash
   git clone https://github.com/broadinstitute/poasta
   ```
 
2. Move into the directory. 

   ```bash
   cd poasta
   ```
 
3. Build using `cargo`. We enable a flag to ensure the compiler uses all features of your machine's CPU. 
   To maximize portability of the binary, however, remove the `RUSTFLAGS="..."` part.
    
   ```bash
   RUSTFLAGS="-C target-cpu=native" cargo build --release
   ```
 
4. The built `poasta` executable will be available in the directory `target/release/`


## Usage

### Creating an alignment from scratch

To create a multiple sequence alignment from scratch, simply give it a FASTA. The FASTA file can be compressed 
with gzip (filename should have a `.gz` extension).

```bash
poasta align -o graph.poasta sequences.fasta
```

This will output the graph to a binary file called `graph.poasta`. POASTA can reuse this file to later align
additional sequences to it.

For complex graphs or specific performance requirements, you can choose different heuristics:

```bash
# Use path-aware heuristic for pangenome graphs with clear reference paths
poasta align -H path -o graph.poasta sequences.fasta
```

### Re-using an earlier alignment

To align additional sequences to an earlier created partial order graph, specify the existing graph using the 
`-g` option.

```bash
poasta align -g graph.poasta -o graph_updated.poasta new_sequences.fasta
```

This will import the graph stored in `graph.poasta`, then align the additional sequences in `new_sequences.fasta` to
this graph, and outputs the updated graph to `graph_updated.poasta`.

### Importing an existing multiple sequence alignment in FASTA format

POASTA can import an existing multiple sequence alignment stored in columnar FASTA format (e.g., those 
created by other tools like `mafft` or `spoa`), create the equivalent partial order graph from the existing alignment,
and then align new sequences to it. To achieve this, specify the FASTA MSA with extension .fa, .fna, or .fasta with
the `-g` option (file is also allowed to be compressed with gzip if it has a `.gz` suffix).

```bash
poasta align -g msa.fasta -o graph_updated.poasta new_sequences.fasta
```

### Other output formats

The default output format is POASTA's binary file format storing the graph data structure.
This the recommended output because you can always convert this binary file
to other formats using `poasta view`.
If you don't need the binary file, however,
you can specify the output format with the `-O` option:

Other supported formats:

* DOT (GraphViz): Specify with `-O dot`
* FASTA MSA: Specify with `-O fasta`
* Graph GFA: Specify with `-O gfa`

For example, to visualize the graph directly with GraphViz:

```bash
poasta align -Odot sequences.fasta | dot -Tpng -o graph.png
```

Note that we did not specify an output file for `poasta align` (we did not use the `-o` option). If no output filename
is given, standard output will be used, so the output can be directly piped to `dot` to create the visualization.

### Using `poasta view` to convert between output formats

By default, POASTA stores the computed graph/MSA in its own binary file format.
To convert a previously computed MSA to other file formats, you can use `poasta view`.
The supported output formats are the same as described above, i.e.:

* DOT (GraphViz): Specify with `-O dot`
* FASTA MSA: Specify with `-O fasta`
* Graph GFA: Specify with `-O gfa`

Example:

```bash
# Convert to GFA
poasta view -Ogfa existing_msa.poasta > poa_graph.gfa

# Convert to FASTA MSA
poasta view -Ofasta existing_msa.poasta > poa_msa.fasta
```
### Inspecting POA graphs

Use `poasta stats` to print statistics for a stored graph. The command
reports fields such as `node_count`, `node_count_with_start`,
`edge_count`, `avg_in_degree`, and `avg_out_degree`.

```bash
poasta stats graph.poasta 2> stats.log
```
Example output:

```
node_count: 124
node_count_with_start: 126
edge_count: 210
avg_in_degree: 1.65
avg_out_degree: 1.65
```

The `contrib/poasta_tools` directory includes plotting tools.
Run `poasta_plot.py` to visualize dynamic programming matrices:

```bash
python contrib/poasta_tools/poasta_plot.py matrix.tsv graph.dot
```

Use `poasta_graphviz_region.py` to render a subgraph:

```bash
python contrib/poasta_tools/poasta_graphviz_region.py \
    graph.dot chr1:1000-2000
```

## Gap penalty configuration

POASTA supports both standard affine gap penalties and two-piece affine gap penalties for more flexible gap modeling.

### Standard affine gap penalties

By default, POASTA uses standard affine gap penalties with the following defaults:
- Gap opening penalty: 6
- Gap extension penalty: 2
- Mismatch penalty: 4

You can customize these penalties:

```bash
poasta align -g 8 -e 3 -n 5 sequences.fasta
```

### Two-piece affine gap penalties

For more sophisticated gap modeling, POASTA supports two-piece affine gap penalties, which use different penalties for short and long gaps. This can better model biological insertion/deletion patterns.

To enable two-piece affine mode, provide comma-separated values for both gap opening and extension penalties:

```bash
# Two-piece affine: short gaps (8,2) and long gaps (24,1)
poasta align -g 8,24 -e 2,1 -n 5 sequences.fasta
```

The first values (8,2) apply to short gaps, while the second values (24,1) apply to longer gaps. The extension penalty for short gaps should be higher than for long gaps to properly model the two-piece behavior.

### Alignment modes

POASTA supports different alignment modes:
- `global`: End-to-end alignment (default)
- `semi-global`: Query aligned globally, but graph gaps are free at ends  
- `ends-free`: Free gaps at sequence beginnings/ends for both query and graph

```bash
poasta align -m ends-free sequences.fasta
```

### A* Heuristics

POASTA uses the A* algorithm for optimal sequence-to-graph alignment. You can choose between different heuristics that affect the search strategy:

- `dijkstra`: No heuristic (equivalent to Dijkstra's algorithm). Explores all paths systematically.
- `mingap`: Minimum gap cost heuristic using bubble structure (default). Provides good performance on most graphs.
- `path`: Path-aware heuristic that uses major paths through the graph to guide the search. Best for graphs with clear path structure.

```bash
# Use path-aware heuristic for complex pangenome graphs
poasta align -H path sequences.fasta

# Use Dijkstra for guaranteed exploration of all options
poasta align -H dijkstra sequences.fasta
```

The heuristic only affects the search order, not the final alignment result - all heuristics find the optimal alignment.

## Aligning sequences with `lasagna`

`lasagna` aligns reads to an existing POA graph described in GFA format and
emits the resulting alignments as GAF records.

```bash
lasagna align graph.gfa reads.fq.gz > alignments.gaf
```

### Arguments

- `graph` – the input GFA graph (must be acyclic)
- `sequences` – FASTA/FASTQ reads to align (gzip supported)
- `-j`, `--num-threads` – number of worker threads
- `-o` – output filename (stdout if omitted)
- `-O` – output format, currently only `gaf`
- `-m` – alignment span (`global`, `semi-global`, `ends-free`)
- `-n` – mismatch penalty
- `-g` – gap open penalty
- `-e` – gap extend penalty

### How it differs from `poasta`

While `poasta` builds and updates partial order graphs from input sequences and
can convert them between several formats, `lasagna` does not modify the graph.
Instead it maps sequences onto a given GFA graph and outputs the alignments in
GAF without producing a POASTA graph file.

## Python visualization helpers

The `contrib/poasta_tools` directory provides small utilities for
inspecting alignment results. `poasta_plot.py` can plot the aligner
state over time while `poasta_graphviz_region.py` extracts and
visualizes subgraphs from a POA graph.

```bash
python contrib/poasta_tools/poasta_plot.py --help
python contrib/poasta_tools/poasta_graphviz_region.py --help
```

These tools depend on `numpy`, `pandas`, `matplotlib`, `networkx`, and
`pygraphviz`.

## Related repositories

* [poa-bench](https://github.com/broadinstitute/poa-bench) - Benchmark POASTA against other POA tools
* [spoa-rs](https://github.com/broadinstitute/spoa-rs) - Rust bindings to SPOA

## Contributing

Run the test suites to ensure everything works as expected.

```bash
cargo test
pytest tests/python_tools
```

Both suites should pass. Set the `SEED` environment variable to force
deterministic behavior. If you measure coverage, tools such as
`cargo tarpaulin` or `pytest --cov` are recommended.
