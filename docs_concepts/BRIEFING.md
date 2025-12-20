# GraphBin: Complete Study Guide and Developer Briefing

## Table of Contents

1. [Introduction](#introduction)
2. [Domain Knowledge: Metagenomics](#domain-knowledge-metagenomics)
3. [Domain Knowledge: Genome Assembly](#domain-knowledge-genome-assembly)
4. [Domain Knowledge: Binning](#domain-knowledge-binning)
5. [The GraphBin Approach](#the-graphbin-approach)
6. [Technical Architecture](#technical-architecture)
7. [Algorithm Details](#algorithm-details)
8. [Supported Assemblers](#supported-assemblers)
9. [Code Organization](#code-organization)
10. [Development Guide](#development-guide)
11. [Research Context](#research-context)

---

## Introduction

**GraphBin** is a metagenomic contig bin refinement tool that leverages assembly graph connectivity information to improve binning results. It represents a novel approach in computational metagenomics by being the first tool to explicitly use assembly graph topology for contig binning refinement.

### What Problem Does GraphBin Solve?

Existing metagenomic binning tools (like MaxBin2, MetaBAT2, CONCOCT, VAMB) use features such as:

- Contig composition (tetranucleotide frequencies, GC content)
- Contig coverage (sequencing depth)

However, they **ignore** the valuable connectivity information present in assembly graphs. GraphBin fills this gap by:

1. **Correcting mis-binned contigs** - fixing errors from initial binning tools
2. **Binning short contigs** - assigning bins to contigs that were too short for initial binners
3. **Utilizing graph topology** - leveraging the fact that connected contigs likely come from the same genome

---

## Domain Knowledge: Metagenomics

### What is Metagenomics?

**Metagenomics** is the study of genetic material recovered directly from environmental samples, rather than from cultured organisms. Unlike traditional genomics (studying one organism), metagenomics analyzes complex mixtures of organisms.

### Key Concepts

#### 1. Microbial Communities

- Environmental samples (soil, ocean water, human gut) contain hundreds to thousands of different species
- Traditional culturing methods can only grow ~1% of environmental microbes
- Metagenomics allows culture-independent study of entire communities

#### 2. Next-Generation Sequencing (NGS)

- Modern DNA sequencing technologies generate millions of short DNA fragments (reads)
- Read lengths: typically 100-300 base pairs (bp) for short-read technologies (Illumina)
- Long-read technologies (PacBio, Oxford Nanopore): 10,000-100,000+ bp but with higher error rates
- From a metagenomic sample, reads come from multiple different organisms

#### 3. The Metagenomic Workflow

```
Environmental Sample
    ↓
DNA Extraction
    ↓
Sequencing (NGS) → Reads (millions of short DNA fragments)
    ↓
Assembly → Contigs (longer assembled sequences)
    ↓
Binning → Bins (groups of contigs from same organism)
    ↓
Analysis (taxonomic classification, functional annotation)
```

---

## Domain Knowledge: Genome Assembly

### Why Assembly?

Individual sequencing reads are too short to analyze meaningfully. Assembly reconstructs longer sequences (contigs) by overlapping reads based on sequence similarity.

### Assembly Approaches

#### 1. De Bruijn Graph (DBG) Approach

**Used by:** SPAdes, metaSPAdes, MEGAHIT, Flye

**How it works:**

- Break reads into k-mers (substrings of length k)
- Example: Read "ATGGCA" with k=3 → k-mers: ATG, TGG, GGC, GCA
- Create nodes for (k-1)-mers: AT, TG, GG, GC, CA
- Create edges representing k-mers
- Paths through the graph represent assembled sequences

**Advantages:**

- Efficient for large datasets
- Good for handling repeats
- Memory efficient with modern implementations

**Graph Structure:**

- Nodes: k-mers or (k-1)-mers
- Edges: overlaps between k-mers
- Contigs: paths through the graph

#### 2. Overlap-Layout-Consensus (OLC) / String Graph Approach

**Used by:** SGA (String Graph Assembler), Canu, Miniasm

**How it works:**

1. **Overlap:** Find all pairwise overlaps between reads
2. **Layout:** Determine the order of reads based on overlaps
3. **Consensus:** Generate consensus sequences from overlapping reads

**String Graph refinement:**

- Modern OLC methods use string graphs (simplified overlap graphs)
- Remove transitive edges (redundant connections)
- More compact representation

**Advantages:**

- Better for long reads
- Preserves more information about read relationships
- Good for repetitive genomes

**Graph Structure:**

- Nodes: reads or contigs
- Edges: overlaps between sequences
- Contigs: paths representing merged reads

### Assembly Graph Outputs

Different assemblers produce different graph formats:

| Assembler         | Graph Format | Additional Files                            |
| ----------------- | ------------ | ------------------------------------------- |
| SPAdes/metaSPAdes | GFA (.gfa)   | contigs.paths - maps contigs to graph paths |
| MEGAHIT           | GFA (.gfa)   | -                                           |
| SGA               | ASQG (.asqg) | -                                           |
| Flye              | GFA (.gfa)   | assembly_info.txt - contig metadata         |
| Canu              | GFA (.gfa)   | -                                           |
| Miniasm           | GFA (.gfa)   | -                                           |

### What is a Contig?

A **contig** (contiguous sequence) is a continuous DNA sequence assembled from overlapping reads. Key properties:

- Length: varies from hundreds to millions of base pairs
- Coverage: average number of reads covering each position
- Represents a path or sub-path through the assembly graph

---

## Domain Knowledge: Binning

### What is Binning?

**Binning** is the process of grouping contigs that originated from the same organism (or genome). It's like sorting puzzle pieces by which puzzle they belong to.

### Why is Binning Challenging?

1. **Sequence similarity:** Related species have similar DNA
2. **Uneven sequencing depth:** Some organisms are more abundant than others
3. **Horizontal gene transfer:** Genes can move between species
4. **Strain variation:** Multiple strains of same species may be present
5. **Short contigs:** Lack enough signal for reliable binning
6. **Chimeric assemblies:** Assembly errors can combine sequences from different organisms

### Traditional Binning Features

#### 1. Composition-based Features

- **Tetranucleotide frequency:** Frequency of 4-nucleotide patterns (256 features)
- **GC content:** Ratio of G and C bases
- Rationale: Different organisms have different genomic "signatures"

#### 2. Coverage-based Features

- **Read depth:** How many reads cover each position
- Rationale: Organisms with higher abundance have higher coverage

### Existing Binning Tools

| Tool       | Method                                   | Strengths                         |
| ---------- | ---------------------------------------- | --------------------------------- |
| MaxBin 2.0 | Expectation-Maximization                 | Good for high-abundance organisms |
| MetaBAT 2  | Probability binning                      | Fast, good recall                 |
| CONCOCT    | Gaussian mixture models                  | Handles strain variation          |
| VAMB       | Variational autoencoders (deep learning) | Good for complex communities      |
| SolidBin   | Semi-supervised deep learning            | Improved precision                |

### Binning Output Format

Standard format (CSV):

```
contig_identifier,bin_identifier
NODE_1_length_458813_cov_136.660185,bin_1
NODE_2_length_409135_cov_127.776630,bin_1
NODE_3_length_346431_cov_35.887787,bin_2
```

Each contig is assigned to exactly one bin (or left unbinned).

---

## The GraphBin Approach

### Core Innovation

GraphBin adds a **third dimension** to binning: **connectivity information from assembly graphs**.

**Key Insight:** If two contigs are connected in the assembly graph, they likely come from the same genome.

### Why Graph Topology Matters

#### Example Scenario:

```
Contig A --edge-- Contig B --edge-- Contig C
  (bin 1)          (bin 2?)         (bin 1)
```

If Contig B is connected to both A and C, and A and C are in bin_1, then B should probably also be in bin_1, even if its composition/coverage features are ambiguous.

### The GraphBin Workflow

```
1. Input:
   - Assembly graph (GFA/ASQG)
   - Contigs (FASTA)
   - Initial binning result (CSV from MaxBin/MetaBAT/etc.)

2. Graph Construction:
   - Parse assembly graph
   - Map contigs to graph nodes
   - Build connectivity network

3. Label Refinement Phase 1:
   - Identify ambiguous vertices (conflicting neighbor labels)
   - Remove labels from contigs with mixed-label neighborhoods

4. Label Propagation:
   - Use semi-supervised learning algorithm
   - Propagate labels from labeled to unlabeled nodes
   - Based on graph connectivity and edge weights

5. Label Refinement Phase 2:
   - Re-check for ambiguous assignments
   - Remove labels that create conflicts

6. Output:
   - Refined binning result (CSV)
   - Binned FASTA files
```

### What GraphBin Does NOT Do

- **Does not perform de novo binning** - requires an initial binning result
- **Does not assemble contigs** - requires pre-assembled contigs and graph
- **Does not evaluate bin quality** - focuses on refinement, not validation

---

## Technical Architecture

### High-Level Component Overview

```
graphbin (CLI entry point)
    ↓
Assembler-specific modules
    ├── graphbin_SPAdes.py
    ├── graphbin_SGA.py
    ├── graphbin_MEGAHIT.py
    ├── graphbin_Flye.py
    ├── graphbin_Canu.py
    └── graphbin_Miniasm.py
    ↓
Parser modules (extract graph structure)
    ├── spades_parser.py
    ├── sga_parser.py
    ├── megahit_parser.py
    └── ... (others)
    ↓
graphbin_Func.py (core logic)
    ↓
Label Propagation Algorithm
    └── labelprop.py
```

### Key Data Structures

#### 1. BidirectionalMap

```python
# Maps node IDs ↔ contig identifiers
contig_names[node_id] = "NODE_1_length_458813_cov_136.660185"
contig_names.inverse["NODE_1_length_458813_cov_136.660185"] = node_id
```

**Purpose:** Quick bidirectional lookup between internal node IDs and contig names

#### 2. Assembly Graph (igraph)

```python
# Directed graph representing contig connectivity
assembly_graph = Graph(directed=True)
assembly_graph.add_vertices(node_count)
assembly_graph.add_edges(edge_list)
```

**Properties:**

- Nodes: contigs
- Edges: connections from assembly process
- Can query neighbors: `assembly_graph.neighbors(node_id, mode="all")`

#### 3. Bins

```python
# List of lists: bins[bin_idx] = [contig_id1, contig_id2, ...]
bins = [[] for x in range(n_bins)]
bins[0] = [1, 5, 7, 12]  # Contigs in bin 0
bins[1] = [2, 3, 9]       # Contigs in bin 1
```

#### 4. Final Bins (Dictionary)

```python
# Maps contig_id → bin_name
final_bins = {
    1: "bin_1",
    5: "bin_1",
    7: "bin_1",
    2: "bin_2",
}
```

---

## Algorithm Details

### Phase 1: Identify and Remove Ambiguous Labels

**Goal:** Remove labels from contigs that have neighbors in multiple different bins.

**Algorithm:**

```python
for each bin:
    for each contig in bin:
        my_bin = current_bin
        neighbors = get_neighbors_from_graph(contig)

        # Check if all labeled neighbors are in the same bin
        neighbors_have_same_label = True
        for neighbor in neighbors:
            if neighbor is in a different bin:
                neighbors_have_same_label = False
                break

        if not neighbors_have_same_label:
            mark contig as ambiguous
            remove label
```

**Why?** These are likely mis-binned contigs or contigs at the boundary between genomes.

### Phase 2: Identify Multi-Level Ambiguous Vertices

**Goal:** Find contigs whose closest labeled neighbors disagree.

**Algorithm:**

```python
for each contig:
    closest_labeled_neighbors = BFS_to_find_closest_labeled_nodes(contig)

    # Check if these neighbors agree
    bin_votes = count_bins_in(closest_labeled_neighbors)

    if len(bin_votes) > 1:  # Disagreement
        mark_as_ambiguous(contig)
        remove_label(contig)
```

**Uses Breadth-First Search (BFS)** to find the nearest labeled contigs in the graph.

### Phase 3: Label Propagation

GraphBin uses a **semi-supervised label propagation algorithm** based on [Zhu and Ghahramani (2002)](http://mlg.eng.cam.ac.uk/zoubin/papers/CMU-CALD-02-107.pdf).

#### Key Concepts

**Semi-supervised learning:** Uses both labeled data (initial bins) and unlabeled data (unbinned/ambiguous contigs).

**Label propagation principle:** Labels "flow" from labeled nodes to unlabeled nodes through graph edges, with labels spreading more easily to nearby/connected nodes.

#### Algorithm Outline

1. **Initialization:**

   ```python
   For labeled nodes: set probability = 1.0 for true label, 0.0 for others
   For unlabeled nodes: set all probabilities = 0.0
   ```

2. **Propagation (iterative):**

   ```python
   for iteration in range(max_iterations):
       for each node:
           # Weighted average of neighbor labels
           new_label_distribution = weighted_sum(
               neighbor_label_distributions,
               edge_weights
           )

           # Normalize by node degree
           new_label_distribution /= node_degree

       # Check convergence
       if changes < diff_threshold:
           break
   ```

3. **Assignment:**
   ```python
   for each unlabeled node:
       assign bin with highest probability
   ```

#### Parameters

- **max_iteration** (default: 100): Maximum propagation iterations
- **diff_threshold** (default: 0.1): Convergence threshold (stop if changes < threshold)

#### Mathematical Foundation

For node $i$ with label distribution $\mathbf{Y}_i$:

$$\mathbf{Y}_i^{(t+1)} = \frac{1}{d_i} \sum_{j \in N(i)} w_{ij} \mathbf{Y}_j^{(t)}$$

Where:

- $d_i$ = degree of node $i$
- $N(i)$ = neighbors of node $i$
- $w_{ij}$ = edge weight between $i$ and $j$
- $t$ = iteration number

### Phase 4: Final Ambiguity Removal

After label propagation, repeat ambiguity detection (Phase 1) to remove any newly created conflicts.

**Why needed?** Label propagation might assign labels that conflict with immediate neighbors.

---

## Supported Assemblers

### Short-Read Assemblers (De Bruijn Graph)

#### 1. SPAdes / metaSPAdes

- **Graph format:** GFA (Graphical Fragment Assembly)
- **Special requirement:** contigs.paths file
- **Why?** Maps contigs to paths through the de Bruijn graph
- **Use case:** Most popular for metagenomic assembly

#### 2. MEGAHIT

- **Graph format:** GFA (converted from FASTG)
- **Special requirement:** None
- **Use case:** Fast, memory-efficient for large metagenomes

### Short-Read Assemblers (String Graph)

#### 3. SGA (String Graph Assembler)

- **Graph format:** ASQG (A-Bruijn/String Graph format)
- **Special requirement:** None
- **Use case:** Alternative approach, less common now

### Long-Read Assemblers

#### 4. Flye

- **Graph format:** GFA
- **Special requirement:** assembly_info.txt
- **Approach:** Repeat graph (variant of de Bruijn)
- **Use case:** Nanopore/PacBio metagenomic assembly

#### 5. Canu

- **Graph format:** GFA
- **Approach:** OLC
- **Use case:** Long-read assembly, older tool

#### 6. Miniasm

- **Graph format:** GFA
- **Approach:** OLC
- **Use case:** Fast long-read assembly

---

## Code Organization

### Directory Structure

```
GraphBin/
├── src/graphbin/          # Main source code
│   ├── cli.py             # Command-line interface (entry point)
│   ├── graphbin_SPAdes.py # SPAdes-specific workflow
│   ├── graphbin_SGA.py    # SGA-specific workflow
│   ├── graphbin_MEGAHIT.py# MEGAHIT-specific workflow
│   ├── graphbin_Flye.py   # Flye-specific workflow
│   ├── graphbin_Canu.py   # Canu-specific workflow
│   ├── graphbin_Miniasm.py# Miniasm-specific workflow
│   ├── graphbin_Func.py   # Core binning logic (shared)
│   ├── parsers/           # Graph parsing modules
│   │   ├── spades_parser.py
│   │   ├── sga_parser.py
│   │   ├── megahit_parser.py
│   │   ├── flye_parser.py
│   │   ├── canu_parser.py
│   │   └── miniasm_parser.py
│   ├── labelpropagation/ # Label propagation algorithm
│   │   └── labelprop.py
│   ├── bidirectionalmap/ # Utility data structure
│   │   └── bidirectionalmap.py
│   ├── support/          # Support scripts
│   │   ├── prep_result.py  # Format initial binning results
│   │   ├── gfa2fasta.py    # Convert GFA to FASTA
│   │   └── visualise_*.py  # Visualization scripts
│   └── utils/            # Alternative implementations (older)
├── tests/                # Test suite
│   ├── test_graphbin.py  # Integration tests
│   ├── test_arguments.py # CLI argument tests
│   ├── test_bidirectional_map.py
│   └── data/             # Test datasets
├── docs/                 # Documentation
├── pyproject.toml        # Project configuration
├── requirements.txt      # Dependencies
└── README.md
```

### Key Modules Explained

#### cli.py

- Entry point for command-line interface
- Uses Click library for argument parsing
- Validates inputs (file existence, parameter ranges)
- Dispatches to appropriate assembler-specific module
- Sets up logging

#### graphbin\_[Assembler].py

Each assembler module follows the same pattern:

1. Parse command-line arguments
2. Read initial binning result
3. Parse assembly graph (assembler-specific)
4. Map contigs to graph nodes
5. Call core GraphBin logic (graphbin_Func.py)
6. Write output files

#### graphbin_Func.py

Core refinement logic (assembler-agnostic):

- `getClosestLabelledVertices()`: BFS to find nearest labeled nodes
- `graphbin_main()`: Main refinement algorithm
  - Phase 1 ambiguity removal
  - Phase 2 multi-level ambiguity removal
  - Label propagation
  - Phase 4 final ambiguity removal

#### Parsers

Each parser handles assembler-specific graph formats:

- **Read graph file:** Extract nodes and edges
- **Map contigs:** Create bidirectional mapping between contig names and node IDs
- **Build igraph:** Construct graph object for analysis
- **Parse metadata:** Handle assembler-specific information (paths, coverage, etc.)

#### labelprop.py

Semi-supervised label propagation implementation:

- `LabelProp` class: Main algorithm
- `Edge` class: Represents graph edges with weights
- `load_data_from_mem()`: Initialize from graph data
- `run()`: Execute label propagation iterations
- Based on Zhu & Ghahramani algorithm

### Data Flow

```
1. CLI Input
   ↓
2. cli.py validates and dispatches
   ↓
3. graphbin_[Assembler].py
   ├── Parser: assembly graph → igraph object
   ├── Parser: initial bins → bins data structure
   └── Calls graphbin_Func.graphbin_main()
       ↓
4. graphbin_Func.py
   ├── Ambiguity detection
   ├── Label propagation (via labelprop.py)
   └── Final refinement
       ↓
5. graphbin_[Assembler].py writes output
   ├── CSV file with refined bins
   └── FASTA files for each bin
```

---

## Development Guide

### Setting Up Development Environment

```bash
# Clone repository
git clone https://github.com/metagentools/GraphBin.git
cd GraphBin

# Install flit (build tool)
pip install flit

# Install GraphBin in development mode
flit install -s --python `which python`

# Install development dependencies
pip install -e ".[dev]"
```

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=graphbin --cov-report=html

# Run specific test
pytest tests/test_graphbin.py::test_graphbin_on_spades_dataset
```

### Code Style

GraphBin follows **PEP 8** with these tools:

```bash
# Format code with black
black src/graphbin/

# Sort imports with isort
isort src/graphbin/

# Check style
flake8 src/graphbin/
```

### Adding Support for a New Assembler

1. **Create parser in `parsers/`:**

   ```python
   def parse_graph(graph_file, ...):
       # Parse graph format
       # Return: assembly_graph, contig_map, node_count

   def get_initial_binning_result(...):
       # Map bins to graph nodes
       # Return: bins list

   def write_output(...):
       # Write refined binning result
   ```

2. **Create main module:**

   ```python
   # graphbin_NewAssembler.py
   def run(args):
       # Parse inputs
       assembly_graph, contigs_map, node_count = parse_graph(...)
       bins = get_initial_binning_result(...)

       # Run core logic
       final_bins, remove_labels, non_isolated = graphbin_main(...)

       # Write output
       write_output(...)
   ```

3. **Update CLI:**

   ```python
   # In cli.py
   @click.option(
       "--assembler",
       type=click.Choice([..., "newassembler"], ...),
       ...
   )

   # In main()
   if assembler.lower() == "newassembler":
       graphbin_NewAssembler.main(args)
   ```

4. **Add tests:**
   ```python
   # In tests/test_graphbin.py
   def test_graphbin_on_newassembler_dataset(tmp_dir):
       ...
   ```

### Key Dependencies

| Library       | Purpose                             | Version  |
| ------------- | ----------------------------------- | -------- |
| python-igraph | Graph data structure and algorithms | Latest   |
| cogent3       | FASTA file parsing                  | Latest   |
| cairocffi     | Graph visualization (optional)      | Latest   |
| click         | Command-line interface              | Latest   |
| pytest        | Testing                             | Dev only |
| black         | Code formatting                     | Dev only |

### Debugging Tips

1. **Enable debug logging:**

   ```python
   # Logs are automatically written to output_folder/graphbin.log
   # Set level to DEBUG for verbose output
   ```

2. **Check intermediate results:**

   ```python
   # Add print statements in graphbin_Func.py
   print(f"Bins after phase 1: {bins}")
   print(f"Removed labels: {remove_labels}")
   ```

3. **Visualize assembly graph:**

   ```bash
   # Use support scripts
   python src/graphbin/support/visualise_result_SPAdes.py \
       --graph assembly_graph.gfa \
       --paths contigs.paths \
       --binned binning_result.csv \
       --output visualization.png
   ```

4. **Test with small datasets:**
   - Use test data in `tests/data/`
   - Create minimal test cases for quick iteration

---

## Research Context

### Original Publication

**Title:** GraphBin: refined binning of metagenomic contigs using assembly graphs

**Authors:** Vijini Mallawaarachchi, Anuradha Wickramarachchi, Yu Lin

**Journal:** Bioinformatics, Volume 36, Issue 11, June 2020

**DOI:** [10.1093/bioinformatics/btaa180](https://doi.org/10.1093/bioinformatics/btaa180)

### Scientific Contribution

**First tool to use assembly graph topology for metagenomic binning refinement.**

Key findings from the paper:

1. Assembly graphs contain valuable connectivity information ignored by existing binners
2. Graph-based refinement can correct mis-binned contigs
3. Graph-based propagation can bin short contigs (often discarded)
4. Works with both de Bruijn and OLC assembly graphs
5. Improves results across multiple initial binning tools

### Related Tools and Concepts

#### MetaWRAP

- Pipeline that combines multiple binning tools
- GraphBin can be used as a refinement step

#### DAS Tool

- Aggregates results from multiple binners
- Different approach: consensus voting vs. graph-based

#### Graph-based binning successors

GraphBin inspired follow-up work:

- **GraphBin2:** Improved version with refined clustering
- Research into other graph-based features

### Current Research Directions

1. **Deep learning on graphs:**

   - Graph neural networks (GNNs) for binning
   - Learn better representations of graph structure

2. **Long-read specific methods:**

   - Exploit unique properties of long-read graphs
   - Handle higher error rates

3. **Strain-level resolution:**

   - Distinguish closely related strains
   - Handle graph complexity from strain variation

4. **Hybrid approaches:**
   - Combine composition, coverage, AND graph topology
   - End-to-end binning (not just refinement)

### Evaluation Metrics

Understanding how GraphBin is evaluated:

#### Precision

Fraction of binned contigs that are correctly binned

```
Precision = True Positives / (True Positives + False Positives)
```

#### Recall

Fraction of all contigs that are correctly binned

```
Recall = True Positives / (True Positives + False Negatives)
```

#### F1 Score

Harmonic mean of precision and recall

```
F1 = 2 × (Precision × Recall) / (Precision + Recall)
```

#### ARI (Adjusted Rand Index)

Measures clustering quality, adjusted for chance

- Range: -1 to 1
- 1 = perfect clustering
- 0 = random clustering

### Funding and Community

- **Funded by:** Chan Zuckerberg Initiative (Essential Open Source Software for Science Grant)
- **Organization:** metagentools (GitHub organization)
- **Community:** Part of the broader computational metagenomics ecosystem

---

## Practical Usage Scenarios

### Scenario 1: Basic Refinement

You have metaSPAdes assembly and MaxBin2 binning:

```bash
# Run MaxBin2
run_MaxBin.pl -contig contigs.fasta -out maxbin_result ...

# Format MaxBin2 output
python prep_result.py --binned maxbin_bins/ --output ./

# Run GraphBin
graphbin --assembler spades \
    --graph assembly_graph_with_scaffolds.gfa \
    --contigs contigs.fasta \
    --paths contigs.paths \
    --binned initial_contig_bins.csv \
    --output refined_bins/
```

### Scenario 2: Long-Read Assembly

Using Flye for long-read assembly:

```bash
# Assemble with Flye
flye --nano-raw reads.fastq --out-dir flye_output --meta

# Bin with VAMB
vamb --fasta flye_output/assembly.fasta ...

# Refine with GraphBin
graphbin --assembler flye \
    --graph flye_output/assembly_graph.gfa \
    --contigs flye_output/assembly.fasta \
    --paths flye_output/assembly_info.txt \
    --binned vamb_bins.csv \
    --output refined_bins/
```

### Scenario 3: Custom Parameters

For complex communities, adjust propagation parameters:

```bash
graphbin --assembler spades \
    --graph assembly_graph.gfa \
    --contigs contigs.fasta \
    --paths contigs.paths \
    --binned initial_bins.csv \
    --output refined/ \
    --max_iteration 200 \      # More iterations
    --diff_threshold 0.05      # Stricter convergence
```

---

## Common Issues and Solutions

### Issue 1: Contigs Not in Graph

**Problem:** Some contigs from binning result not found in assembly graph.

**Cause:** Different assembler outputs or filtering.

**Solution:**

- Ensure using same assembly for both graphing and binning
- Check that graph file includes all contigs
- Verify contig naming consistency

### Issue 2: No Refinement Observed

**Problem:** Output identical to input.

**Possible causes:**

1. Very high-quality initial binning (nothing to refine)
2. Disconnected graph (isolated contigs)
3. Ambiguous neighborhoods everywhere

**Solution:**

- Check graph connectivity
- Visualize graph to understand structure
- Try different diff_threshold values

### Issue 3: Memory Issues

**Problem:** Out of memory error.

**Cause:** Very large metagenomic assembly graphs.

**Solutions:**

- Filter short contigs before binning
- Use more memory-efficient assembler (MEGAHIT)
- Increase system memory
- Split into smaller samples

### Issue 4: MEGAHIT FASTG Conversion

**Problem:** MEGAHIT outputs FASTG, not GFA.

**Solution:**

```bash
# Use fastg2gfa converter
# See: https://github.com/lh3/gfa1
```

---

## Future Development Ideas

### Technical Improvements

1. **Parallelization:** Multi-threading for label propagation
2. **GPU acceleration:** For large graphs
3. **Streaming algorithms:** Handle graphs too large for memory
4. **Better convergence detection:** Adaptive diff_threshold

### Algorithmic Extensions

1. **Edge weighting:** Use coverage/composition similarity
2. **Multiple propagation rounds:** With different parameters
3. **Consensus mechanisms:** Aggregate multiple refinement runs
4. **Strain-aware refinement:** Handle within-species variation

### Tool Integration

1. **Direct tool interfaces:** Call MaxBin/MetaBAT directly
2. **Pipeline integration:** MetaWRAP, Snakemake workflows
3. **Visualization dashboard:** Interactive graph exploration
4. **Quality metrics:** Built-in bin quality assessment

---

## Glossary

| Term                  | Definition                                                  |
| --------------------- | ----------------------------------------------------------- |
| **Assembly**          | Process of reconstructing longer sequences from short reads |
| **Bin**               | Group of contigs from the same organism                     |
| **Binning**           | Grouping contigs by organism of origin                      |
| **Contig**            | Contiguous sequence assembled from overlapping reads        |
| **Coverage**          | Average number of reads covering each position              |
| **De Bruijn Graph**   | Graph where nodes are k-mers and edges represent overlaps   |
| **GFA**               | Graphical Fragment Assembly format                          |
| **K-mer**             | Substring of length k                                       |
| **Label Propagation** | Semi-supervised algorithm that spreads labels through graph |
| **Metagenome**        | Genetic material from environmental sample                  |
| **NGS**               | Next-Generation Sequencing                                  |
| **OLC**               | Overlap-Layout-Consensus assembly approach                  |
| **Read**              | Short DNA sequence fragment from sequencer                  |
| **Refinement**        | Improving existing binning result                           |
| **String Graph**      | Simplified overlap graph for assembly                       |

---

## References and Further Reading

### Papers

1. **GraphBin paper:** Mallawaarachchi et al., Bioinformatics 2020
2. **Label propagation:** Zhu & Ghahramani, CMU-CALD-02-107, 2002
3. **metaSPAdes:** Nurk et al., Genome Research 2017
4. **MetaBAT 2:** Kang et al., PeerJ 2019

### Tools

- [metaSPAdes](http://cab.spbu.ru/software/spades/)
- [MEGAHIT](https://github.com/voutcn/megahit)
- [MaxBin2](https://sourceforge.net/projects/maxbin2/)
- [MetaBAT2](https://bitbucket.org/berkeleylab/metabat)

### Courses and Tutorials

- [Metagenomics Workshop](https://github.com/ngs-docs/angus)
- [Assembly Tutorial](https://github.com/rrwick/Bandage/wiki)

---

## Quick Reference Card

### Installation

```bash
conda install -c bioconda graphbin
# OR
pip install graphbin
```

### Basic Usage

```bash
graphbin --assembler spades \
    --graph assembly_graph.gfa \
    --contigs contigs.fasta \
    --paths contigs.paths \
    --binned initial_bins.csv \
    --output results/
```

### Key Parameters

- `--max_iteration`: Max label propagation iterations (default: 100)
- `--diff_threshold`: Convergence threshold (default: 0.1)
- `--delimiter`: CSV delimiter (default: comma)

### Output Files

- `graphbin_output.csv`: Refined binning result
- `graphbin_bin_*.fasta`: FASTA files for each bin
- `graphbin.log`: Detailed log file

### Getting Help

- Documentation: https://graphbin.readthedocs.io/
- Issues: https://github.com/metagentools/GraphBin/issues
- Paper: https://doi.org/10.1093/bioinformatics/btaa180

---

**Last Updated:** December 2025  
**GraphBin Version:** 1.7.4  
**Maintained by:** Vijini Mallawaarachchi, Anuradha Wickramarachchi, Yu Lin
