# Domain Knowledge: Genome Assembly - Deep Dive

## Table of Contents

1. [Introduction to Genome Assembly](#introduction-to-genome-assembly)
2. [The Assembly Problem](#the-assembly-problem)
3. [De Bruijn Graph Assembly](#de-bruijn-graph-assembly)
4. [Overlap-Layout-Consensus Assembly](#overlap-layout-consensus-assembly)
5. [String Graph Assembly](#string-graph-assembly)
6. [Assembly Graph Structures](#assembly-graph-structures)
7. [Metagenomic Assembly Challenges](#metagenomic-assembly-challenges)
8. [Assembly Evaluation](#assembly-evaluation)
9. [Assembler Implementations](#assembler-implementations)
10. [Best Practices](#best-practices)

---

## Introduction to Genome Assembly

### The Fundamental Problem

Genome sequencing produces millions of short DNA fragments (reads), but we need the complete genome sequence. **Assembly** is the computational process of reconstructing longer sequences from these overlapping fragments.

**Analogy:** Imagine trying to reconstruct a book from millions of small, overlapping snippets of text, where:

- Some pages appear multiple times (coverage)
- Some identical sections appear in different chapters (repeats)
- Some snippets contain typos (sequencing errors)
- You don't know the order of snippets (assembly challenge)

### Historical Context

**Pre-Sequencing Era:**

- Genetic maps based on recombination
- Physical maps using restriction enzymes
- Limited to small regions

**Sanger Era (1977-2005):**

- First genomes: φX174 bacteriophage (1977), 5,386 bp
- Human Genome Project (1990-2003): $3 billion, 13 years
- Clone-by-clone approach: Break genome into manageable pieces
- Hierarchical sequencing strategy

**NGS Era (2005-Present):**

- Whole genome shotgun sequencing
- Computational assembly becomes central
- Can sequence and assemble genomes in days/weeks
- Metagenomic assembly: Multiple genomes simultaneously

### Why Assembly is Necessary

**Direct read analysis limitations:**

1. **Too short:** Reads (100-300 bp) much shorter than genes (1-10 kb)
2. **No context:** Can't determine gene order, structure
3. **Fragmented annotation:** Partial genes not informative
4. **Comparative genomics:** Need complete sequences for comparison
5. **Variant calling:** Need reference for mapping

**Assembly benefits:**

- Longer sequences (contigs/scaffolds)
- Gene prediction and annotation
- Structural variant detection
- Comparative analysis
- Complete genome reconstruction (in ideal cases)

---

## The Assembly Problem

### Mathematical Formulation

**Shortest Superstring Problem:**
Given a set of strings (reads), find the shortest string that contains all of them as substrings.

```
Input: {ATGC, TGCA, GCAT}
Output: ATGCAT (contains all three)
```

**Challenge:** NP-hard problem (computationally intractable for large inputs)

**Graph Formulation:**
Assembly can be represented as finding paths through graphs:

- **Nodes:** Sequences (reads or k-mers)
- **Edges:** Overlaps between sequences
- **Assembly:** Path that visits all nodes/edges

### Key Challenges

#### 1. Repeats

**Types of Repeats:**

**Exact Repeats:**

```
Genome:  ...ATCG[REPEAT]AAAA[REPEAT]GGGG...
Reads:        TCGREPEAT
              REPEATAAAA
              REPEATGGGG
```

Problem: Can't determine which repeat instance read came from.

**Near-Identical Repeats:**

```
Repeat A: ATCGATCGATCGATCG
Repeat B: ATCGATCGATTGATCG (one SNP difference)
```

Even harder: Sequencing errors can mask differences.

**Common Repeat Classes:**

- **Transposable elements:** Mobile DNA, can be thousands of copies
- **Tandem repeats:** Direct adjacency (e.g., microsatellites)
- **Inverted repeats:** Palindromic sequences
- **Segmental duplications:** Large (>1 kb) duplicated regions
- **rRNA operons:** Multiple copies of ribosomal RNA genes

**Impact on Assembly:**

- Breaks contigs (gaps in assembly)
- Collapses repeats (assembly shorter than true genome)
- Mis-assemblies (wrong connections)
- Ambiguous paths in graph

#### 2. Sequencing Errors

**Error Types:**

**Substitution Errors:**

```
True:  ATCGATCG
Read:  ATCGATGG (C→G substitution)
```

**Insertion Errors:**

```
True:  ATCGATCG
Read:  ATCGAATCG (extra A)
```

**Deletion Errors:**

```
True:  ATCGATCG
Read:  ATCGTCG (T deleted)
```

**Error Rates by Technology:**

- Illumina: ~0.1-1% (mostly substitutions)
- PacBio CLR: ~10-15% (mostly insertions/deletions)
- PacBio HiFi: ~0.1-1% (after consensus)
- Nanopore: ~5-15% (all types, improving)

**Impact on Assembly:**

- False branches in graph (error creates new k-mer)
- Increased memory usage
- Need for error correction
- "Tips" and "bubbles" in graph

#### 3. Coverage Variation

**Uneven Coverage:**

```
Coverage: ||||||||||||  |||  ||||||||||||||||||  |  |||||
Genome:   ─────────────────────────────────────────────────
          Easy         Hard  Easy              Hard Easy
```

**Causes:**

- **GC bias:** High/low GC regions sequence poorly
- **PCR bias:** Amplification not uniform
- **DNA secondary structure:** Hairpins resist sequencing
- **Biological variation:** Copy number variants, ploidy

**Metagenomic-Specific:**

- **Abundance:** Dominant organisms 1000× more covered than rare
- **Differential lysis:** Some organisms harder to extract from
- **Host DNA:** Can swamp microbial DNA

**Impact on Assembly:**

- Low coverage: Gaps, fragmentation
- High coverage: Computational burden
- Variable coverage: Ambiguous graph regions

#### 4. Polyploidy and Heterozygosity

**Diploid/Polyploid Organisms:**

```
Chromosome 1a: ATCGATCGATCG
Chromosome 1b: ATCGATCGCTCG (SNP at position 9)
```

**Challenges:**

- Should haplotypes be assembled separately or together?
- Heterozygous SNPs create bubbles in graph
- Structural variants between haplotypes

**Metagenomic Context:**

- Strains: Multiple closely related strains in community
- Each strain like a "haplotype"
- GraphBin helps distinguish at binning stage

#### 5. Computational Complexity

**Resource Requirements:**

**Memory:**

- Store graph structure
- Human genome: 10-200 GB RAM depending on algorithm
- Metagenome: 100-500+ GB RAM

**Time:**

- Graph construction: Hours to days
- Graph simplification: Hours
- Path finding: Can be exponential

**Scalability:**

- Metagenomes: Gigabases to terabases
- Thousands of genomes simultaneously
- Need efficient data structures and algorithms

---

## De Bruijn Graph Assembly

### Fundamental Concepts

#### K-mers

**Definition:** A k-mer is a substring of length k.

**Example:**

```
Sequence: ATCGATCG
3-mers:   ATC, TCG, CGA, GAT, ATC, TCG
4-mers:   ATCG, TCGA, CGAT, GATC, ATCG
5-mers:   ATCGA, TCGAT, CGATC, GATCG
```

**Key Property:** Two k-mers overlap by k-1 bases if they're consecutive in sequence.

```
K-mer 1: ATCG
K-mer 2:  TCGA
Overlap:  TCG (k-1 = 3 bases)
```

#### De Bruijn Graph Construction

**Classic Definition:**

- **Nodes:** (k-1)-mers
- **Edges:** k-mers
- **Edge from node A to B:** k-mer starting with A, ending with B

**Example (k=3):**

```
Sequence: ATCGATCG

K-mers (k=3): ATC, TCG, CGA, GAT, ATC, TCG

(k-1)-mers: AT, TC, CG, GA, AT, TC

Graph:
   AT → TC → CG
         ↓    ↓
        GA → AT → TC
         ↑_________|
```

**Alternative Definition (Used by SPAdes, MEGAHIT):**

- **Nodes:** k-mers
- **Edges:** (k-1) overlap
- Mathematically equivalent, different implementation

#### Why De Bruijn Graphs?

**Advantages over Overlap Graphs:**

1. **Computational Efficiency:**

   - Overlap graph: O(n²) comparisons (all reads vs. all reads)
   - De Bruijn graph: O(n) (linear in number of k-mers)

2. **Memory Efficiency:**

   - Store k-mers, not full reads
   - Hash table or bloom filter data structures
   - Can process very large datasets

3. **Natural Repeat Handling:**
   - Repeats automatically merge in graph
   - Create branching structure
   - Explicit representation of ambiguity

**Disadvantages:**

- Information loss (read connectivity)
- K-mer size parameter critical
- Sensitive to errors (error creates new k-mer)

### Algorithm Steps

#### 1. K-mer Counting

```python
def count_kmers(reads, k):
    kmer_counts = {}
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    return kmer_counts
```

**Purpose:**

- Count occurrence of each k-mer
- Filter rare k-mers (likely errors)
- Estimate coverage

**Typical Threshold:**

- Keep k-mers appearing ≥2-3 times
- Removes most error k-mers

#### 2. Graph Construction

```python
def build_debruijn_graph(kmers, k):
    graph = {}  # adjacency list
    for kmer in kmers:
        prefix = kmer[:k-1]  # First k-1 bases
        suffix = kmer[1:]    # Last k-1 bases
        if prefix not in graph:
            graph[prefix] = []
        graph[prefix].append(suffix)
    return graph
```

**Result:** Directed graph where paths represent potential sequences

#### 3. Error Correction

**Tips:** Dead-end branches, usually errors

```
Main path: ────────────────
            ↓
Error tip: ───── (short dead end)
```

**Action:** Remove tips shorter than threshold

**Bubbles:** Parallel paths, can be errors or variants

```
        ╭──── Path A ────╮
Node ──┤                 ├─── Node
        ╰──── Path B ────╯
```

**Action:**

- If paths very similar: Merge (error correction)
- If paths different: Keep (true variant)

**Low Coverage Regions:**

- Remove nodes/edges with coverage < threshold
- Balances error removal vs. losing real data

#### 4. Graph Simplification

**Unambiguous Path Compression:**

```
Before: A → B → C → D (4 nodes, linear)
After:  A → BCD (2 nodes, one edge labeled "BCD")
```

**Purpose:**

- Reduce graph complexity
- Create longer contigs from unambiguous regions
- Easier visualization and downstream analysis

**Resulting Structure:**

- Nodes: Compressed contigs
- Edges: Ambiguous connections
- Contigs: Simple paths in simplified graph

#### 5. Path Traversal / Scaffolding

**Challenge:** Multiple valid paths through graph

**Approaches:**

**Greedy:**

- Follow highest coverage path
- Simple, fast
- May miss optimal solution

**Paired-End Information:**

- Two reads separated by known distance
- Resolves ambiguities
- Links contigs into scaffolds

**Coverage-Based:**

- Follow paths maintaining consistent coverage
- Works well for single genomes
- Challenging for metagenomes

**Graph Segmentation:**

- Identify connected components
- Assemble each separately
- Useful for metagenomes

### K-mer Size Selection

**Trade-offs:**

#### Small k (e.g., k=21-31)

**Advantages:**

- More overlaps found
- Better connectivity
- Longer contigs in low-coverage regions

**Disadvantages:**

- Repeats collapse
- More false connections
- Larger graph (more k-mers)

#### Large k (e.g., k=101-127)

**Advantages:**

- Resolves more repeats
- Fewer false connections
- More specificity

**Disadvantages:**

- Fewer overlaps
- Fragmented assembly
- Requires higher coverage

**Practical Approach:**

- **Multiple k-mer sizes:** Use several k values, merge results
- **Adaptive k:** Start large, decrease for unassembled regions
- **SPAdes approach:** Iteratively increase k, build graph at each stage

### SPAdes De Bruijn Graph Specifics

**Multi-Sized De Bruijn Graph:**

```
k=21 → k=33 → k=55 → k=77 → k=99 → k=127
 ↓      ↓      ↓      ↓      ↓      ↓
G21 → G33  → G55  → G77  → G99  → G127
       ↑      ↑      ↑      ↑      ↑
       Map contigs between graphs
```

**Process:**

1. Start with small k (good connectivity)
2. Map contigs to next larger k graph
3. Extend contigs using larger k
4. Resolve ambiguities at each stage
5. Final assembly combines information from all k

**Advantages:**

- Gets benefits of both small and large k
- Better repeat resolution
- Longer, more accurate contigs

### MEGAHIT De Bruijn Graph Specifics

**Succinct De Bruijn Graph:**

- Space-efficient representation
- Stores only essential information
- Can handle large metagenomes with limited RAM

**Iterative Approach:**

- Similar to SPAdes: Multiple k values
- More aggressive memory optimization
- Trade-off: Slightly lower quality for much lower memory

**Typical Parameters:**

```
--k-min 21 --k-max 141 --k-step 12
Creates graphs at k = 21, 33, 45, 57, ..., 141
```

---

## Overlap-Layout-Consensus Assembly

### Fundamental Concepts

**Core Idea:** Find overlaps between reads, determine their layout, and generate consensus.

**Three Phases:**

#### 1. Overlap Phase

**Goal:** Find all significant overlaps between reads

**All-vs-All Comparison:**

```
Read 1: ATCGATCGATCG
Read 2:       TCGATCGAAAA
        ──────────── (overlap)
```

**Challenges:**

- **Computational:** O(n²) comparisons naively
- **Sensitivity:** Need to detect overlaps despite errors
- **Specificity:** Avoid false overlaps from repeats

**Optimization Strategies:**

**Seed-and-Extend:**

1. Find short exact matches (seeds)
2. Extend to full overlap
3. Common approach: BWT, FM-index, minimizers

**Minhash:**

- Sketch reads with representative k-mers
- Compare sketches (much faster)
- Candidates for full overlap check

**Overlap Criteria:**

```
Minimum overlap length: 30-50 bp typically
Maximum error rate: 5-15% depending on technology
```

#### 2. Layout Phase

**Goal:** Determine the order and orientation of reads

**Overlap Graph:**

- **Nodes:** Reads
- **Edges:** Overlaps between reads
- **Layout:** Path through graph visiting all nodes

**Example:**

```
Reads:    R1      R2      R3      R4
          ──────
              ──────
                  ──────
                      ──────
Layout:   R1 → R2 → R3 → R4
```

**Challenges:**

**Repeat-Induced Tangles:**

```
        ╭─── R2 ───╮
R1 ────┤           ├──── R4
        ╰─── R3 ───╯
(R2 and R3 from different repeat copies)
```

**Solution Strategies:**

- **Best Overlap:** For each read, keep only best overlap(s)
- **Transitive Reduction:** Remove redundant edges
- **Coverage Analysis:** Use coverage to detect repeats

#### 3. Consensus Phase

**Goal:** Generate consensus sequence from aligned reads

**Multiple Sequence Alignment:**

```
Read 1: ATCGATCG
Read 2: ATCGATTG (error at position 7)
Read 3: ATCGATCG
Read 4: ATCGATCG
        ||||||||
Consensus: ATCGATCG (majority vote)
```

**Considerations:**

- **Quality scores:** Weight by base quality
- **Coverage:** More reads → more confidence
- **Errors:** Identify and correct
- **Variants:** True biological variation vs. errors

### String Graph Assembly

**Motivation:** Overlap graphs are redundant

**Transitive Edges:**

```
A ────→ B ────→ C
 ╰──────────────╯ (redundant if A→B→C already connected)
```

**String Graph:**

- Remove transitive edges
- Simplify overlap graph
- More compact representation
- Easier path finding

**Construction:**

1. Build overlap graph
2. Identify transitive edges
3. Remove transitive edges
4. **Contigs:** Simple paths in string graph

**SGA Algorithm:**

- Uses FM-index for efficient overlap finding
- String graph for layout
- Focuses on high-quality, error-corrected reads

### OLC vs. De Bruijn

| Aspect                 | OLC                              | De Bruijn                           |
| ---------------------- | -------------------------------- | ----------------------------------- |
| **Computational Cost** | O(n²) (optimized to ~O(n log n)) | O(n)                                |
| **Memory**             | High (store reads + overlaps)    | Lower (k-mers)                      |
| **Read Length**        | Better for long reads            | Originally for short reads          |
| **Repeat Handling**    | Explicit in graph                | Automatic merging                   |
| **Information Loss**   | Minimal (keeps read pairs)       | Some (k-mer decomposition)          |
| **Error Sensitivity**  | Lower (longer overlaps)          | Higher (single error affects k-mer) |
| **Typical Use**        | Long-read assemblers             | Short-read assemblers               |

---

## String Graph Assembly

### Motivation and Theory

**Problem with Overlap Graphs:**

- Many redundant edges (transitive)
- Graph size can be huge
- Path-finding computationally expensive

**String Graph Solution:**

- Represents same information
- Dramatically fewer edges
- Provably correct (Myers 2005)

**Formal Definition:**

**Transitive Edge:**
Edge v → w is transitive if there exists path v → ... → w through other vertices with equal or greater overlap.

```
Example:
v ─────100bp─────→ w
 ╰──50bp→ x ──60bp→╯

v → w is transitive (removed)
Keep: v → x → w
```

### SGA (String Graph Assembler)

**Specific Implementation:**

#### 1. Read Preprocessing

**Error Correction:**

- **K-mer based:** Correct to most common k-mer variant
- **Overlap based:** Align reads, correct low-quality discrepancies
- Critical because errors create false branches

**Read Filtering:**

- Remove low-quality reads
- Remove duplicates
- Trim low-quality ends

#### 2. Overlap Computation

**FM-Index:**

- Data structure for efficient substring search
- Based on Burrows-Wheeler Transform (BWT)
- Enables fast overlap finding

**Process:**

```
1. Build FM-index of all reads
2. For each read:
   a. Query FM-index with suffixes
   b. Find all significant overlaps
   c. Store overlap information
```

**Irreducible Overlaps:**

- Overlap that's not containment
- Read not completely contained in another
- These form string graph edges

#### 3. String Graph Construction

**Algorithm:**

```python
def build_string_graph(overlaps):
    graph = {}

    # Add all overlap edges
    for overlap in overlaps:
        add_edge(graph, overlap)

    # Remove transitive edges
    for vertex in graph:
        for edge1 in outgoing_edges(vertex):
            for edge2 in outgoing_edges(edge1.dest):
                if is_transitive(edge1, edge2):
                    mark_for_removal(edge1 or edge2)

    remove_marked_edges(graph)
    return graph
```

**Transitive Reduction:**

- Identify three-vertex paths (v → x → w)
- Check if direct edge v → w exists
- Remove transitive edge
- Preserve graph connectivity

#### 4. Graph Simplification

**Tip Removal:**

```
Main: ─────────────────
       ↓
Tip:  ───── (short branch)
```

**Bubble Popping:**

```
     ╭──── Variant A ────╮
Node┤                    ├─Node
     ╰──── Variant B ────╯
```

**Criteria:**

- Similar length
- High sequence similarity
- Likely sequencing errors or SNPs

#### 5. Contig Extraction

**Simple Path Traversal:**

- Start at unambiguous node
- Follow until branching or end
- Output as contig

**Paired-End Scaffolding:**

- Use read pairs to link contigs
- Estimate gap sizes
- Create scaffolds

### ASQG Format

**ASQG:** Assembly Sequence Graph (used by SGA)

**Structure:**

```
HT      VT:i:1  # Header: Version
VT      seq1    ATCGATCGATCG...  # Vertex: Read
VT      seq2    GCTAGCTAGCTA...
ED      seq1 seq2 10 20 30 40 0 0  # Edge: Overlap
```

**Edge Format:**

```
ED  read1  read2  start1  end1  start2  end2  orient  score
```

**Information:**

- Which reads overlap
- Position of overlap on each read
- Orientation (forward/reverse complement)
- Overlap quality

---

## Assembly Graph Structures

### Graph Formats

#### GFA (Graphical Fragment Assembly)

**Most Common Format:** Used by SPAdes, MEGAHIT, Flye, Canu, Miniasm

**GFA 1.0 Structure:**

```
H       VN:Z:1.0
S       1       ATCGATCG        LN:i:8  RC:i:10
S       2       GCTAGCTA        LN:i:8  RC:i:15
L       1       +       2       +       4M
P       path1   1+,2+   4M
```

**Line Types:**

**H (Header):**

```
H       VN:Z:1.0        # Version
```

**S (Segment):**

```
S       <name>  <sequence>      [tags]
S       contig1 ATCGATCG        LN:i:8  RC:i:10

Tags:
- LN:i:length
- RC:i:read_count (coverage)
- KC:i:kmer_count
- FC:i:fragment_count
```

**L (Link):**

```
L       <from>  <from_orient>   <to>    <to_orient>     <overlap>
L       1       +       2       +       4M

Orientations: + (forward), - (reverse complement)
Overlap: CIGAR-like (M=match)
```

**P (Path):**

```
P       <name>  <segments>      <overlaps>
P       contig1 1+,2+,3-        4M,5M

Describes a path through segments
Used by SPAdes for contig paths
```

**GFA 2.0:**

- More features (gaps, fragments)
- Less widely adopted
- Backward compatible

#### SPAdes-Specific Graph Features

**Assembly Graph:**

```
assembly_graph.gfa          # Main graph
assembly_graph.fastg        # Alternative format (deprecated)
```

**Contigs.paths:**

```
NODE_1_length_1000_cov_10.5
1+,2+,3-,4+
NODE_2_length_500_cov_5.2
5-,6+
;
```

**Format:**

- Contig name
- Path through graph (segment IDs + orientations)
- Semicolon separates entries

**Why Needed:**

- GFA has segments (k-mers or unitigs)
- Contigs are paths through segments
- Maps contigs to graph structure
- Essential for GraphBin!

#### ASQG Format Details

**Used by:** SGA

**Structure:**

```
HT      VT:i:1
VT      read1   SEQUENCE...     QV:Z:!!!...
VT      read2   SEQUENCE...     QV:Z:!!!...
ED      read1 read2 10 100 10 100 0 99  ER:i:0
```

**Advantages:**

- Includes quality values
- Explicit overlap coordinates
- Error rate information

**Disadvantages:**

- Less compact than GFA
- Less tool support
- Specific to SGA

### Graph Topology Patterns

#### Linear Regions (Unambiguous)

```
A ──→ B ──→ C ──→ D
```

**Characteristics:**

- One path
- Easy to assemble
- High confidence
- Becomes contigs

**Biological Meaning:**

- Unique sequences
- No repeats
- Or repeats longer than reads can span

#### Branches (Ambiguity)

```
     ╭──→ B
A ──┤
     ╰──→ C
```

**Causes:**

1. **Repeats:** A appears multiple times in genome
2. **Strain variation:** Different strains/haplotypes
3. **Errors:** Sequencing/assembly errors

**Resolution:**

- Paired-end info
- Coverage analysis
- Long reads spanning branch
- GraphBin-type approaches!

#### Bubbles (Variants or Errors)

```
     ╭──→ B ──╮
A ──┤         ├──→ D
     ╰──→ C ──╯
```

**Interpretation:**

1. **SNP bubble:** B and C differ by single nucleotide
2. **Indel bubble:** B and C differ by insertion/deletion
3. **Error bubble:** One path is error

**Tools often:**

- Pop small bubbles (errors)
- Keep larger bubbles (true variants)

#### Tips (Dead Ends)

```
Main: A ──→ B ──→ C ──→ D
       ↓
Tip:   E ──→ F
```

**Causes:**

- Sequencing errors
- Coverage drop-off
- End of linear chromosome
- Unconnected contigs

**Handling:**

- Remove short tips
- Keep longer tips (may be real)

#### Cycles (Circular DNA)

```
A ──→ B ──→ C ──→ D
↑                  ↓
╰──────────────────╯
```

**Biological:**

- Bacterial chromosomes (often circular)
- Plasmids
- Viral genomes

**Assembly Challenge:**

- Where to start/end?
- Need to detect circularity

---

## Metagenomic Assembly Challenges

### Multiple Genomes Simultaneously

**Complexity:**

- Single genome: ~1-10 Mb, one set of repeats
- Metagenome: 1-100 Gb, thousands of repeat sets
- Graph becomes highly interconnected

**Shared Sequences:**

```
Species A: ...ATCG[CONSERVED GENE]GCTA...
Species B: ...TTAA[CONSERVED GENE]CCGG...
Species C: ...GGCC[CONSERVED GENE]AATT...

Graph:    ╭─ A
          │
Gene ────┼─ B
          │
          ╰─ C  (impossible to resolve without additional info)
```

**Example:** 16S rRNA gene

- Present in all bacteria
- Highly conserved regions
- Creates massive tangles in graph

### Strain Variation

**Closely Related Strains:**

```
Strain A: ATCGATCGATCG  (SNP at position 7)
Strain B: ATCGATTGATCG
Similarity: >99.9%
```

**Graph Result:**

- Small bubbles everywhere
- Should we keep strains separate or merge?
- Depends on application!

**Challenges:**

1. **Coverage confusion:** Hard to distinguish strains vs. coverage variation
2. **Chimeras:** Assembly may create artificial hybrid genomes
3. **Fragmentation:** Overly cautious → breaks at every SNP
4. **Mis-assemblies:** Too aggressive → wrong combinations

### Uneven Abundance

**Typical Distribution:**

- 10 dominant species: 80% of reads
- 100 moderate species: 15% of reads
- 1000+ rare species: 5% of reads

**Coverage Range:** 0.01× to 1000× in same sample!

**Impact:**

**High Abundance Species:**

- Very high coverage (100-1000×)
- Good assembly quality
- May overshadow others

**Low Abundance Species:**

- Low coverage (0.1-10×)
- Fragmented assembly
- May be missed entirely

**Computational:**

- Need to handle huge coverage variation
- Memory for high-coverage nodes
- Sensitivity for low-coverage nodes

### Chimerism

**Assembly Chimeras:**
Artificial sequences combining pieces from different organisms.

**Causes:**

1. **Shared Sequences:**

```
Organism A: ...AAAA[SHARED]TTTT...
Organism B: ...GGGG[SHARED]CCCC...

Chimera:    ...AAAA[SHARED]CCCC...  (wrong!)
```

2. **Strain Mixing:**

```
Strain 1: ATCG[Region A from S1][Region B from S1]GCTA
Strain 2: ATCG[Region A from S2][Region B from S2]GCTA

Chimera:  ATCG[Region A from S1][Region B from S2]GCTA
```

**Detection:**

- Coverage discontinuities
- Composition shifts (GC content, tetranucleotide)
- Taxonomic inconsistencies
- Read pair evidence

**Prevention:**

- Conservative assembly parameters
- Use coverage information
- Binning before/during assembly

### Metagenomic Assembler Strategies

**metaSPAdes:**

- Modified de Bruijn graph traversal
- Coverage-aware graph simplification
- Keeps alternative paths (don't merge strains)
- Iterative assembly

**MEGAHIT:**

- Memory-efficient de Bruijn graph
- Handles very large metagenomes
- Succinct data structures
- Multiple k-mer sizes

**metaFlye:**

- Repeat graph approach
- Built for long reads
- Better repeat resolution
- Strain-aware assembly

---

## Assembly Evaluation

### Contiguity Metrics

#### N50

**Definition:** Median contig length (weighted by length)

**Calculation:**

1. Order contigs by length (longest first)
2. Calculate cumulative length
3. N50 = length of contig at 50% cumulative length

**Example:**

```
Contigs: 100, 90, 80, 70, 60, 50, 40, 30, 20, 10 bp
Total: 550 bp
Cumulative: 100, 190, 270, 340, 400, 450, 490, 520, 540, 550
50% of 550 = 275 bp
N50 = 80 bp (contig where cumulative exceeds 275)
```

**Interpretation:**

- Higher is better
- Half of assembly is in contigs ≥ N50
- Commonly reported metric
- Can be misleading (longer contigs may be mis-assemblies)

#### L50

**Definition:** Number of contigs that make up 50% of assembly

**Example (same data):**

```
L50 = 3 (need contigs 100 + 90 + 80 to reach 50%)
```

**Interpretation:**

- Lower is better
- Measures how many contigs needed
- Complements N50

#### NG50

**Definition:** Like N50 but relative to genome size (if known)

**Better for:**

- Comparing assemblies of same genome
- Accounts for different total assembly sizes

#### Other N-statistics

- **N90:** 90% of assembly in contigs ≥ this length
- **N10:** 10% of assembly in contigs ≥ this length
- **auN:** Area under Nx curve (single metric)

### Completeness Metrics

#### Total Assembly Size

```
Assembly size vs. expected genome size(s)
```

**Issues:**

- Collapsed repeats → smaller
- Mis-assemblies → larger
- Contamination → larger

#### BUSCO (Benchmarking Universal Single-Copy Orthologs)

**Concept:**

- Genes expected to be present in single copy
- Taxon-specific gene sets (bacteria, fungi, plants, etc.)
- Check how many present in assembly

**Scores:**

- Complete: Gene found, complete
- Fragmented: Gene found, partial
- Missing: Gene not found
- Duplicated: Gene found multiple times (potential issue)

**Example:**

```
BUSCO Results:
C:95%[S:90%,D:5%],F:3%,M:2%,n:124

C = Complete (95%)
S = Single-copy (90%)
D = Duplicated (5%)
F = Fragmented (3%)
M = Missing (2%)
n = Total BUSCOs checked (124)
```

**Interpretation:**

- High C: Good completeness
- High D: Potential contamination or redundancy
- High F/M: Fragmented or incomplete assembly

#### CheckM (for Bacterial/Archaeal Genomes)

**Approach:**

- Similar to BUSCO but specifically for prokaryotes
- Uses marker gene sets
- Also estimates contamination

**Metrics:**

- **Completeness:** % of expected genes found
- **Contamination:** % of genes appearing multiple times (should be single-copy)
- **Strain heterogeneity:** Variation within genome

**Quality Tiers:**

```
High quality:     >90% complete, <5% contamination
Medium quality:   ≥50% complete, <10% contamination
Low quality:      <50% complete or >10% contamination
```

### Accuracy Metrics

#### Mis-assembly Detection

**Reference-Based (if reference available):**

- Align assembly to reference
- Identify:
  - Relocations
  - Translocations
  - Inversions
  - Structural variants

**Reference-Free:**

- **Read alignment:** Align reads back to assembly
  - Proper pairs: Both mates align correctly
  - Discordant pairs: Wrong distance/orientation
  - Clipped reads: Soft-clipping indicates issues
- **Coverage:** Sudden coverage changes → chimera
- **Composition:** GC%, tetranucleotide shifts

**Tools:**

- QUAST: With/without reference
- REAPR: Read-based evaluation
- ALE: Assembly Likelihood Estimation

#### Indel/SNP Error Rate

**Method:**

- Align reads to assembly
- Call variants
- High SNP rate → assembly errors

**Expected:**

- Should match sequencing error rate
- Excess indicates assembly issues

### Metagenomic-Specific Evaluation

#### Per-Species Metrics

**Challenges:**

- Don't know true number of species
- No reference genomes for most

**Approaches:**

1. **Binning + CheckM:** Evaluate each bin separately
2. **Simulations:** Create synthetic metagenomes with known composition
3. **Mock communities:** Sequencing of known mixtures

#### Assembly Graph Metrics

**Relevant for GraphBin:**

**Graph Connectivity:**

- Number of connected components
- Proportion of isolated contigs
- Graph diameter

**Graph Complexity:**

- Number of ambiguous branching points
- Bubble count
- Tip count

**Implications:**

- More connected → more information for GraphBin
- Highly fragmented → limited graph benefit
- Many bubbles → strain variation or errors

---

## Assembler Implementations

### SPAdes / metaSPAdes

**Algorithm:** Multi-k de Bruijn graph

**Key Features:**

- Iterative k-mer approach
- Paired-end scaffolding
- Error correction
- Repeat resolution using paired reads

**metaSPAdes Modifications:**

- Coverage-aware graph traversal
- Keep alternative paths (strains)
- Modified repeat resolution
- Chimera detection

**Performance:**

- Memory: Moderate to high
- Time: Moderate
- Quality: High for complex communities

**Best For:**

- Complex metagenomes
- When high quality is priority
- Sufficient computational resources available

### MEGAHIT

**Algorithm:** Succinct de Bruijn graph

**Key Features:**

- Memory-efficient data structures
- Multiple k-mer sizes
- Fast execution
- Scales to very large datasets

**Data Structures:**

- Succinct representation of graph
- Compressed k-mer storage
- Minimal memory footprint

**Performance:**

- Memory: Low
- Time: Fast
- Quality: Good (slightly lower than metaSPAdes)

**Best For:**

- Very large metagenomes (>100 Gb)
- Limited computational resources
- Speed is priority

### SGA

**Algorithm:** String graph (OLC-based)

**Key Features:**

- Overlap-based
- FM-index for efficiency
- String graph simplification
- Good for longer reads

**Error Correction:**

- K-mer based
- Overlap based
- Built into pipeline

**Performance:**

- Memory: Moderate
- Time: Moderate to slow
- Quality: High

**Best For:**

- High-quality reads
- When overlap-based approach preferred
- Comparative studies

### Flye / metaFlye

**Algorithm:** Repeat graph

**Key Features:**

- Designed for long reads
- Handles high error rates
- Repeat-aware graph construction
- Polishing with long/short reads

**Repeat Graph:**

- Explicitly represents repeat structure
- Different from de Bruijn graph
- Better repeat resolution

**Performance:**

- Memory: Moderate to high
- Time: Moderate
- Quality: Excellent with long reads

**Best For:**

- PacBio, Nanopore data
- Complex metagenomes with good long-read coverage
- When complete genomes desired

### Canu

**Algorithm:** OLC for long reads

**Key Features:**

- Correction → Trimming → Assembly pipeline
- Sophisticated error correction
- Handles high error rates
- Repeat-aware

**Stages:**

1. **Correction:** Fix sequencing errors using overlaps
2. **Trimming:** Remove low-quality regions
3. **Assembly:** OLC with corrected, trimmed reads

**Performance:**

- Memory: High
- Time: Slow (thorough correction)
- Quality: Excellent for single genomes

**Best For:**

- High-quality genome assemblies
- When computation time not limiting
- Complex genomes with lots of repeats

### Miniasm

**Algorithm:** Fast OLC

**Key Features:**

- Minimal error correction
- Fast overlap finding
- String graph assembly
- Requires polishing

**Philosophy:**

- Assembly first, polish later
- Speed over initial accuracy
- Relies on downstream polishing (Racon, Pilon)

**Performance:**

- Memory: Moderate
- Time: Very fast
- Quality: Draft (needs polishing)

**Best For:**

- Quick assembly for exploration
- When will polish anyway
- Time-sensitive applications

---

## Best Practices

### Parameter Selection

#### K-mer Size

**Rules of Thumb:**

- **Short reads:** k=31-127 (odd numbers)
- **Long reads:** k=1000+ for overlap finding
- **Multiple k:** Use 3-5 different k values

**Considerations:**

- Coverage: Higher coverage → can use larger k
- Error rate: Higher errors → use smaller k
- Genome complexity: More repeats → try larger k

#### Coverage Cutoff

**Purpose:** Remove error k-mers

**Typical Values:**

- **High coverage (>50×):** cutoff = 3-5
- **Medium coverage (20-50×):** cutoff = 2-3
- **Low coverage (<20×):** cutoff = 1-2

**Auto-detection:**

- Many assemblers detect automatically
- Based on k-mer frequency histogram

#### Error Correction

**When to Use:**

- Always for Illumina (low error rate makes it safe)
- Carefully for long reads (may over-correct real variants)

**Iterative:**

- Correct reads
- Assemble
- Polish with raw reads
- Better than single-pass

### Quality Control

**Pre-Assembly:**

1. **Read QC:** Check quality (FastQC)
2. **Trimming:** Remove adapters, low quality
3. **Filtering:** Remove short reads, contamination
4. **Error correction:** If appropriate

**Post-Assembly:**

1. **Contiguity:** N50, L50
2. **Completeness:** BUSCO, CheckM
3. **Accuracy:** Map reads back, check for errors
4. **Contamination:** Taxonomic classification of contigs

### Computational Resources

**Memory Estimation:**

```
Rule of thumb:
Illumina: 1 Gb RAM per 1 Gb of reads (minimum)
Better: 2-5 Gb RAM per 1 Gb of reads

Metagenome example:
100 Gb reads → 200-500 Gb RAM recommended
```

**Time:**

- Assembly: Hours to days
- Error correction: Often the slowest step
- Plan for 2-10× the shortest tool's time

**Parallelization:**

- Use all available cores
- Many assemblers support multi-threading
- Some steps inherently serial

### Troubleshooting

**Fragmented Assembly:**

- Increase coverage
- Use multiple k-mer sizes
- Try different assembler
- Long reads for scaffolding

**Mis-assemblies:**

- More conservative parameters
- Better error correction
- Validate with read pairs
- Use reference if available

**Out of Memory:**

- Use MEGAHIT (more memory efficient)
- Filter rare k-mers more aggressively
- Downsample reads
- Use cluster/cloud computing

**Too Slow:**

- Use faster assembler (MEGAHIT, Miniasm)
- Reduce k-mer sizes tried
- Downsample (if high coverage)
- Parallelize on cluster

---

## Summary

Genome assembly is the critical first step in analyzing sequencing data, converting short reads into longer, more useful sequences. Two main approaches dominate:

**De Bruijn Graph Assembly:**

- Efficient, scalable
- Used by most modern short-read assemblers
- Naturally handles repeats (as branches)
- SPAdes, MEGAHIT use this approach

**Overlap-Layout-Consensus / String Graph Assembly:**

- More direct representation of read relationships
- Better for long reads
- String graphs simplify overlap graphs
- SGA, Canu, Miniasm use variants of this

**For Metagenomics:**

- Challenges multiply: multiple genomes, strain variation, uneven coverage
- Specialized assemblers needed (metaSPAdes, MEGAHIT, metaFlye)
- Assembly graphs become highly complex
- **This is where GraphBin becomes valuable!**

**GraphBin Context:**

- Works with assembly graphs from any assembler
- Uses graph topology to improve binning
- Requires understanding of graph structure
- Benefits from high-quality assemblies

Understanding assembly helps developers:

- Parse different graph formats (GFA, ASQG)
- Interpret graph structure (repeats, strain variation)
- Understand limitations and artifacts
- Make informed algorithmic choices

---

## Additional Resources

### Papers

**Foundational:**

- Pevzner et al. (2001): "An Eulerian path approach to DNA fragment assembly"
- Myers (2005): "The fragment assembly string graph"
- Zerbino & Birney (2008): "Velvet: Algorithms for de novo short read assembly"

**Modern Assemblers:**

- Bankevich et al. (2012): "SPAdes: A new genome assembly algorithm"
- Li et al. (2015): "MEGAHIT: An ultra-fast single-node solution"
- Nurk et al. (2017): "metaSPAdes: A new versatile metagenomic assembler"
- Kolmogorov et al. (2019): "Assembly of long, error-prone reads using repeat graphs"

### Software

**Assemblers:**

- SPAdes: http://cab.spbu.ru/software/spades/
- MEGAHIT: https://github.com/voutcn/megahit
- Flye: https://github.com/fenderglass/Flye
- SGA: https://github.com/jts/sga

**Evaluation:**

- QUAST: http://bioinf.spbau.ru/quast
- BUSCO: https://busco.ezlab.org/
- CheckM: https://ecogenomics.github.io/CheckM/

**Visualization:**

- Bandage: https://rrwick.github.io/Bandage/
- AGB: https://github.com/almiheenko/AGB

### Tutorials

- SPAdes tutorial: http://cab.spbu.ru/software/spades/
- Flye tutorial: https://github.com/fenderglass/Flye/wiki
- Assembly workshop: https://github.com/rrwick/SLING2023-Workshop

---

**Last Updated:** December 2025  
**Part of:** GraphBin Developer Documentation
