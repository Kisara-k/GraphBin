# Bioinformatics Foundations: A Comprehensive Guide

## Table of Contents

1. [Introduction](#introduction)
2. [Sequence Alignment](#sequence-alignment)
3. [Sequence Motifs and Patterns](#sequence-motifs-and-patterns)
4. [Phylogenetic Analysis](#phylogenetic-analysis)
5. [Genomic Databases and Resources](#genomic-databases-and-resources)
6. [File Formats in Bioinformatics](#file-formats-in-bioinformatics)
7. [Genome Assembly](#genome-assembly)
8. [Read Mapping and Variant Analysis](#read-mapping-and-variant-analysis)
9. [Computational Tools and Programming](#computational-tools-and-programming)
10. [Summary and Learning Resources](#summary-and-learning-resources)

---

## Introduction

Bioinformatics is the interdisciplinary field that develops and applies computational methods to analyze and interpret biological data, particularly molecular sequence data. As genomic technologies generate unprecedented volumes of data, bioinformatics has become essential for extracting biological insights.

### What is Bioinformatics?

**Definition:** The application of computational and statistical techniques to understand biological data.

**Core Areas:**

- **Sequence Analysis:** Comparing, aligning, and analyzing DNA, RNA, and protein sequences
- **Genome Analysis:** Assembly, annotation, and comparative genomics
- **Structural Bioinformatics:** Protein structure prediction and analysis
- **Functional Genomics:** Understanding gene function and regulation
- **Systems Biology:** Modeling biological systems and networks
- **Phylogenetics:** Evolutionary relationships and tree construction

### Why Bioinformatics Matters

**The Data Explosion:**

- Human genome: 3 billion base pairs
- A single Illumina NovaSeq run: up to 6 terabases
- Thousands of genomes sequenced daily
- Without computational tools, this data would be meaningless

**Key Applications:**

- Identifying disease-causing mutations
- Drug discovery and development
- Understanding evolution and biodiversity
- Agricultural improvements
- Personalized medicine
- Infectious disease tracking (e.g., COVID-19 surveillance)
- Metagenomics and microbiome studies

### The Bioinformatics Workflow

**Typical Pipeline:**

1. **Data Generation:** Sequencing, mass spectrometry, etc.
2. **Quality Control:** Assess and filter data quality
3. **Preprocessing:** Trim, normalize, error correction
4. **Analysis:** Alignment, assembly, annotation
5. **Interpretation:** Biological meaning, statistical significance
6. **Visualization:** Present results clearly
7. **Validation:** Experimental confirmation

### Interdisciplinary Nature

Bioinformatics sits at the intersection of:

- **Biology:** Understanding biological questions and context
- **Computer Science:** Algorithms, data structures, software engineering
- **Statistics:** Experimental design, hypothesis testing, modeling
- **Mathematics:** Optimization, graph theory, probability
- **Domain Knowledge:** Chemistry, physics, medicine

### Historical Context

**Key Milestones:**

- **1965:** First protein sequence database (Margaret Dayhoff)
- **1970:** Needleman-Wunsch algorithm (global alignment)
- **1981:** Smith-Waterman algorithm (local alignment)
- **1982:** GenBank established
- **1990:** BLAST algorithm published
- **1995:** First bacterial genome sequenced (_H. influenzae_)
- **2001:** Human genome draft published
- **2005:** Next-generation sequencing introduced
- **2012:** ENCODE project reveals genome functionality
- **2020s:** Single-cell and spatial omics era

### Skills for Bioinformatics

**Essential:**

- Programming (Python, R, Bash)
- Statistics and probability
- Molecular biology fundamentals
- Command-line proficiency
- Algorithm understanding

**Important:**

- Database management (SQL)
- High-performance computing
- Version control (Git)
- Scientific writing and visualization
- Domain-specific knowledge

---

## Sequence Alignment

Sequence alignment is the fundamental operation in bioinformatics, identifying regions of similarity between biological sequences that may indicate functional, structural, or evolutionary relationships.

### Why Align Sequences?

**Purposes:**

- Identify homologous genes or proteins
- Infer evolutionary relationships
- Predict protein structure and function
- Identify conserved functional domains
- Detect mutations and variants
- Assemble genomes from sequencing reads

### Pairwise Alignment

Aligning two sequences to identify similarities.

#### Scoring Schemes

**Match/Mismatch:**

- Simple scheme: +1 for match, -1 for mismatch
- More sophisticated: Substitution matrices

**Substitution Matrices:**

**For DNA:**

- Simple scoring often sufficient
- Transition/transversion weighting

**For Proteins:**

**PAM (Point Accepted Mutation):**

- Based on evolutionary distance
- PAM1: 1% of amino acids changed
- PAM250: More distant relationships
- Higher PAM numbers = more divergent sequences

**BLOSUM (Blocks Substitution Matrix):**

- Based on conserved regions in aligned sequences
- BLOSUM62: Most commonly used (62% identity)
- Lower BLOSUM numbers = more divergent sequences
- Higher BLOSUM numbers = more similar sequences

**Example BLOSUM62 excerpt:**

```
     A   R   N   D   C   Q   E   G   H
A    4  -1  -2  -2   0  -1  -1   0  -2
R   -1   5   0  -2  -3   1   0  -2   0
N   -2   0   6   1  -3   0   0   0   1
```

**Gap Penalties:**

Gaps represent insertions or deletions (indels) in evolution.

- **Linear Gap Penalty:** Each gap costs the same
  - Score = -g × (gap length)
- **Affine Gap Penalty:** Opening a gap costs more than extending it
  - Score = -d - (e × (gap length - 1))
  - d = gap opening penalty
  - e = gap extension penalty
  - Models biology better (one indel event creates multiple gaps)

---

#### Global Alignment: Needleman-Wunsch Algorithm

**Purpose:** Align entire sequences end-to-end

**Algorithm:**

1. **Initialize:** Create matrix (m+1) × (n+1) for sequences of length m and n
2. **Fill first row/column:** Cumulative gap penalties
3. **Fill matrix:** For each cell (i,j):
   ```
   Score(i,j) = max(
       Score(i-1,j-1) + s(xi, yj),  # Diagonal (match/mismatch)
       Score(i-1,j) + gap_penalty,   # Up (gap in sequence y)
       Score(i,j-1) + gap_penalty    # Left (gap in sequence x)
   )
   ```
4. **Traceback:** Start from bottom-right, follow path to top-left

**Example:**

```
Sequence 1: ACGT
Sequence 2: AGT

     -   A   C   G   T
-    0  -1  -2  -3  -4
A   -1   1   0  -1  -2
G   -2   0   0   1   0
T   -3  -1  -1   0   2

Alignment:
A C G T
A - G T
```

**Complexity:** O(m × n) time and space

**When to Use:**

- Sequences of similar length
- Expected to be similar along entire length
- Homologous sequences

---

#### Local Alignment: Smith-Waterman Algorithm

**Purpose:** Find best matching subsequences

**Algorithm:**

Similar to Needleman-Wunsch, but:

1. **Negative scores set to 0** (can start alignment anywhere)
2. **Start traceback from maximum score** (not necessarily corner)
3. **Stop traceback at 0** (alignment can end anywhere)

**Example:**

```
Sequence 1: ACGTACGT
Sequence 2: CGTA

Best local alignment:
CGTA
CGTA
```

**When to Use:**

- Sequences of very different lengths
- Only part of sequences expected to be similar
- Domain identification in proteins
- Database searching

**Complexity:** O(m × n) time and space

---

#### Heuristic Alignment: BLAST

Smith-Waterman is too slow for database searches. BLAST (Basic Local Alignment Search Tool) uses heuristics for speed.

**Algorithm Overview:**

1. **Seeding:**

   - Break query into words (k-mers, typically k=3 for proteins, k=11 for DNA)
   - Find exact matches in database
   - Example: "ACGTACGT" → "ACG", "CGT", "GTA", etc.

2. **Extension:**

   - Extend seeds in both directions
   - Continue while score increases
   - Stop when score drops below threshold

3. **Evaluation:**
   - Calculate statistical significance (E-value)
   - Report significant matches

**BLAST Variants:**

| Program       | Query                   | Database                | Use Case                           |
| ------------- | ----------------------- | ----------------------- | ---------------------------------- |
| **blastn**    | Nucleotide              | Nucleotide              | Find similar DNA/RNA sequences     |
| **blastp**    | Protein                 | Protein                 | Find similar proteins              |
| **blastx**    | Nucleotide (translated) | Protein                 | Find proteins from DNA query       |
| **tblastn**   | Protein                 | Nucleotide (translated) | Find DNA encoding similar proteins |
| **tblastx**   | Nucleotide (translated) | Nucleotide (translated) | Compare sequences at protein level |
| **megablast** | Nucleotide              | Nucleotide              | Very similar sequences (faster)    |

**E-value (Expect Value):**

Number of alignments with this score expected by chance given database size.

- E-value = K × m × n × e^(-λS)
  - K, λ: statistical parameters
  - m: query length
  - n: database size
  - S: alignment score

**Interpretation:**

- E < 10^-50: Identical or nearly identical sequences
- E < 10^-10: Highly significant homology
- E < 10^-5: Significant homology (commonly used threshold)
- E < 10^-3: Suggestive of homology
- E > 0.01: Likely due to chance

**Bit Score:**

Normalized score independent of database size and search parameters.

- Higher bit score = better alignment
- Can compare across different searches

**Running BLAST:**

```bash
# Command-line BLAST
blastp -query protein.fasta \
       -db nr \
       -out results.txt \
       -evalue 1e-5 \
       -outfmt 6 \
       -num_threads 8

# Output format 6 (tabular):
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

**BLAST Parameters:**

- **Word size:** Smaller = more sensitive, slower
- **E-value threshold:** Lower = more stringent
- **Max target sequences:** Number of results to return
- **Scoring matrix:** BLOSUM62 (default for proteins)

**Alternatives to BLAST:**

- **DIAMOND:** 100-20,000× faster than BLAST, similar sensitivity
- **MMseqs2:** Ultra-fast homology search
- **HMMER:** Profile HMM-based search (more sensitive)

---

### Multiple Sequence Alignment (MSA)

Aligning three or more sequences simultaneously.

#### Why MSA?

- Identify conserved regions across multiple species
- Build phylogenetic trees
- Predict secondary structure
- Profile searches (PSSM, HMM)
- Functional site identification

#### Algorithms

**Progressive Alignment:**

Most common approach (used by ClustalW, MUSCLE, MAFFT)

**Steps:**

1. **Calculate pairwise distances** between all sequences
2. **Build guide tree** (clustering based on similarity)
3. **Progressive alignment** following tree:
   - Align most similar sequences first
   - Add more distant sequences progressively
   - Once a gap, always a gap

**Example:**

```
Seq1: ACGT
Seq2: ACCT
Seq3: AGGT

1. Pairwise distances → Guide tree: ((Seq1, Seq2), Seq3)
2. Align Seq1 and Seq2:
   ACGT
   ACCT
3. Align result with Seq3:
   ACGT
   ACCT
   AGGT
```

**Limitations:**

- Greedy algorithm (not guaranteed optimal)
- Early mistakes propagate
- Depends on guide tree quality

---

**Iterative Refinement:**

Improve initial alignment through iteration.

1. Create initial alignment
2. Remove one or more sequences
3. Realign to remaining alignment
4. Repeat until convergence

Used by MUSCLE, MAFFT (iterative mode)

---

**Consistency-Based Methods:**

Use pairwise alignment information to guide MSA.

- **T-Coffee:** Combines multiple pairwise alignments
- More accurate but slower

---

#### MSA Tools

**ClustalW / Clustal Omega:**

- Most widely used historically
- Progressive alignment
- Clustal Omega: Faster, more accurate
- Produces phylogenetic trees

**MUSCLE:**

- Fast and accurate
- Iterative refinement
- Good for large alignments

**MAFFT:**

- Very fast
- Multiple algorithms (FFT-based, iterative)
- Recommended for large datasets
- Supports very large alignments (10,000+ sequences)

**T-Coffee:**

- Very accurate
- Consistency-based
- Slower than others
- Good for divergent sequences

**PRANK:**

- Phylogeny-aware alignment
- Distinguishes insertions from deletions
- Better for evolutionary analysis

**Comparison:**

| Tool              | Speed     | Accuracy       | Best For              |
| ----------------- | --------- | -------------- | --------------------- |
| **Clustal Omega** | Medium    | Good           | General purpose       |
| **MUSCLE**        | Fast      | Good           | Medium-large datasets |
| **MAFFT**         | Very Fast | Good-Excellent | Large datasets        |
| **T-Coffee**      | Slow      | Excellent      | Small, divergent      |
| **PRANK**         | Slow      | Good           | Evolutionary studies  |

---

#### Viewing and Editing MSA

**Visualization Tools:**

- **JalView:** Feature-rich, annotations, editing
- **MEGA:** MSA + phylogenetics
- **Geneious:** Commercial, comprehensive
- **UGENE:** Free, multi-platform
- **WebLogo:** Sequence logos (conservation visualization)

**Alignment Quality Assessment:**

**Conserved Columns:**

- High conservation = likely functional
- Star (\*) in Clustal format = fully conserved
- Colon (:) = strongly similar
- Period (.) = weakly similar

**Trimming:**

Remove poorly aligned regions before phylogenetic analysis.

- **trimAl:** Automated trimming
- **Gblocks:** Conservative trimming
- Manual inspection recommended

---

### Profile Searches and Hidden Markov Models

For detecting remote homologs beyond BLAST sensitivity.

#### Position-Specific Scoring Matrices (PSSM)

**PSI-BLAST (Position-Specific Iterative BLAST):**

1. Run initial BLAST
2. Build PSSM from significant hits
3. Search database with PSSM
4. Iterate (update PSSM with new hits)
5. Converges after 3-5 iterations

**Advantages:**

- Detects remote homologs
- More sensitive than BLAST

**Risks:**

- Can drift to unrelated sequences
- Requires careful threshold setting

---

#### Hidden Markov Models (HMMs)

Probabilistic models of sequence families.

**Profile HMMs:**

- Model entire MSA
- Each position has emission probabilities
- Allow insertions and deletions
- More sophisticated than PSSM

**HMMER:**

- Build HMMs from MSA
- Search sequences against HMM
- Search HMM against database
- Very sensitive

**Pfam Database:**

- Protein family database
- Each family represented by HMM
- Used for domain identification

**Example Workflow:**

```bash
# Build HMM from alignment
hmmbuild model.hmm alignment.sto

# Search sequences with HMM
hmmsearch model.hmm sequences.fasta

# Search sequence against Pfam
hmmscan --domtblout results.txt pfam_database.hmm protein.fasta
```

**When to Use:**

- Remote homology detection
- Domain identification
- Protein family classification
- More sensitive than BLAST for divergent sequences

---

### Practical Considerations

**Choosing Alignment Method:**

- **Two similar sequences:** Needleman-Wunsch (global) or Smith-Waterman (local)
- **Database search:** BLAST, DIAMOND, or MMseqs2
- **Multiple sequences (small):** MUSCLE or MAFFT
- **Multiple sequences (large):** MAFFT
- **Remote homology:** PSI-BLAST or HMMER
- **Domain identification:** HMMER against Pfam

**Alignment Quality:**

- **Coverage:** How much of sequences aligned?
- **Identity:** Percentage of identical residues
- **Similarity:** Including conservative substitutions
- **Gaps:** Too many gaps suggest poor alignment

**Common Pitfalls:**

- Aligning non-homologous sequences (always get _some_ alignment)
- Not considering biology (alignment != homology)
- Ignoring statistical significance
- Over-interpreting low-confidence regions

---

## Sequence Motifs and Patterns

Sequence motifs are short, recurring patterns in biological sequences that often have functional or structural significance.

### Types of Motifs

**DNA Motifs:**

- **Transcription Factor Binding Sites (TFBS):** Where regulatory proteins bind
- **Promoter elements:** TATA box, CAAT box, GC box
- **Splice sites:** GT-AG rule, branch point
- **Ribosome binding sites:** Shine-Dalgarno (prokaryotes), Kozak (eukaryotes)
- **Restriction enzyme sites:** Specific recognition sequences
- **CpG islands:** Regulatory regions in mammals

**Protein Motifs:**

- **Functional domains:** Catalytic sites, binding sites
- **Structural motifs:** Helix-turn-helix, zinc finger, leucine zipper
- **Localization signals:** Nuclear localization (NLS), export signals
- **Post-translational modification sites:** Phosphorylation, glycosylation
- **Protein-protein interaction:** SH2, SH3 domains

**RNA Motifs:**

- **Stem-loops:** Secondary structure elements
- **Riboswitches:** Regulatory elements
- **Internal ribosome entry sites (IRES)**

---

### Representing Motifs

#### Consensus Sequences

Simple representation showing most common base/amino acid at each position.

**IUPAC Ambiguity Codes:**

```
DNA:
R = A or G (puRine)
Y = C or T (pYrimidine)
M = A or C (aMino)
K = G or T (Keto)
S = G or C (Strong)
W = A or T (Weak)
H = A or C or T (not G)
B = C or G or T (not A)
V = A or C or G (not T)
D = A or G or T (not C)
N = any base

Protein:
B = D or N (aspartic acid or asparagine)
Z = E or Q (glutamic acid or glutamine)
X = any amino acid
```

**Example:**

TATA box consensus: TATAWAAW (W = A or T)

**Limitations:**

- Loses information about position-specific preferences
- Doesn't capture varying conservation levels
- Binary (included or not)

---

#### Position Weight Matrices (PWM)

Quantitative representation of motif conservation at each position.

**Also Called:**

- Position-Specific Scoring Matrix (PSSM)
- Position-Specific Weight Matrix (PSWM)

**Construction:**

1. **Count matrix:** Count each nucleotide at each position from aligned sequences

```
Example: 5 sequences
ATCG
ATGG
ATCG
ATGG
ATCG

Count Matrix:
     Pos1  Pos2  Pos3  Pos4
A     5     0     0     0
C     0     0     4     0
G     0     0     1     5
T     0     5     0     0
```

2. **Frequency matrix:** Normalize by number of sequences

```
Frequency Matrix:
     Pos1  Pos2  Pos3  Pos4
A    1.0   0.0   0.0   0.0
C    0.0   0.0   0.8   0.0
G    0.0   0.0   0.2   1.0
T    0.0   1.0   0.0   0.0
```

3. **Weight matrix:** Convert to log-odds scores

```
Score = log2(observed frequency / background frequency)

If background = 0.25 for all bases:

Weight Matrix:
     Pos1  Pos2  Pos3  Pos4
A    2.0  -inf  -inf  -inf
C   -inf  -inf   1.68 -inf
G   -inf  -inf  -0.32  2.0
T   -inf   2.0  -inf  -inf
```

**Pseudocounts:**

Add small value to avoid zero frequencies (undefined log).

```
Adjusted frequency = (count + pseudocount) / (total + 4*pseudocount)
```

**Scoring Sequences:**

Sum scores across all positions.

```
Sequence: ATCG
Score = 2.0 + 2.0 + 1.68 + 2.0 = 7.68
```

**Threshold:**

Define cutoff score for motif presence.

---

#### Sequence Logos

Visual representation of motif conservation.

**Features:**

- Height of letter = information content (bits)
- Taller = more conserved
- Stack height = total information at position
- Maximum 2 bits for DNA, ~4.3 bits for protein

**Information Content:**

```
I = Σ(p × log2(p/background))

Where p = frequency of base/amino acid
```

**Example:**

```
      C
     TC    G
A    TGC   G
TTTTTCCGGGG
Pos: 1234567
```

**Tools:**

- **WebLogo:** Most popular, web-based
- **Seq2Logo:** Supports protein logos
- **ggseqlogo:** R package for customization

---

### Motif Discovery

Finding motifs in sequences without prior knowledge.

#### De Novo Motif Discovery

**Problem:**

Given set of sequences, find overrepresented patterns (motifs).

**Challenges:**

- Motif location unknown
- Motif length unknown
- May have multiple motifs
- Background noise
- Degenerate positions

---

#### Algorithms

**Expectation-Maximization (EM):**

**MEME (Multiple EM for Motif Elicitation):**

- Most widely used motif discovery tool
- Finds multiple motifs in sequences
- Handles variable motif widths
- Provides statistical significance

**Algorithm:**

1. **Initialize:** Random motif starting positions
2. **Expectation:** Estimate motif PWM from current sites
3. **Maximization:** Find best match positions for PWM
4. **Iterate:** Until convergence
5. **Repeat:** Find additional motifs

**Usage:**

```bash
meme sequences.fasta \
     -dna \
     -nmotifs 5 \
     -minw 6 \
     -maxw 20 \
     -revcomp

# Output: HTML report with motifs, logos, E-values
```

**Parameters:**

- **-nmotifs:** Number of motifs to find
- **-minw/-maxw:** Motif width range
- **-revcomp:** Search both strands
- **-mod:** Distribution model (oops, zoops, anr)

---

**Gibbs Sampling:**

Probabilistic approach using sampling.

**Tools:**

- **AlignACE:** Gibbs sampling implementation
- **MotifSampler:** Gibbs sampling with higher-order background

**Advantages:**

- Can escape local optima better than EM
- Good for weak motifs

---

**Enumeration Methods:**

**DREME (Discriminative Regular Expression Motif Elicitation):**

- Fast motif discovery
- Uses discriminative approach (positive vs negative sequences)
- Finds short, ungapped motifs
- Part of MEME Suite

**STREME:**

- Successor to DREME
- Faster
- Sensitive enrichment

---

**Word-Based Methods:**

Count k-mer frequencies, identify overrepresented words.

**Advantages:**

- Very fast
- Exact (no probabilistic sampling)

**Limitations:**

- Limited to short, exact motifs
- Doesn't handle degeneracy well

---

#### Motif Scanning

Search sequences for known motifs.

**FIMO (Find Individual Motif Occurrences):**

- Part of MEME Suite
- Scans sequences with PWM
- Provides p-values and q-values
- Fast and sensitive

**Usage:**

```bash
fimo --thresh 1e-4 motif.meme sequences.fasta

# Output: Positions, scores, p-values for matches
```

**MAST (Motif Alignment and Search Tool):**

- Searches for combinations of motifs
- Ranks sequences by motif content
- Part of MEME Suite

---

**HOMER (Hypergeometric Optimization of Motif EnRichment):**

- Designed for ChIP-seq peak sequences
- De novo discovery and known motif scanning
- Highly optimized for regulatory sequences

```bash
findMotifsGenome.pl peaks.bed genome output_dir
```

---

#### Comparative Approaches

**Phylogenetic Footprinting:**

Find conserved motifs across related species.

- Alignment of promoter regions
- Identify conserved blocks
- Higher confidence than single species

**Tools:**

- **PhyloGibbs:** Gibbs sampling with phylogeny
- **PhyMe:** Phylogenetic motif finding
- **PRIORITY:** Uses evolutionary conservation

---

### Motif Databases

**JASPAR:**

- Open-access transcription factor database
- Curated TFBS profiles
- Vertebrate focus (but includes other organisms)
- PWMs freely downloadable

**TRANSFAC:**

- Commercial transcription factor database
- Comprehensive but requires license
- Public version available (limited)

**Rfam:**

- RNA motif database
- Covariance models
- Non-coding RNA families

**PROSITE:**

- Protein motif database
- Patterns and profiles
- Functional sites and domains

**ELM (Eukaryotic Linear Motifs):**

- Short linear motifs in proteins
- Regulatory sites
- Interaction domains

---

### RNA Structure Motifs

RNA forms complex secondary and tertiary structures.

#### Secondary Structure Prediction

**Thermodynamic Methods:**

Based on free energy minimization.

**Minimum Free Energy (MFE):**

- Find structure with lowest free energy
- Tools: RNAfold (ViennaRNA package), mfold

**Partition Function:**

- Consider ensemble of structures
- Base pairing probabilities
- More realistic than single MFE structure

**Algorithms:**

- **Nussinov:** Simple, no energy model
- **Zuker:** MFE with energy parameters
- **McCaskill:** Partition function

---

#### Structure Motifs

**Common RNA Structures:**

- **Stem-loops (hairpins):** Base-paired stem + unpaired loop
- **Bulges:** Unpaired bases on one strand
- **Internal loops:** Unpaired bases on both strands
- **Multi-branch loops:** Junction of multiple stems
- **Pseudoknots:** Non-nested base pairing (complex)

**Functional RNA Motifs:**

- **Shine-Dalgarno sequence:** AGGAGGU (ribosome binding)
- **Kozak sequence:** GCCRCCATGG (translation start in eukaryotes)
- **Poly-A signal:** AAUAAA (polyadenylation)
- **Iron response element (IRE):** Regulates iron homeostasis
- **Riboswitches:** Ligand-binding regulatory elements

---

### Practical Applications

#### Promoter Prediction

**Bacterial Promoters:**

- -10 box (Pribnow box): TATAAT
- -35 box: TTGACA
- Spacing: ~17 bp

**Eukaryotic Promoters:**

- TATA box: TATAWAAW (~25-30 bp upstream of TSS)
- CAAT box: CCAAT
- GC box: GGGCGG
- More complex, often TATA-less

**Tools:**

- **PromoterInspector:** Eukaryotic promoters
- **BPROM:** Bacterial promoters
- **ProScan:** Multiple organisms

---

#### Splice Site Prediction

**Consensus Sequences:**

- **5' splice site (donor):** MAG|GURAGU (| = splice point)
- **3' splice site (acceptor):** YAG|G
- **Branch point:** YURAY (~20-50 bp upstream of 3' site)

**Tools:**

- **MaxEntScan:** Maximum entropy model
- **SpliceAI:** Deep learning-based
- **NetGene2:** Neural network

---

#### Binding Site Prediction

**Protein-DNA Binding:**

- Use PWMs from databases (JASPAR, TRANSFAC)
- Scan sequences for matches
- Consider evolutionary conservation

**Protein-Protein Interaction Sites:**

- SH2, SH3 domains
- PDZ domains
- WW domains
- Short linear motifs (SLiMs)

**Tools:**

- **ScanProsite:** Scan for PROSITE patterns
- **InterProScan:** Comprehensive domain scan

---

### Challenges and Considerations

**Motif Discovery:**

- **Weak signals:** Biological motifs often degenerate
- **Background noise:** Many sequences, few true motifs
- **Motif length:** Must specify or search range
- **Multiple motifs:** Can interfere with each other

**Statistical Significance:**

- Adjust for multiple testing
- Consider background nucleotide composition
- Validation essential (experimental or conservation)

**False Positives:**

- High-scoring matches may not be functional
- Context matters (chromatin, 3D structure)
- Validation required

**Best Practices:**

- Use multiple discovery methods
- Validate with orthogonal data
- Consider evolutionary conservation
- Experimental validation when possible
- Be skeptical of marginal predictions

---

## Phylogenetic Analysis

Phylogenetics reconstructs evolutionary relationships between organisms or genes, representing them as evolutionary trees (phylogenies).

### Why Phylogenetic Analysis?

**Applications:**

- Understand evolutionary history
- Classify organisms
- Track disease outbreaks (epidemiology)
- Guide drug/vaccine development
- Predict gene function by homology
- Conservation planning
- Forensics

### Terminology

**Basic Concepts:**

- **Phylogeny (Phylogenetic Tree):** Hypothesis of evolutionary relationships
- **Topology:** Branching pattern of tree
- **Branch:** Line connecting nodes
- **Node:** Point on tree (internal or terminal)
- **Terminal node (leaf):** Represents taxon (species, gene, sequence)
- **Internal node:** Represents common ancestor
- **Root:** Most recent common ancestor of all taxa
- **Branch length:** Evolutionary distance (time, mutations)

**Tree Types:**

- **Rooted:** Direction of evolution specified (common ancestor at root)
- **Unrooted:** Shows relationships but not evolutionary direction
- **Bifurcating:** Each node splits into two branches
- **Multifurcating (polytomy):** Node splits into >2 branches (unresolved)

**Tree Representations:**

```
Rooted tree:
        ┌─ A
    ┌───┤
    │   └─ B
────┤
    │   ┌─ C
    └───┤
        └─ D

Unrooted tree:
    A
     \
      \
       *─── C
      /
     /
    B    D
```

---

### Phylogenetic Methods

#### Distance-Based Methods

Calculate pairwise distances, construct tree from distance matrix.

**Distance Calculation:**

**p-distance:**

```
d = number of differences / sequence length
```

**Jukes-Cantor (JC69):**

Corrects for multiple substitutions at same site.

```
d = -3/4 × ln(1 - 4p/3)
```

Assumes equal substitution rates.

**Kimura 2-parameter (K2P):**

```
d = -1/2 × ln(1 - 2P - Q) - 1/4 × ln(1 - 2Q)
```

P = transitions, Q = transversions
Accounts for transition/transversion bias.

**More Complex Models:**

- **GTR (General Time Reversible):** Most general substitution model
- **HKY85:** Unequal base frequencies
- **Gamma distribution:** Rate variation across sites

---

**UPGMA (Unweighted Pair Group Method with Arithmetic Mean):**

**Algorithm:**

1. Find closest pair of taxa
2. Join them, calculate average distance to others
3. Treat as single unit
4. Repeat until all joined

**Characteristics:**

- Produces rooted tree
- Assumes molecular clock (equal rates)
- Simple and fast
- Often inaccurate (clock assumption violated)

**When to Use:**

- Very closely related sequences
- Quick approximation
- Guide tree for MSA

---

**Neighbor-Joining (NJ):**

**Algorithm:**

1. Calculate distance matrix
2. Find pair that minimizes total branch length
3. Join pair, recalculate distances
4. Repeat until tree complete

**Characteristics:**

- Produces unrooted tree
- No molecular clock assumption
- Fast (polynomial time)
- Widely used
- Reasonably accurate for close to moderate distances

**When to Use:**

- Large datasets (fast)
- First-pass analysis
- Bootstrap resampling (need speed)

**Formula:**

```
Net divergence:
r_ij = d_ij - (u_i + u_j)

where u_i = average distance from i to all other taxa
```

Join taxa i and j with minimum r_ij.

---

#### Character-Based Methods

Use sequence data directly (not distance matrix).

**Maximum Parsimony (MP):**

**Principle:** Find tree requiring fewest evolutionary changes (most parsimonious).

**Algorithm:**

1. Generate possible trees
2. Count number of changes for each tree
3. Choose tree(s) with minimum changes

**Informative Sites:**

Only parsimony-informative sites matter:

```
Site type:        Informative?
A A A A           No (constant)
A A A T           No (one variant)
A A T T           Yes (two variants, ≥2 each)
A T T C           Yes
```

**Tree Searching:**

- **Exhaustive:** All possible trees (small datasets only)
- **Branch-and-bound:** Optimal but slow
- **Heuristic:** Stepwise addition, tree bisection-reconnection (TBR)

**Characteristics:**

- No explicit evolutionary model
- Works well for closely related sequences
- Can be inconsistent (long-branch attraction)
- Multiple equally parsimonious trees common

**When to Use:**

- Closely related sequences
- Morphological data
- Few sequences

---

**Maximum Likelihood (ML):**

**Principle:** Find tree most likely to produce observed data given evolutionary model.

**Likelihood Calculation:**

```
L = P(Data | Tree, Model)
```

For each tree topology and branch lengths, calculate likelihood.

**Steps:**

1. Choose evolutionary model
2. Propose tree topology and branch lengths
3. Calculate likelihood
4. Optimize branch lengths
5. Compare trees
6. Select tree with highest likelihood

**Evolutionary Models:**

Must specify substitution model:

- **JC69:** Equal rates, equal frequencies
- **K80 (K2P):** Different transition/transversion rates
- **HKY85:** Unequal base frequencies
- **GTR:** Most general (6 substitution rates)

Plus additional parameters:

- **+I:** Proportion of invariant sites
- **+G:** Gamma-distributed rate variation (α parameter)
- **+F:** Empirical base frequencies

Example: **GTR+G+I** (complex, realistic model)

**Model Selection:**

Use information criteria to choose best model:

- **AIC (Akaike Information Criterion)**
- **BIC (Bayesian Information Criterion)**
- **Likelihood Ratio Test (LRT)**

Tools: jModelTest (DNA), ProtTest (protein)

**Characteristics:**

- Statistical framework
- Handles complex evolutionary scenarios
- Computationally intensive
- Most widely used for publication-quality trees

**When to Use:**

- Most phylogenetic analyses
- Need statistical support
- Sufficient computational resources

---

**Bayesian Inference:**

**Principle:** Calculate posterior probability distribution of trees.

**Bayes' Theorem:**

```
P(Tree | Data) = P(Data | Tree) × P(Tree) / P(Data)

Posterior = Likelihood × Prior / Evidence
```

**Markov Chain Monte Carlo (MCMC):**

Sample trees proportionally to their probability:

1. Start with random tree
2. Propose modification
3. Accept/reject based on posterior probability
4. Repeat millions of times
5. Collect samples after burn-in
6. Summarize (consensus tree, posterior probabilities)

**Characteristics:**

- Provides posterior probabilities (branch support)
- Requires prior distributions
- Very computationally intensive
- Multiple chains for convergence checking

**When to Use:**

- Publication-quality trees
- Need uncertainty estimates
- Complex models
- Have computational resources and time

---

### Phylogenetic Software

| Software     | Method           | Speed     | Features                                |
| ------------ | ---------------- | --------- | --------------------------------------- |
| **MEGA**     | Distance, MP, ML | Medium    | GUI, user-friendly, visualization       |
| **RAxML**    | ML               | Fast      | Large datasets, bootstrap, multicore    |
| **IQ-TREE**  | ML               | Very Fast | Model selection, ultrafast bootstrap    |
| **PhyML**    | ML               | Fast      | Good for moderate datasets              |
| **MrBayes**  | Bayesian         | Slow      | Posterior probabilities, complex models |
| **BEAST**    | Bayesian         | Slow      | Time-calibrated trees, molecular dating |
| **FastTree** | ML (approximate) | Very Fast | Huge datasets (10,000+ taxa)            |
| **TNT**      | MP               | Fast      | Large parsimony analyses                |

**Recommendations:**

- **Quick analysis:** FastTree, neighbor-joining
- **Standard ML:** IQ-TREE or RAxML
- **Bayesian:** MrBayes or BEAST
- **Beginners:** MEGA (GUI)
- **Large datasets:** IQ-TREE or FastTree

---

### Tree Evaluation and Support

**Bootstrap Analysis:**

Assess confidence in tree topology.

**Procedure:**

1. Resample alignment with replacement
2. Reconstruct tree from resampled data
3. Repeat 100-1000 times
4. Count how often each clade appears

**Interpretation:**

- **>95%:** Strong support
- **75-95%:** Moderate support
- **<70%:** Weak support
- **<50%:** Not supported

**Types:**

- **Standard bootstrap:** Slow (ML, MP)
- **Ultrafast bootstrap (UFBoot):** Much faster, IQ-TREE
- **Approximate likelihood ratio test (aLRT):** Alternative to bootstrap

**Bayesian Posterior Probabilities:**

- **>0.95:** Strong support
- Generally higher than bootstrap values
- Can be overconfident

---

### Tree Rooting

Unrooted trees show relationships but not evolutionary direction.

**Methods:**

**Outgroup Rooting:**

- Include distantly related taxon
- Root on branch leading to outgroup
- Most common method

**Midpoint Rooting:**

- Place root at midpoint of longest path
- Assumes molecular clock
- Use when no outgroup available

**Molecular Clock:**

- Use calibration points (fossils, geological events)
- Estimate divergence times
- Requires special software (BEAST, r8s)

---

### Tree Visualization and Manipulation

**Visualization Tools:**

**FigTree:**

- Popular, easy to use
- Annotate, color, format trees
- Export publication-quality figures

**iTOL (Interactive Tree of Life):**

- Web-based
- Rich annotations
- Handles very large trees

**Dendroscope:**

- Compare multiple trees
- Tanglegrams

**ggtree (R):**

- Programmatic tree visualization
- Highly customizable
- Integrated with R/Bioconductor

**ETE Toolkit (Python):**

- Python library
- Programmatic manipulation
- Custom visualizations

---

**Tree Manipulation:**

**Newick Format:**

Standard text representation of trees.

```
((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);

Structure:
(taxon:branch_length, taxon:branch_length):branch_length;
```

**Nexus Format:**

More complex, includes alignment and tree.

```
#NEXUS
BEGIN TREES;
    TREE tree1 = ((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);
END;
```

---

### Common Phylogenetic Analyses

#### Molecular Clock Analysis

Estimate divergence times.

**Assumptions:**

- Constant rate of evolution (strict clock)
- Or relaxed clock (rate varies among branches)

**Calibration:**

- Fossil evidence
- Geological events
- Known divergence dates

**Tools:**

- **BEAST:** Bayesian, relaxed clock
- **r8s:** ML and penalized likelihood
- **TimeTree:** Database of divergence times

---

#### Ancestral State Reconstruction

Infer characteristics of ancestral organisms.

**Methods:**

- **Parsimony:** Minimize changes
- **ML:** Maximum likelihood estimation
- **Bayesian:** Posterior probabilities

**Applications:**

- Trait evolution
- Biogeography
- Protein ancestral sequences

**Tools:**

- **PAML:** ML ancestral reconstruction
- **Mesquite:** Character evolution
- **ape (R):** Ancestral state estimation

---

#### Species Tree vs Gene Tree

**Gene Tree:** Phylogeny of specific gene
**Species Tree:** Phylogeny of organisms

**Why Different?**

- **Incomplete lineage sorting:** Polymorphism maintained across speciation
- **Horizontal gene transfer:** Especially in prokaryotes
- **Gene duplication/loss:** Paralogs vs orthologs
- **Recombination:** Different parts of genome have different histories

**Species Tree Methods:**

- **Concatenation:** Combine multiple genes (assumes same tree)
- **Coalescent methods:** Account for gene tree discordance
  - ASTRAL: Fast, accurate
  - SVDquartets: SNP-based
  - \*BEAST: Bayesian co-estimation

---

### Phylogenomics

Phylogenetics using genome-scale data.

**Approaches:**

- **Core genome:** Genes present in all taxa
- **Pangenome:** All genes in any taxa
- **SNPs:** Single nucleotide polymorphisms
- **Whole genome alignment:** Entire genomes

**Advantages:**

- More data = better resolution
- Can resolve difficult relationships
- Detect hybridization and horizontal transfer

**Challenges:**

- Computational demands
- Gene tree discordance
- Model violations
- Alignment ambiguity

**Tools:**

- **OrthoFinder:** Identify orthologs
- **BUSCO:** Single-copy orthologs
- **RAxML-NG:** Large-scale ML
- **IQ-TREE:** Phylogenomic analysis

---

### Practical Workflow

**Standard Phylogenetic Analysis:**

1. **Collect sequences:**

   - From databases or experiments
   - Include outgroup

2. **Multiple sequence alignment:**

   - MAFFT, MUSCLE, or Clustal Omega
   - Check alignment quality

3. **Trim alignment:**

   - Remove poorly aligned regions
   - trimAl or Gblocks

4. **Model selection:**

   - jModelTest, ModelFinder (in IQ-TREE)

5. **Tree reconstruction:**

   - ML (IQ-TREE, RAxML) or Bayesian (MrBayes)
   - Include bootstrap or posterior probabilities

6. **Root tree:**

   - Outgroup or midpoint

7. **Visualize:**

   - FigTree, iTOL, or ggtree

8. **Interpret:**
   - Assess support values
   - Consider biological context

---

### Common Pitfalls

**Poor Alignment:**

- Garbage in, garbage out
- Inspect and manually refine if needed
- Exclude ambiguous regions

**Wrong Model:**

- Simple models can mislead
- Always perform model selection
- More parameters ≠ always better (overfitting)

**Insufficient Data:**

- Short sequences give uncertain trees
- Consider more genes or longer regions

**Long-Branch Attraction:**

- Fast-evolving lineages artifactually group together
- Affects parsimony especially
- Use ML with appropriate model

**Taxon Sampling:**

- Missing intermediate taxa problematic
- More taxa often better than more data per taxon

**Overinterpreting Support:**

- High bootstrap ≠ correct tree
- Consider alternative topologies
- Biological validation important

---

## Genomic Databases and Resources

Biological databases are essential infrastructure for bioinformatics, providing curated collections of sequences, structures, functions, and annotations.

### Primary Sequence Databases

These contain raw sequence data submitted by researchers.

#### GenBank (NCBI)

**Description:**

- Comprehensive nucleotide sequence database
- Part of International Nucleotide Sequence Database Collaboration (INSDC)
- Over 250 million sequences
- Public submissions from worldwide researchers

**Components:**

- **GenBank:** Main annotated sequence database
- **dbEST:** Expressed Sequence Tags
- **dbGSS:** Genome Survey Sequences
- **SRA (Sequence Read Archive):** Raw sequencing reads

**Access:**

- Web: https://www.ncbi.nlm.nih.gov/genbank/
- FTP: Bulk downloads
- E-utilities API: Programmatic access
- Entrez: Search interface

**Accession Numbers:**

```
Format examples:
U12345    - Single sequence (1-5 letters + 5-8 digits)
AABCD01000000 - WGS project
NM_123456 - RefSeq (curated)
```

**File Formats:**

- GenBank flat file format (.gb, .gbk)
- FASTA (.fasta, .fa)
- XML

---

#### ENA (European Nucleotide Archive)

**Description:**

- European counterpart to GenBank
- Part of INSDC
- Hosted by EMBL-EBI
- Synchronized with GenBank and DDBJ

**Unique Features:**

- Stronger in European submissions
- Read files and assemblies
- Excellent metadata

**Access:**

- Web: https://www.ebi.ac.uk/ena
- FTP and API available

---

#### DDBJ (DNA Data Bank of Japan)

**Description:**

- Japanese member of INSDC
- Synchronized with GenBank and ENA
- Strong in Asian submissions

---

### Curated and Reference Databases

Higher quality, manually curated subset of sequences.

#### RefSeq (Reference Sequence Database)

**Description:**

- NCBI's curated non-redundant sequence database
- One representative sequence per gene/genome
- Manually reviewed annotations
- High quality

**Accession Prefixes:**

```
NM_  - mRNA (curated)
NP_  - Protein (curated)
NR_  - Non-coding RNA (curated)
NC_  - Complete chromosome/genome (curated)
XM_  - mRNA (predicted)
XP_  - Protein (predicted)
```

**Advantages:**

- High quality annotations
- Stable identifiers
- Standard reference
- Well-curated

**Use Cases:**

- Reference genomes
- Annotation transfer
- Clinical variant interpretation

---

#### Ensembl

**Description:**

- Genome database for vertebrates and model organisms
- Comprehensive annotations
- Comparative genomics
- Extensive visualization tools

**Content:**

- Genome assemblies
- Gene predictions
- Regulatory features
- Variation data
- Comparative genomics

**Tools:**

- **Genome browser:** Visualize genes and features
- **BioMart:** Data mining and bulk retrieval
- **Variant Effect Predictor (VEP):** Predict variant consequences
- **REST API:** Programmatic access

**Bacteria:**

- **Ensembl Bacteria:** Bacterial and archaeal genomes
- **Ensembl Fungi, Protists, Plants, Metazoa:** Other domains

**Access:**

- Web: https://www.ensembl.org
- FTP: ftp://ftp.ensembl.org
- BioMart and REST API

---

#### UCSC Genome Browser

**Description:**

- Comprehensive genome annotation database
- Excellent visualization
- Multiple species (focus on human)
- Rich track data

**Features:**

- Gene annotations
- Comparative genomics (conservation tracks)
- Regulation data (ENCODE)
- Variation databases
- Custom tracks

**Tools:**

- **BLAT:** Fast sequence alignment
- **Table Browser:** Query and download annotations
- **Genome Browser:** Visualization
- **LiftOver:** Convert coordinates between assemblies

**Access:**

- Web: https://genome.ucsc.edu
- MySQL database: Direct access
- FTP downloads

---

### Protein Databases

#### UniProt (Universal Protein Resource)

**Description:**

- Comprehensive protein sequence and functional information
- Most important protein database
- Combines Swiss-Prot (curated) and TrEMBL (automated)

**Components:**

**Swiss-Prot:**

- Manually annotated and reviewed
- High quality
- ~560,000 entries
- Functional information, PTMs, variants

**TrEMBL:**

- Automatically annotated
- Computationally analyzed
- > 200 million entries
- Not manually reviewed

**UniRef:**

- Clustered sequences (UniRef100, 90, 50)
- Reduce redundancy
- Faster searches

**UniParc:**

- Non-redundant archive
- All protein sequences
- No annotation

**Accession Numbers:**

```
P12345    - Swiss-Prot
Q8N123    - Swiss-Prot
A0A123ABC4 - TrEMBL
```

**Information Provided:**

- Protein function
- Catalytic activity
- Cofactors and PTMs
- Subcellular location
- Disease associations
- 3D structure links
- Sequence features
- Publications

**Access:**

- Web: https://www.uniprot.org
- FTP: Bulk downloads
- REST API
- SPARQL endpoint

---

#### PDB (Protein Data Bank)

**Description:**

- Repository of 3D structural data
- X-ray crystallography, NMR, cryo-EM
- > 200,000 structures

**File Formats:**

- PDB: Legacy text format
- mmCIF: Modern structured format
- PDBML: XML format

**Related Resources:**

- **RCSB PDB:** US portal
- **PDBe:** European portal
- **PDBj:** Japanese portal
- **AlphaFold DB:** Predicted structures (millions)

**Access:**

- Web: https://www.rcsb.org
- FTP downloads
- REST API

---

### Functional Databases

#### KEGG (Kyoto Encyclopedia of Genes and Genomes)

**Description:**

- Metabolic pathways and biochemical reactions
- Functional genomics resource
- Widely used for pathway analysis

**Databases:**

- **PATHWAY:** Metabolic, signaling, disease pathways
- **GENES:** Gene catalogs
- **ORTHOLOGY (KO):** Functional orthologs
- **COMPOUND:** Chemical substances
- **REACTION:** Biochemical reactions
- **ENZYME:** Enzyme nomenclature (EC numbers)
- **BRITE:** Hierarchical classifications
- **MODULE:** Functional units in pathways

**Identifiers:**

```
K00001    - KO (KEGG Orthology)
ko00010   - Glycolysis pathway
C00002    - ATP
R00001    - Reaction
M00001    - Module
```

**Use Cases:**

- Pathway enrichment analysis
- Metabolic reconstruction
- Functional annotation
- Comparative genomics

**Access:**

- Web: https://www.kegg.jp
- REST API
- Commercial license required for some uses

---

#### GO (Gene Ontology)

**Description:**

- Controlled vocabulary for gene function
- Three ontologies (BP, MF, CC)
- Species-independent
- Hierarchical structure

**Three Ontologies:**

1. **Biological Process (BP):** "What" gene product does

   - Example: GO:0006260 (DNA replication)

2. **Molecular Function (MF):** Biochemical activity

   - Example: GO:0003677 (DNA binding)

3. **Cellular Component (CC):** Where gene product is located
   - Example: GO:0005737 (cytoplasm)

**Structure:**

- Directed acyclic graph (DAG)
- Terms connected by relationships:
  - is_a (stricter relationship)
  - part_of
  - regulates

**GO Annotations:**

- Link gene products to GO terms
- Evidence codes indicate confidence
- Species-specific databases (GOA)

**Use Cases:**

- Functional annotation
- Gene set enrichment analysis
- Compare gene lists
- Predict function

**Access:**

- Web: http://geneontology.org
- Downloads: OBO and OWL formats
- QuickGO: Search and browse

---

#### Pfam (Protein Families)

**Description:**

- Database of protein families and domains
- Profile Hidden Markov Models (HMMs)
- ~20,000 families
- High quality, manually curated

**Content:**

- Protein domain definitions
- Multiple sequence alignments
- HMM profiles
- Functional annotations

**Use Cases:**

- Domain identification
- Protein classification
- Functional prediction
- Evolutionary analysis

**Access:**

- Web: http://pfam.xfam.org
- InterPro integration
- HMM files downloadable

---

#### InterPro

**Description:**

- Integrates multiple protein signature databases
- Unified protein family resource
- Functional and structural classification

**Integrated Databases:**

- Pfam
- PROSITE
- PRINTS
- ProDom
- SMART
- TIGRFAMs
- PIRSF
- SUPERFAMILY
- Gene3D
- PANTHER

**Advantages:**

- One search, multiple databases
- Comprehensive coverage
- Hierarchical organization
- GO term mapping

**Tools:**

- **InterProScan:** Scan sequences against all databases
- Command-line and web interface

**Access:**

- Web: https://www.ebi.ac.uk/interpro/
- InterProScan standalone

---

### Genomic Variation Databases

#### dbSNP

**Description:**

- Single Nucleotide Polymorphism database (NCBI)
- Short variants (SNPs, indels, microsatellites)
- Millions of human variants

**Identifiers:**

```
rs123456 - RefSNP ID
```

**Content:**

- Variant position and alleles
- Population frequencies
- Clinical significance
- Functional consequences

---

#### ClinVar

**Description:**

- Clinical significance of genetic variants
- Links genotypes to phenotypes
- Clinical interpretation

**Classifications:**

- Pathogenic
- Likely pathogenic
- Uncertain significance (VUS)
- Likely benign
- Benign

**Use Cases:**

- Clinical diagnostics
- Variant interpretation
- Research on disease variants

---

#### gnomAD (Genome Aggregation Database)

**Description:**

- Population genetics resource
- Allele frequencies from >140,000 individuals
- Diverse populations

**Use Cases:**

- Filter common variants
- Estimate pathogenicity
- Population genetics

**Access:**

- Web browser: https://gnomad.broadinstitute.org
- Downloads available

---

### Taxonomic Databases

#### NCBI Taxonomy

**Description:**

- Hierarchical classification of organisms
- Taxonomic IDs (TaxIDs)
- Links to sequence data

**Structure:**

```
Superkingdom → Phylum → Class → Order → Family → Genus → Species
```

**TaxIDs:**

```
9606  - Homo sapiens
562   - Escherichia coli
```

**Access:**

- Web: https://www.ncbi.nlm.nih.gov/taxonomy
- E-utilities API
- Taxonomy browser

---

#### GTDB (Genome Taxonomy Database)

**Description:**

- Phylogeny-based taxonomy
- Bacterial and archaeal genomes
- Alternative to NCBI taxonomy
- Based on genome trees

**Advantages:**

- Consistent with phylogeny
- Regular updates
- Standardized nomenclature
- Better representation of microbial diversity

**Use Cases:**

- Microbial taxonomy
- Metagenomics classification
- Phylogenetic analysis

**Access:**

- Web: https://gtdb.ecogenomic.org
- Downloads: Trees, metadata, taxonomies

---

### Specialized Databases

#### MicrobeDB/IMG (Integrated Microbial Genomes)

- Microbial genomes
- Comparative analysis tools
- Functional annotations

#### PATRIC (Pathosystems Resource Integration Center)

- Bacterial pathogens
- Comparative genomics
- Antibiotic resistance data

#### FungiDB, PlasmoDB, VectorBase

- Eukaryotic pathogen databases
- Part of VEuPathDB consortium
- Integrated tools and genomics

#### miRBase

- microRNA database
- Sequences and annotations
- Nomenclature authority

#### Rfam

- Non-coding RNA families
- Covariance models
- Functional RNAs

---

### Metagenomic Databases

#### MGnify (formerly EBI Metagenomics)

- Metagenomic and metatranscriptomic data
- Automated analysis pipelines
- Taxonomic and functional analysis

#### NCBI SRA (Sequence Read Archive)

- Raw sequencing data
- Includes metagenomics
- Largest repository

#### IMG/M (Integrated Microbial Genomes & Microbiomes)

- Metagenome analysis
- Comparative tools
- Functional profiling

---

### Database Search and Retrieval

#### NCBI Entrez

**Description:**

- Unified search system across NCBI databases
- Links between related records
- PubMed, GenBank, Protein, etc.

**Tools:**

- Web interface
- E-utilities (programmatic access)

**E-utilities:**

```bash
# Example: Fetch sequence by accession
esearch -db nucleotide -query "NM_000546" | efetch -format fasta

# Search and count
esearch -db protein -query "human[ORGN] AND kinase[TITLE]" | efetch -format fasta
```

---

#### SRA Toolkit

**Description:**

- Access SRA data
- Convert SRA format to FASTQ
- Download and process sequencing data

**Tools:**

```bash
# Download and convert to FASTQ
fastq-dump SRR123456

# Paired-end reads
fastq-dump --split-files SRR123456

# Compressed output
fastq-dump --gzip SRR123456
```

---

### Data Submission

#### Why Submit?

- Publication requirement
- Data sharing (reproducibility)
- Community resource
- Citation credit

#### Major Repositories:

**Sequences:**

- GenBank/ENA/DDBJ (INSDC)
- Submit to one, synchronized to all

**Raw Reads:**

- SRA/ENA/DRA
- FASTQ files

**Assemblies:**

- GenBank (WGS)
- ENA

**Proteins:**

- UniProt (via proteome submissions)

**Structures:**

- PDB (wwPDB)

**Functional Genomics:**

- GEO (Gene Expression Omnibus) - NCBI
- ArrayExpress - EBI

**Submission Tools:**

- **GenBank:** BankIt, Sequin, tbl2asn
- **SRA:** SRA submission portal
- **ENA:** Webin
- **Programmatic:** APIs available

---

### Best Practices for Database Use

**Search Strategies:**

- Use specific keywords and filters
- Boolean operators (AND, OR, NOT)
- Field-specific searches
- Refine iteratively

**Data Quality:**

- Check sequence quality and completeness
- Verify annotations
- Check publication date (recent = better)
- Prefer curated over automated

**Citations:**

- Cite database and version
- Include accession numbers
- Note access date (for web resources)

**Programmatic Access:**

- APIs more reliable than web scraping
- Respect rate limits
- Cache results
- Check terms of use

---

## File Formats in Bioinformatics

Bioinformatics uses many standardized file formats for storing and exchanging biological data. Understanding these formats is essential for data analysis.

### Sequence Formats

#### FASTA Format

**Description:**

- Simple text format for nucleotide or protein sequences
- Most widely used sequence format
- Human-readable

**Structure:**

```
>sequence_identifier description
SEQUENCE
DATA
ON
MULTIPLE
LINES
```

**Example:**

```
>NM_000546.6 Homo sapiens tumor protein p53 (TP53), transcript variant 1, mRNA
GATGGGATTGGGGTTTTCCCCTCCCATGTGCTCAAGACTGGCGCTAAAAGTTTTGAGCT
TCTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACT
TTGCGTTCGGGCTGGGAGCGTGCTTTCCACGACGGTGACACGCTTCCCTGGATTGGCAG
```

**Header Line:**

- Starts with `>`
- First word (no spaces) = identifier
- Rest = description (optional)

**Sequence Lines:**

- Can be wrapped at any length (typically 60-80 characters)
- No line length requirement
- Case usually irrelevant (but lowercase often indicates masked repeats)

**Parsing:**

```python
# Python example
from Bio import SeqIO

for record in SeqIO.parse("sequences.fasta", "fasta"):
    print(f"ID: {record.id}")
    print(f"Length: {len(record.seq)}")
    print(f"Sequence: {record.seq[:50]}...")
```

**Multi-FASTA:**

Multiple sequences in one file, each with `>` header.

---

#### FASTQ Format

**Description:**

- FASTA with quality scores
- Standard for sequencing data (Illumina, Nanopore, PacBio)
- 4 lines per sequence

**Structure:**

```
@sequence_identifier description
SEQUENCE
+optional_identifier (usually empty or same as line 1)
QUALITY_SCORES
```

**Example:**

```
@SRR123456.1 HWI-ST:123:C0MHFACXX:1:1101:1234:2000 length=100
GATTTGGGGTTCAGCCAAGGATGTTTGATGGGATTGGGGTTTTCCCCTCCCATGTGCTC
+
IIIIIIIIHHHHHHHHHHHHHHHHHHGGGGGGGGGGGGGGGGGFFFFFFFFFFFFFCCC
```

**Quality Encoding (Phred Scores):**

ASCII character encodes quality score:

- **Phred+33 (Sanger, Illumina 1.8+):** Most common

  - ASCII 33-126 → Q0-93
  - `!` = Q0, `"` = Q1, ..., `I` = Q40

- **Phred+64 (Old Illumina):** Deprecated
  - ASCII 64-126 → Q0-62

**Quality Score Interpretation:**

```
Q = -10 × log10(P)

Q10 = 10% error rate (90% accuracy)
Q20 = 1% error rate (99% accuracy)
Q30 = 0.1% error rate (99.9% accuracy)
Q40 = 0.01% error rate (99.99% accuracy)
```

**Common Quality Thresholds:**

- Q20: Acceptable for most purposes
- Q30: High quality (Illumina target)
- Q40: Very high quality

**Compressed FASTQ:**

Usually gzipped (`.fastq.gz`) to save space:

```bash
# View compressed FASTQ
zcat file.fastq.gz | head -n 4

# Convert to FASTA
zcat file.fastq.gz | paste - - - - | cut -f 1,2 | sed 's/@/>/' > file.fasta
```

---

### Alignment Formats

#### SAM (Sequence Alignment/Map)

**Description:**

- Standard format for sequence alignments
- Text-based
- Maps reads to reference genome

**Structure:**

1. **Header section:** Lines starting with `@`
2. **Alignment section:** One line per aligned read

**Header Lines:**

```
@HD VN:1.6 SO:coordinate
@SQ SN:chr1 LN:248956422
@RG ID:sample1 SM:sample1 PL:ILLUMINA
@PG ID:bwa PN:bwa VN:0.7.17
```

- `@HD`: Header line
- `@SQ`: Reference sequence
- `@RG`: Read group
- `@PG`: Program

**Alignment Fields (11 mandatory):**

```
QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAGS]
```

| Field | Description               | Example  |
| ----- | ------------------------- | -------- |
| QNAME | Query name                | read_001 |
| FLAG  | Bitwise flags             | 163      |
| RNAME | Reference name            | chr1     |
| POS   | 1-based leftmost position | 12345    |
| MAPQ  | Mapping quality           | 60       |
| CIGAR | Alignment descriptor      | 75M      |
| RNEXT | Reference of mate         | =        |
| PNEXT | Position of mate          | 12500    |
| TLEN  | Template length           | 230      |
| SEQ   | Sequence                  | ACGT...  |
| QUAL  | Quality string            | IIII...  |

**FLAG Field (Bit Flags):**

```
1    = read paired
2    = properly aligned
4    = read unmapped
8    = mate unmapped
16   = reverse strand
32   = mate reverse strand
64   = first in pair
128  = second in pair
256  = secondary alignment
512  = fails QC
1024 = PCR duplicate
2048 = supplementary alignment
```

Example: FLAG=163 = 1+2+32+128 (paired, proper, mate reverse, second read)

**CIGAR String:**

Describes alignment: Matches, insertions, deletions.

- `M` = Match/mismatch
- `I` = Insertion (to reference)
- `D` = Deletion (from reference)
- `N` = Skipped region (intron)
- `S` = Soft clipping
- `H` = Hard clipping
- `P` = Padding
- `=` = Exact match
- `X` = Mismatch

Example: `70M5I25M` = 70 matches, 5 insertion, 25 matches

**Example SAM Line:**

```
read1 163 chr1 100 60 75M = 250 225 ACGTACGTACGT... IIIIIIIIIIII... NM:i:0 MD:Z:75 AS:i:75
```

---

#### BAM (Binary Alignment/Map)

**Description:**

- Binary version of SAM
- Compressed (smaller file size)
- Faster to process
- Requires special tools to view

**Advantages:**

- ~5-10× smaller than SAM
- Indexed for fast random access
- Sorted for efficient queries

**Tools:**

```bash
# View BAM file (samtools required)
samtools view file.bam | less

# BAM to SAM
samtools view -h file.bam > file.sam

# SAM to BAM
samtools view -bS file.sam > file.bam

# Sort BAM
samtools sort file.bam -o file.sorted.bam

# Index BAM (creates .bai file)
samtools index file.sorted.bam

# View specific region
samtools view file.bam chr1:1000-2000

# Basic statistics
samtools flagstat file.bam
```

---

#### CRAM Format

**Description:**

- Compressed BAM
- Reference-based compression
- ~50% smaller than BAM
- Requires reference genome

**When to Use:**

- Long-term storage
- Archival purposes
- Bandwidth-limited transfers

**Conversion:**

```bash
# BAM to CRAM
samtools view -C -T reference.fa input.bam > output.cram

# CRAM to BAM
samtools view -b -T reference.fa input.cram > output.bam
```

---

### Variant Call Format (VCF)

**Description:**

- Standard for storing genetic variations
- SNPs, indels, structural variants
- Text-based

**Structure:**

1. **Meta-information lines:** Start with `##`
2. **Header line:** Starts with `#CHROM`
3. **Data lines:** One per variant

**Header:**

```
##fileformat=VCFv4.2
##reference=hg38.fa
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT Sample1
```

**Data Fields:**

| Field          | Description              | Example      |
| -------------- | ------------------------ | ------------ |
| CHROM          | Chromosome               | chr1         |
| POS            | Position (1-based)       | 12345        |
| ID             | Variant ID (e.g., dbSNP) | rs123456     |
| REF            | Reference allele         | A            |
| ALT            | Alternate allele(s)      | G            |
| QUAL           | Quality score            | 100.5        |
| FILTER         | Filter status            | PASS         |
| INFO           | Additional information   | DP=50;AF=0.5 |
| FORMAT         | Genotype format          | GT:DP:GQ     |
| Sample columns | Per-sample data          | 0/1:30:99    |

**Genotype (GT) Field:**

```
0/0 = Homozygous reference
0/1 = Heterozygous
1/1 = Homozygous alternate
./. = Missing
1|0 = Phased heterozygous (first haplotype has alt)
```

**Example VCF Line:**

```
chr1 12345 rs123 A G 100 PASS DP=50;AF=0.5 GT:DP:GQ 0/1:30:99
```

**BCF Format:**

- Binary VCF
- Smaller, faster
- Use bcftools to manipulate

```bash
# VCF to BCF
bcftools view -O b file.vcf > file.bcf

# BCF to VCF
bcftools view file.bcf > file.vcf

# Index VCF
bcftools index file.vcf.gz

# Query region
bcftools view file.vcf.gz chr1:1000-2000
```

---

### Annotation Formats

#### GFF (General Feature Format)

**Description:**

- Tab-delimited format for genomic features
- Genes, exons, CDS, etc.

**GFF3 (Current Version):**

**9 Fields:**

```
seqid source type start end score strand frame attributes
```

**Example:**

```
chr1 MAKER gene 1000 2000 . + . ID=gene001;Name=ABC1
chr1 MAKER mRNA 1000 2000 . + . ID=mRNA001;Parent=gene001
chr1 MAKER exon 1000 1200 . + . ID=exon001;Parent=mRNA001
chr1 MAKER CDS  1050 1200 . + 0 ID=cds001;Parent=mRNA001
chr1 MAKER exon 1800 2000 . + . ID=exon002;Parent=mRNA001
chr1 MAKER CDS  1800 1950 . + 0 ID=cds002;Parent=mRNA001
```

**Key Features:**

- **1-based, inclusive** coordinates
- **Hierarchical relationships:** Parent-child via ID/Parent
- **Strand:** + (forward), - (reverse), . (not stranded)
- **Frame:** 0, 1, or 2 for CDS (reading frame offset)

---

#### GTF (Gene Transfer Format)

**Description:**

- Similar to GFF2
- Stricter format
- Used by many RNA-seq tools (Cufflinks, StringTie)

**Format:**

Similar to GFF but with specific attribute structure:

```
chr1 StringTie transcript 1000 2000 . + . gene_id "GENE001"; transcript_id "TX001";
chr1 StringTie exon      1000 1200 . + . gene_id "GENE001"; transcript_id "TX001";
```

**Differences from GFF3:**

- Uses `gene_id` and `transcript_id` (not ID/Parent)
- Less flexible attribute format
- Some tools require GTF, not GFF3

---

#### BED (Browser Extensible Data)

**Description:**

- Simple tab-delimited format
- Genomic regions
- Used for visualization and region operations

**Minimum 3 Fields:**

```
chrom chromStart chromEnd
```

**Extended BED (12 fields):**

```
chr1 1000 2000 feature1 100 + 1000 2000 255,0,0 2 200,200 0,800
```

| Column | Field          | Description             |
| ------ | -------------- | ----------------------- |
| 1      | chrom          | Chromosome              |
| 2      | chromStart     | Start (0-based)         |
| 3      | chromEnd       | End (exclusive)         |
| 4      | name           | Feature name            |
| 5      | score          | 0-1000                  |
| 6      | strand         | +, -, or .              |
| 7-8    | thickStart/End | Drawing thick part      |
| 9      | itemRgb        | RGB color               |
| 10     | blockCount     | Number of blocks        |
| 11     | blockSizes     | Block sizes             |
| 12     | blockStarts    | Block starts (relative) |

**Important:**

- **0-based, half-open** coordinates (different from GFF!)
- chr1:0-1000 = first 1000 bases
- Makes interval arithmetic easier

**Tools:**

```bash
# BEDTools: powerful suite for BED operations
bedtools intersect -a file1.bed -b file2.bed
bedtools merge -i file.bed
bedtools coverage -a genes.bed -b reads.bam
```

---

### Assembly and Graph Formats

#### FASTA (for Assemblies)

Contigs/scaffolds as FASTA sequences.

**Example:**

```
>contig_1 length=5000 coverage=50.0
ACGTACGTACGTACGT...
>contig_2 length=3000 coverage=45.0
TGCATGCATGCATGCA...
```

---

#### AGP (A Golden Path)

**Description:**

- Describes assembly structure
- How contigs are arranged into scaffolds/chromosomes
- Tab-delimited

**Example:**

```
chr1 1    1000  1 W contig1 1 1000 +
chr1 1001 1100  2 N 100     scaffold yes
chr1 1101 2000  3 W contig2 1 900  -
```

---

#### GFA (Graphical Fragment Assembly)

**Description:**

- Represents assembly graph
- Used by modern assemblers (SPAdes, Flye, Canu, etc.)
- Stores sequences and connections

**Line Types:**

- `S`: Segment (sequence)
- `L`: Link (connection)
- `P`: Path
- `C`: Containment

**Example:**

```
S segment1 ACGTACGT
S segment2 TGCATGCA
L segment1 + segment2 + 4M
P path1 segment1+,segment2+ *
```

**Relevance:**

- **GraphBin** works with assembly graphs
- GFA is key format for metagenome assembly

---

### Phylogenetic Formats

#### Newick Format

**Description:**

- Standard tree representation
- Compact text format

**Structure:**

```
((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);
```

- Parentheses = clades
- Commas = separate taxa/clades
- Colons = branch lengths
- Semicolon = end

---

#### NEXUS Format

**Description:**

- More complex than Newick
- Can store alignment + tree + metadata

**Structure:**

```
#NEXUS
BEGIN DATA;
    DIMENSIONS NTAX=4 NCHAR=100;
    FORMAT DATATYPE=DNA;
    MATRIX
    Seq1 ACGT...
    Seq2 ACGT...
    ;
END;
BEGIN TREES;
    TREE tree1 = ((A,B),(C,D));
END;
```

---

### Specialized Formats

#### HMM Formats

**Stockholm Format:**

- Multiple sequence alignment
- Rich annotation
- Used by Pfam, Rfam

**HMMER Format:**

- Profile HMM definition
- Binary and text versions

---

#### Protein Structure Formats

**PDB Format:**

- Atomic coordinates
- Legacy text format
- 80-character lines

**mmCIF Format:**

- Modern structured format
- Replaces PDB
- No column limitations

---

### Format Conversion Tools

**SeqKit:**

```bash
# FASTQ to FASTA
seqkit fq2fa reads.fastq.gz -o sequences.fasta

# Statistics
seqkit stats *.fasta

# Filter by length
seqkit seq -m 1000 sequences.fasta
```

**SAMtools/BCFtools:**

- SAM ↔ BAM ↔ CRAM
- VCF ↔ BCF

**BEDTools:**

- BED format operations
- Convert to other formats

**GFFREAD (part of Cufflinks):**

```bash
# GFF3 to GTF
gffread annotation.gff3 -T -o annotation.gtf

# Extract sequences
gffread -w transcripts.fa -g genome.fa annotation.gff3
```

---

### Best Practices

**Choosing Formats:**

- **Storage:** Use compressed binary (BAM, BCF)
- **Processing:** Binary formats faster
- **Sharing:** Text formats more portable
- **Archives:** CRAM for long-term

**Validation:**

- Check file integrity
- Validate format compliance
- Use format-specific validators

**Compression:**

- GZIP for text files (`.gz`)
- Block GZIP (BGZIP) for random access (`.gz` with index)
- Never compress already-compressed formats (BAM, CRAM)

**Documentation:**

- Include file format in metadata
- Specify version (VCF4.2, GFF3, etc.)
- Document any non-standard usage

---

## Genome Assembly

Genome assembly reconstructs complete genome sequences from millions of short sequencing reads or fewer long reads.

### The Assembly Problem

**Challenge:** Reassemble genome from fragmented sequencing data

**Complications:**

- Repeats longer than read length (ambiguous placement)
- Sequencing errors
- Uneven coverage
- Polyploidy and heterozygosity
- Large genome size

**Analogy:** Assembling a jigsaw puzzle where:

- Pieces are randomly sampled
- Some pieces are damaged (errors)
- Many pieces look identical (repeats)
- You don't know what the picture looks like

---

### Assembly Strategies

#### De Bruijn Graph Assembly

**Used by:** SPAdes, MEGAHIT, Velvet, IDBA, ABySS

**Algorithm:**

1. **Break reads into k-mers** (k-length substrings)

   - Read: ACGTACGT (k=4)
   - K-mers: ACGT, CGTA, GTAC, TACG, ACGT

2. **Build De Bruijn graph:**

   - Nodes = (k-1)-mers
   - Edges = k-mers
   - Edge from node A to B if they overlap by k-1

3. **Find Eulerian path** through graph
   - Path visiting every edge once
   - Represents genome sequence

**Example:**

```
K=3, reads: ACGT, CGTA

K-mers: ACG, CGT, GTA
(k-1)-mers: AC, CG, GT, TA

Graph:
AC → CG → GT → TA

Assembly: ACGTA
```

**Advantages:**

- Memory efficient
- Handles large datasets
- Works well with short reads
- Fast

**Limitations:**

- K-mer selection critical
- Errors create bubbles
- Repeats create tangles

---

#### Overlap-Layout-Consensus (OLC)

**Used by:** Canu, Flye, Miniasm, HGAP, Falcon

**Algorithm:**

1. **Overlap:** Find all pairwise overlaps between reads
2. **Layout:** Build overlap graph
   - Nodes = reads
   - Edges = overlaps
3. **Consensus:** Generate consensus sequence for each contig

**Advantages:**

- Works well with long reads
- Preserves read structure
- Handles repeats better

**Limitations:**

- Computationally expensive (O(n²) overlaps)
- Memory intensive
- Slow for many reads

**Optimizations:**

- MinHash for fast overlap detection
- Sparse overlap graphs
- Hierarchical assembly

---

### Assembly Metrics

**N50:**

- Sort contigs by length (longest first)
- N50 = length of contig where 50% of total bases are in contigs ≥ this length
- Higher = better
- **Not** average or median!

**Example:**

```
Contigs: 100 bp, 200 bp, 300 bp, 400 bp, 500 bp
Total: 1500 bp

Cumulative:
500 (500), 400 (900), 300 (1200), ...
50% of 1500 = 750

N50 = 400 bp (cumulative sum passes 750)
```

**L50:**

- Number of contigs ≥ N50
- Lower = better (fewer contigs)

**Other Metrics:**

- **N90, N75:** Similar to N50 but at different percentages
- **Total length:** Should match genome size
- **Number of contigs:** Fewer = better (more contiguous)
- **Longest contig:** Indicator of best assembly
- **GC content:** Should match expected
- **Coverage:** Sequencing depth

**BUSCO (Benchmarking Universal Single-Copy Orthologs):**

- Assess completeness
- Percent of expected single-copy genes found
- Complete, fragmented, missing
- > 95% complete = high quality

```bash
busco -i assembly.fasta -l bacteria_odb10 -o busco_results -m genome
```

**QUAST (Quality Assessment Tool):**

- Comprehensive assembly evaluation
- With or without reference genome
- Generates detailed reports

```bash
quast.py assembly.fasta -r reference.fasta -g genes.gff -o quast_output
```

---

### Short-Read Assemblers

#### SPAdes

**Description:**

- De Bruijn graph-based
- Multiple k-mer approach
- Excellent for bacterial genomes
- Also: metaSPAdes (metagenomes), rnaSPAdes (transcriptomes)

**Usage:**

```bash
spades.py -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o output_dir
```

**Features:**

- Automatic k-mer selection
- Error correction
- Paired-end support
- Good for isolates

---

#### MEGAHIT

**Description:**

- Fast, memory-efficient
- Designed for metagenomes
- Handles large datasets

**Usage:**

```bash
megahit -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o output_dir
```

**Advantages:**

- Very fast
- Low memory (<10 GB for large metagenomes)
- Good for complex communities

---

### Long-Read Assemblers

#### Flye

**Description:**

- De Bruijn-like approach for long reads
- Works with PacBio and Nanopore
- Repeat-resolution

**Usage:**

```bash
flye --nano-raw reads.fastq.gz --out-dir output_dir --threads 32
```

**Features:**

- Handles high error rates
- Resolves repeats
- Generates assembly graph (GFA format)

---

#### Canu

**Description:**

- OLC assembler
- Optimized for long reads
- Three stages: correction, trimming, assembly

**Usage:**

```bash
canu -p prefix -d output_dir genomeSize=5m -nanopore reads.fastq.gz
```

**Features:**

- High-quality assemblies
- Error correction included
- Computational intensive

---

#### Miniasm

**Description:**

- Ultra-fast assembler
- No error correction
- Quick draft assemblies

**Usage:**

```bash
minimap2 -x ava-ont reads.fastq.gz reads.fastq.gz | \
miniasm -f reads.fastq.gz - > assembly.gfa
```

**Advantages:**

- Very fast (<1 hour for bacterial genome)
- Minimal resources
- Good for quick assessment

**Limitations:**

- Requires polishing (Racon, Medaka)

---

### Hybrid Assembly

Combine short and long reads for best results.

**Advantages:**

- Long reads for structure (spans repeats)
- Short reads for accuracy

**Strategies:**

1. **Long-read assembly + short-read polishing**
2. **Integrate both in assembly process**

**Tools:**

- **Unicycler:** Automated hybrid assembly (bacteria)
- **MaSuRCA:** Mega-reads approach
- **SPAdes (--pacbio/--nanopore):** Hybrid mode

**Polishing:**

After long-read assembly, improve accuracy:

```bash
# Pilon (short reads)
bwa mem assembly.fasta reads_R1.fq reads_R2.fq | samtools sort > aligned.bam
samtools index aligned.bam
pilon --genome assembly.fasta --frags aligned.bam --output polished

# Racon (long reads)
minimap2 -ax map-ont assembly.fasta reads.fastq > aligned.sam
racon reads.fastq aligned.sam assembly.fasta > polished.fasta

# Medaka (Nanopore)
medaka_consensus -i reads.fastq -d assembly.fasta -o medaka_output
```

---

### Metagenome Assembly

Special considerations for mixed communities.

**Challenges:**

- Multiple organisms
- Varying abundances (rare species)
- Strain variation
- Higher complexity

**Assemblers:**

- **metaSPAdes:** SPAdes for metagenomes
- **MEGAHIT:** Fast, memory-efficient
- **metaFlye:** Long-read metagenome assembly
- **IDBA-UD:** Uneven depth

**Binning:**

After assembly, group contigs by organism.

**Relevant to GraphBin:**

- **GraphBin:** Uses assembly graph for improved binning
- Alternative to traditional coverage/composition binning
- Leverages graph connectivity

**Binning Tools:**

- **MetaBAT2:** Coverage and composition
- **MaxBin2:** EM algorithm
- **CONCOCT:** Composition and coverage
- **GraphBin:** Graph-based refinement (this project!)

---

### Assembly Graph Formats

**GFA (Graphical Fragment Assembly):**

- Standard format
- Stores assembly graph structure
- Used by SPAdes, Flye, Canu

**Example:**

```
S   contig1   ACGTACGT
S   contig2   TGCATGCA
L   contig1   +   contig2   +   4M
```

**GraphBin Usage:**

- Takes assembly graph (GFA)
- Initial binning results
- Refines bins using graph connectivity

---

### Best Practices

**Pre-Assembly:**

1. **Quality control:**

   - FastQC for quality assessment
   - Trim adapters (Trimmomatic, fastp)
   - Remove low-quality reads

2. **Error correction:**

   - Built into many assemblers
   - Or use dedicated tools (Lighter, BFC)

3. **Coverage estimation:**
   - Calculate sequencing depth
   - Aim for 50-100× for bacteria, higher for eukaryotes

**Assembly:**

1. **Choose appropriate assembler:**

   - Read type (short/long)
   - Genome type (isolate/metagenome)
   - Available resources

2. **Try multiple parameters:**

   - Different k-mers
   - Coverage cutoffs

3. **Compare assemblies:**
   - Use QUAST
   - Select best by N50, completeness

**Post-Assembly:**

1. **Polish:**

   - Long reads: Racon, Medaka
   - Short reads: Pilon
   - Multiple rounds

2. **Validate:**

   - BUSCO for completeness
   - Map reads back (should be >95%)
   - Check for contamination

3. **Scaffold:**
   - Use additional data (Hi-C, optical maps)
   - Tools: SSPACE, BESST

---

## Read Mapping and Variant Analysis

### Read Mapping

Aligning sequencing reads to a reference genome to identify where they originated.

#### Why Map Reads?

**Applications:**

- Variant calling (SNPs, indels)
- Gene expression quantification (RNA-seq)
- ChIP-seq peak detection
- Metagenomics classification
- Quality control

---

#### Mapping Algorithms

**Exact Matching:**

- Fast but inflexible
- No mismatches or gaps
- Rarely used for real data

**Seed-and-Extend:**

- Find exact matches (seeds)
- Extend allowing mismatches
- Used by BWA, Bowtie

**Burrows-Wheeler Transform (BWT):**

- Compress reference for fast searching
- Memory efficient
- Used by BWA, Bowtie2

---

#### Popular Mappers

**BWA (Burrows-Wheeler Aligner):**

Best for short Illumina reads.

**Algorithms:**

- **BWA-ALN:** Original, up to 100 bp
- **BWA-MEM:** Modern, >70 bp, recommended
- **BWA-SW:** Long reads (deprecated)

**Usage:**

```bash
# Index reference
bwa index reference.fasta

# Map paired-end reads
bwa mem reference.fasta reads_R1.fq reads_R2.fq > aligned.sam

# Complete pipeline
bwa mem -t 8 reference.fasta reads_R1.fq reads_R2.fq | \
samtools sort -o aligned.sorted.bam
samtools index aligned.sorted.bam
```

---

**Bowtie2:**

Fast and accurate for short reads.

```bash
# Build index
bowtie2-build reference.fasta ref_index

# Map reads
bowtie2 -x ref_index -1 reads_R1.fq -2 reads_R2.fq -S aligned.sam
```

---

**Minimap2:**

Versatile mapper for long reads and transcripts.

**Presets:**

- `-x sr`: Short reads (Illumina)
- `-x map-ont`: Oxford Nanopore
- `-x map-pb`: PacBio
- `-x splice`: RNA-seq (spliced alignment)
- `-x asm5/asm10`: Assembly-to-assembly

```bash
# Map Nanopore reads
minimap2 -ax map-ont reference.fasta reads.fastq > aligned.sam

# Map RNA-seq
minimap2 -ax splice reference.fasta rna_reads.fastq > aligned.sam
```

---

**STAR (Spliced Transcripts Alignment to a Reference):**

Specialized for RNA-seq with splicing.

```bash
# Generate genome index
STAR --runMode genomeGenerate --genomeDir genome_dir \
     --genomeFastaFiles reference.fasta --sjdbGTFfile annotation.gtf

# Map RNA-seq reads
STAR --genomeDir genome_dir --readFilesIn reads_R1.fq reads_R2.fq \
     --outFileNamePrefix sample_ --outSAMtype BAM SortedByCoordinate
```

---

#### Mapping Quality (MAPQ)

Phred-scaled probability that mapping is incorrect:

```
MAPQ = -10 × log10(P)

MAPQ 60 = 1 in 1,000,000 chance of error
MAPQ 30 = 1 in 1,000
MAPQ 20 = 1 in 100
MAPQ 0  = Multiple equal positions (unmappable)
```

**Filtering:**

```bash
# Keep only high-quality mappings (MAPQ ≥ 30)
samtools view -q 30 -b input.bam > filtered.bam
```

---

### Variant Calling

Identifying differences between sequenced sample and reference genome.

#### Variant Types

- **SNP (Single Nucleotide Polymorphism):** Single base change
- **Indel:** Insertion or deletion
- **MNP:** Multiple adjacent nucleotide changes
- **Structural Variant (SV):** Large deletions, duplications, inversions, translocations

---

#### Variant Calling Workflow

1. **Map reads** to reference
2. **Mark duplicates** (PCR/optical)
3. **Call variants** (identify differences)
4. **Filter variants** (quality control)
5. **Annotate variants** (predict effects)

---

#### Variant Callers

**GATK (Genome Analysis Toolkit):**

Gold standard for human genetics.

**HaplotypeCaller:**

```bash
# Call variants
gatk HaplotypeCaller \
    -R reference.fasta \
    -I input.bam \
    -O variants.vcf.gz

# Joint calling (multiple samples)
gatk GenotypeGVCFs \
    -R reference.fasta \
    -V sample1.g.vcf.gz \
    -V sample2.g.vcf.gz \
    -O cohort.vcf.gz
```

**Best Practices:**

- Base quality score recalibration (BQSR)
- Indel realignment (older versions)
- Variant quality score recalibration (VQSR)

---

**FreeBayes:**

Fast, haplotype-based variant caller.

```bash
freebayes -f reference.fasta input.bam > variants.vcf
```

**Advantages:**

- Simple to use
- No training data required
- Works on any organism

---

**BCFtools mpileup/call:**

Fast and simple variant calling.

```bash
bcftools mpileup -f reference.fasta input.bam | \
bcftools call -mv -Oz -o variants.vcf.gz
```

---

**DeepVariant:**

Deep learning-based variant caller.

```bash
deepvariant \
    --model_type=WGS \
    --ref=reference.fasta \
    --reads=input.bam \
    --output_vcf=variants.vcf.gz
```

**Advantages:**

- High accuracy
- Trained on diverse datasets
- Works across organisms

---

#### Variant Filtering

**Hard Filtering:**

```bash
# BCFtools
bcftools filter -i 'QUAL>30 && DP>10' variants.vcf.gz > filtered.vcf

# GATK
gatk VariantFiltration \
    -R reference.fasta \
    -V variants.vcf.gz \
    -O filtered.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
    --filter-name "basic_filter"
```

**Common Filters:**

- **Quality (QUAL):** >20 or >30
- **Depth (DP):** >10 (too low) and <3× mean (too high, repeats)
- **Mapping quality (MQ):** >40
- **Strand bias:** Variants on only one strand suspicious
- **Allele frequency:** For heterozygous, expect ~0.5

---

#### Variant Annotation

Predict effects of variants on genes and proteins.

**VEP (Variant Effect Predictor):**

```bash
vep -i variants.vcf -o annotated.vcf --cache --force_overwrite \
    --species homo_sapiens --assembly GRCh38
```

**SnpEff:**

```bash
# Download database
java -jar snpEff.jar download GRCh38.86

# Annotate
java -jar snpEff.jar GRCh38.86 variants.vcf > annotated.vcf
```

**Annotations:**

- Gene name
- Transcript ID
- Consequence (missense, nonsense, frameshift, etc.)
- Protein change (p.Ala123Thr)
- Conservation scores
- Population frequencies
- Clinical significance

**Functional Predictions:**

- **SIFT:** Predicts deleterious amino acid substitutions
- **PolyPhen-2:** Predicts impact on protein function
- **CADD:** Combined annotation-dependent depletion score
- **GERP:** Evolutionary conservation

---

### Coverage Analysis

**Calculate Coverage:**

```bash
# Per-base coverage
samtools depth input.bam > coverage.txt

# Average coverage
samtools depth input.bam | awk '{sum+=$3} END {print sum/NR}'

# BEDTools coverage
bedtools genomecov -ibam input.bam -bg > coverage.bedgraph
```

**Visualize:**

- IGV (Integrative Genomics Viewer)
- UCSC Genome Browser
- JBrowse

---

## Computational Tools and Programming

### Command-Line Environment

**Unix/Linux Basics:**

Essential for bioinformatics.

**Key Commands:**

```bash
# Navigation
cd, ls, pwd, mkdir, rm, mv, cp

# File viewing
cat, less, head, tail, grep

# File manipulation
cut, sort, uniq, wc, awk, sed

# Compression
gzip, gunzip, tar

# Pipes and redirection
command1 | command2
command > output.txt
command >> append.txt
```

**Example Workflow:**

```bash
# Count number of sequences in FASTQ
zcat reads.fastq.gz | echo $((`wc -l`/4))

# Extract sequence IDs
grep "^>" sequences.fasta | cut -d' ' -f1

# Find files
find . -name "*.fasta" -type f
```

---

### Programming Languages

#### Python

Most popular for bioinformatics.

**Key Libraries:**

- **Biopython:** Sequence analysis, file I/O
- **NumPy/Pandas:** Data manipulation
- **Matplotlib/Seaborn:** Visualization
- **SciPy/scikit-learn:** Statistical analysis, machine learning
- **NetworkX:** Graph analysis (relevant to GraphBin!)

**Example:**

```python
from Bio import SeqIO

# Parse FASTA
for record in SeqIO.parse("sequences.fasta", "fasta"):
    print(f"{record.id}: {len(record.seq)} bp")

# Calculate GC content
def gc_content(seq):
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq) * 100
```

---

#### R

Statistical analysis and visualization.

**Key Packages:**

- **Bioconductor:** Genomics packages
- **ggplot2:** Beautiful plots
- **dplyr/tidyr:** Data manipulation
- **DESeq2:** RNA-seq differential expression
- **edgeR:** Another RNA-seq tool
- **GenomicRanges:** Genomic intervals

**Example:**

```r
library(Biostrings)

# Read sequences
seqs <- readDNAStringSet("sequences.fasta")

# GC content
gc_content <- letterFrequency(seqs, "GC", as.prob=TRUE)

# Plot
hist(gc_content, main="GC Content Distribution")
```

---

#### Bash Scripting

Automate workflows.

**Example:**

```bash
#!/bin/bash

# Quality control and mapping pipeline

# Variables
REF="reference.fasta"
THREADS=8

# Loop through samples
for SAMPLE in sample1 sample2 sample3; do
    echo "Processing $SAMPLE"

    # QC
    fastqc ${SAMPLE}_R1.fq ${SAMPLE}_R2.fq

    # Map
    bwa mem -t $THREADS $REF ${SAMPLE}_R1.fq ${SAMPLE}_R2.fq | \
        samtools sort -o ${SAMPLE}.bam

    # Index
    samtools index ${SAMPLE}.bam

    # Stats
    samtools flagstat ${SAMPLE}.bam > ${SAMPLE}.stats.txt
done
```

---

### Workflow Management

**Snakemake:**

Python-based workflow manager.

**Example:**

```python
rule all:
    input: "results/final_report.txt"

rule map_reads:
    input:
        ref="reference.fasta",
        r1="reads_R1.fq",
        r2="reads_R2.fq"
    output: "mapped.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | "
        "samtools sort -o {output}"

rule call_variants:
    input: "mapped.bam"
    output: "variants.vcf"
    shell:
        "bcftools mpileup -f {input.ref} {input} | "
        "bcftools call -mv -o {output}"
```

**Other Tools:**

- **Nextflow:** Portable, scalable
- **CWL (Common Workflow Language):** Standard format
- **WDL (Workflow Description Language):** Used by GATK

---

### High-Performance Computing

**Parallel Processing:**

```bash
# GNU Parallel
parallel -j 8 "bwa mem ref.fa {}.fq | samtools sort -o {}.bam" ::: sample*

# Job arrays (SLURM)
#SBATCH --array=1-100
```

**Considerations:**

- CPU vs memory vs I/O bottlenecks
- Parallelize embarrassingly parallel tasks
- Use appropriate data structures
- Profile code to find bottlenecks

---

### Version Control

**Git:**

Essential for collaboration and reproducibility.

```bash
# Initialize repository
git init
git add script.py
git commit -m "Add analysis script"

# Push to GitHub
git remote add origin https://github.com/user/repo.git
git push -u origin main

# Collaborate
git pull
git branch new-feature
git checkout new-feature
```

---

### Best Practices

**Reproducibility:**

- Document everything
- Use version control
- Save parameters and software versions
- Use workflow managers
- Create environments (conda, Docker)

**Code Quality:**

- Comment your code
- Use meaningful variable names
- Modularize (functions, modules)
- Test your code
- Follow style guides (PEP8 for Python)

**Data Management:**

- Organize directory structure
- Use consistent naming
- Backup important data
- Document file formats
- Keep raw data read-only

---

## Summary and Learning Resources

### Key Takeaways

**Bioinformatics Foundations:**

1. **Sequence Alignment:** Fundamental operation for comparing sequences

   - Pairwise (Needleman-Wunsch, Smith-Waterman, BLAST)
   - Multiple (MUSCLE, MAFFT, Clustal)
   - Profile searches (PSI-BLAST, HMMER)

2. **Motifs and Patterns:** Identifying functional elements

   - Consensus sequences, PWMs, sequence logos
   - Discovery (MEME, DREME) and scanning (FIMO)
   - Applications in regulation and function prediction

3. **Phylogenetics:** Understanding evolutionary relationships

   - Distance methods (NJ) vs character-based (ML, Bayesian)
   - Tree evaluation (bootstrap, posterior probabilities)
   - Software: IQ-TREE, RAxML, MrBayes

4. **Databases:** Essential resources for bioinformatics

   - Primary: GenBank, ENA, DDBJ
   - Curated: RefSeq, UniProt, Ensembl
   - Functional: KEGG, GO, Pfam, InterPro
   - Specialized: GTDB (microbial taxonomy)

5. **File Formats:** Standard data representations

   - Sequences: FASTA, FASTQ
   - Alignments: SAM/BAM/CRAM
   - Variants: VCF/BCF
   - Annotations: GFF3/GTF, BED
   - Graphs: GFA (assembly graphs)

6. **Genome Assembly:** Reconstructing genomes from reads

   - De Bruijn graphs (SPAdes, MEGAHIT)
   - OLC (Canu, Flye)
   - Metrics: N50, BUSCO
   - Metagenome assembly and binning

7. **Read Mapping and Variants:** From reads to biological insights

   - Mappers: BWA, Bowtie2, Minimap2, STAR
   - Variant calling: GATK, FreeBayes, BCFtools
   - Annotation: VEP, SnpEff
   - Quality control throughout

8. **Computational Skills:** Tools for analysis
   - Command line (Unix/Linux)
   - Programming (Python, R, Bash)
   - Workflows (Snakemake, Nextflow)
   - Version control (Git)
   - HPC and parallelization

---

### Learning Pathways

**Beginner:**

1. Learn Unix/Linux basics
2. Understand molecular biology fundamentals
3. Learn one programming language (Python or R)
4. Practice with small datasets
5. Use online resources and tutorials

**Intermediate:**

1. Master alignment algorithms and tools
2. Understand statistical concepts
3. Work with real datasets
4. Learn workflow management
5. Contribute to open-source projects

**Advanced:**

1. Develop new algorithms
2. Optimize computational performance
3. Machine learning applications
4. Method development and benchmarking
5. Publish and share tools

---

### Online Resources

**Courses:**

- **Rosalind:** Bioinformatics problem-solving platform
- **Coursera:** Bioinformatics specialization (UC San Diego)
- **edX:** Multiple bioinformatics courses
- **Canadian Bioinformatics Workshops:** Comprehensive modules
- **EMBL-EBI Training:** Free courses and webinars

**Tutorials:**

- **NCBI Education:** Database tutorials
- **Galaxy Training:** Workflow-based tutorials
- **Software Carpentry:** Programming for scientists
- **Bioconductor:** R/Bioconductor workshops

**Books:**

- "Bioinformatics Data Skills" by Vince Buffalo
- "Biological Sequence Analysis" by Durbin et al.
- "Bioinformatics Algorithms" by Compeau & Pevzner
- "Statistical Methods in Bioinformatics" by Ewens & Grant

**Communities:**

- **Biostars:** Q&A forum
- **SEQanswers/SEQannotations:** Sequencing forums
- **Reddit:** r/bioinformatics
- **Twitter:** #bioinformatics community
- **GitHub:** Open-source projects

**Documentation:**

- Software manuals (GATK, Bioconductor, etc.)
- NCBI Handbook
- Ensembl documentation
- Tool-specific wikis and forums

---

### Staying Current

**Journals:**

- _Bioinformatics_
- _Nucleic Acids Research_
- _Genome Biology_
- _Nature Methods_
- _PLOS Computational Biology_
- _Briefings in Bioinformatics_

**Preprint Servers:**

- **bioRxiv:** Biology preprints
- **arXiv q-bio:** Quantitative biology

**Conferences:**

- ISMB (Intelligent Systems for Molecular Biology)
- RECOMB (Research in Computational Molecular Biology)
- ASHG (American Society of Human Genetics)
- AGBT (Advances in Genome Biology and Technology)

**Software Updates:**

- Follow tool GitHub repositories
- Subscribe to mailing lists
- Check for new versions regularly

---

### Practical Advice

**For Students:**

- Start with simple analyses
- Don't be afraid to ask questions
- Practice regularly
- Work on real projects
- Collaborate with biologists
- Present your work

**For Researchers:**

- Document everything
- Make analyses reproducible
- Share code and data
- Validate computational results
- Keep learning new methods
- Contribute to the community

**For Tool Developers:**

- Write clear documentation
- Provide examples and tutorials
- Use standard formats
- Make code open-source
- Benchmark against existing tools
- Respond to user feedback

---

### The Future of Bioinformatics

**Emerging Trends:**

- **Single-cell and spatial omics:** Higher resolution data
- **Long-read technologies:** Complete, accurate genomes
- **AI and machine learning:** Improved predictions
- **Cloud computing:** Scalable analyses
- **Multi-omics integration:** Systems-level understanding
- **Real-time analysis:** Nanopore adaptive sampling
- **Reproducible research:** Containers, workflows
- **Data sharing:** FAIR principles (Findable, Accessible, Interoperable, Reusable)

**Skills in Demand:**

- Deep learning and AI
- Graph algorithms (like GraphBin!)
- Cloud computing (AWS, Google Cloud)
- Workflow development
- Statistical genetics
- Data visualization
- Large-scale data management

---

### Conclusion

Bioinformatics is a rapidly evolving field at the intersection of biology, computer science, statistics, and mathematics. The foundations covered in this document—sequence alignment, phylogenetics, databases, file formats, assembly, variant analysis, and computational tools—provide the essential knowledge for working in genomics and bioinformatics.

Success in bioinformatics requires:

- **Solid computational skills:** Programming, command line, algorithms
- **Biological knowledge:** Understanding the questions and context
- **Statistical thinking:** Proper analysis and interpretation
- **Continuous learning:** Field evolves rapidly
- **Collaboration:** Work with diverse experts
- **Attention to detail:** Small errors can have big impacts

As sequencing technologies continue to advance and biological questions become more complex, bioinformatics will play an increasingly central role in biology and medicine. Tools like **GraphBin**, which leverage assembly graph structure for improved metagenome binning, exemplify how computational innovation drives biological discovery.

Whether you're analyzing genomes, studying gene expression, classifying microbes, or developing new computational methods, these foundational concepts and tools will serve as the building blocks for your work in bioinformatics.

---

**Document Information:**

- **Created:** 2025
- **Purpose:** Comprehensive bioinformatics foundations guide
- **Scope:** Sequence analysis, phylogenetics, databases, assembly, variant calling, computation
- **Target Audience:** Students, researchers, bioinformaticians
- **Related:** GraphBin documentation (metagenome binning using assembly graphs)

---
