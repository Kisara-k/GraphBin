# Domain Knowledge: Genomics - Deep Dive

## Table of Contents

1. [Introduction to Genomics](#introduction-to-genomics)
2. [DNA Structure and Function](#dna-structure-and-function)
3. [Genome Organization](#genome-organization)
4. [Gene Structure and Expression](#gene-structure-and-expression)
5. [Genome Evolution and Variation](#genome-evolution-and-variation)
6. [Microbial Genomics](#microbial-genomics)
7. [Comparative Genomics](#comparative-genomics)
8. [Functional Genomics](#functional-genomics)
9. [Pangenomics](#pangenomics)
10. [Genomic Technologies](#genomic-technologies)
11. [Bioinformatics Foundations](#bioinformatics-foundations)
12. [Applications in Metagenomics](#applications-in-metagenomics)

---

## Introduction to Genomics

### Defining Genomics

**Genomics** is the comprehensive study of genomes - the complete set of DNA (or RNA in some viruses) of an organism, including all of its genes and non-coding sequences.

**Etymology:**

- **Genome** = Gene + Chromosome (coined 1920)
- **Genomics** = Genome + -omics (suffix meaning "comprehensive study of")

**Scope:**

- **Structural genomics:** Organization and content of genomes
- **Functional genomics:** How genomes function
- **Comparative genomics:** Evolution and relationships
- **Metagenomics:** Genomes in communities (GraphBin's domain!)

### Historical Milestones

#### Pre-Genomic Era

**1866:** Gregor Mendel

- Laws of inheritance
- Concept of discrete "factors" (genes)

**1869:** Friedrich Miescher

- Isolated "nuclein" (DNA) from white blood cells
- Didn't know its function

**1944:** Oswald Avery, Colin MacLeod, Maclyn McCarty

- Proved DNA is genetic material
- "Transforming principle" experiment

**1953:** James Watson, Francis Crick (with Rosalind Franklin's data)

- DNA double helix structure
- Explained heredity mechanism
- Foundation for molecular biology

**1977:** Frederick Sanger

- First genome sequenced: φX174 bacteriophage (5,386 bp)
- Sanger sequencing method
- Won second Nobel Prize

#### Genomic Era Begins

**1995:** First free-living organism genome

- _Haemophilus influenzae_ (1.8 Mb)
- Craig Venter and team
- Whole-genome shotgun approach

**1996:** First eukaryote genome

- _Saccharomyces cerevisiae_ (yeast, 12 Mb)
- International collaboration
- Model organism

**1998:** First animal genome

- _Caenorhabditis elegans_ (worm, 97 Mb)
- Multicellular organism
- Complex development

**2000-2003:** Human Genome Project

- 3 billion base pairs
- ~20,000-25,000 genes (fewer than expected!)
- $3 billion, 13 years
- Public consortium + Celera (private)
- "Book of life" decoded

#### Modern Genomics Era

**2007:** Personal genome era

- James Watson's genome sequenced
- Craig Venter's genome
- Cost dropping rapidly

**2008:** 1000 Genomes Project

- Human genetic variation
- Population genomics
- Common variants catalog

**2010s:** Genome sequencing routine

- Thousands of genomes sequenced
- Clinical applications
- Agricultural genomics
- Microbiome projects

**Present:** Post-Genomic Era

- $100-1000 per genome
- Clinical genomics standard
- Personalized medicine
- CRISPR gene editing
- Synthetic biology

### Why Genomics Matters

#### Scientific Understanding

**Life's Blueprint:**

- How organisms develop
- How they function
- How they evolve
- How they interact

**Universal Principles:**

- All life uses same genetic code (with minor variations)
- Reveals common ancestry
- Evolutionary relationships
- Conservation of mechanisms

#### Medical Applications

**Disease Understanding:**

- Genetic disorders (cystic fibrosis, sickle cell)
- Cancer genomics (personalized treatment)
- Infectious diseases (pathogen identification)
- Pharmacogenomics (drug response)

**Diagnosis and Treatment:**

- Rapid pathogen identification
- Antibiotic resistance detection
- Targeted therapies
- Preventive medicine

#### Agricultural Applications

**Crop Improvement:**

- Higher yields
- Disease resistance
- Drought tolerance
- Nutritional enhancement

**Livestock:**

- Breeding selection
- Disease resistance
- Production traits

#### Industrial Applications

**Biotechnology:**

- Enzyme production
- Biofuels
- Biomaterials
- Biosensors

**Synthetic Biology:**

- Designed organisms
- Metabolic engineering
- Drug production

---

## DNA Structure and Function

### Chemical Structure

#### The Building Blocks

**Nucleotides:** Basic units of DNA

**Components:**

1. **Sugar:** Deoxyribose (5-carbon sugar)
2. **Phosphate group:** Connects nucleotides
3. **Nitrogenous base:** Carries genetic information

**Four Bases:**

**Purines (double ring):**

- **Adenine (A)**
- **Guanine (G)**

**Pyrimidines (single ring):**

- **Cytosine (C)**
- **Thymine (T)**

**In RNA:** Thymine replaced by **Uracil (U)**

#### Base Pairing Rules

**Watson-Crick Pairing:**

```
A ≡ T  (Adenine pairs with Thymine, 2 hydrogen bonds)
G ≡ C  (Guanine pairs with Cytosine, 3 hydrogen bonds)
```

**Why These Pairs?**

1. **Size complementarity:** Purine + Pyrimidine = constant width
2. **Hydrogen bonding:** Specific number and geometry
3. **Chemical stability:** GC stronger (3 bonds vs. 2)

**Wobble Pairing:**

- Less stringent rules in some contexts (tRNA)
- G can pair with U

#### The Double Helix

**Structure:**

```
    5' ────→ 3'  (One strand)

    A ≡ T
    T ≡ A
    G ≡ C
    C ≡ G
    A ≡ T

    3' ←──── 5'  (Complementary, antiparallel)
```

**Key Features:**

**Antiparallel Strands:**

- One strand: 5' → 3'
- Other strand: 3' → 5'
- Critical for replication and function

**Major and Minor Grooves:**

- Major groove: Wider, more accessible
- Minor groove: Narrower
- Proteins bind preferentially to grooves

**Right-Handed Helix:**

- B-form DNA (most common in cells)
- ~10 base pairs per turn
- ~3.4 nm pitch

**Alternative Forms:**

- **A-DNA:** Right-handed, shorter, wider (in dehydration)
- **Z-DNA:** Left-handed, zigzag (in high salt, GC-rich)

### DNA Properties

#### Melting and Annealing

**Denaturation (Melting):**

```
Double-stranded DNA ─heat→ Two single strands
(dsDNA)                    (ssDNA)
```

**Melting Temperature (Tm):**

- Temperature at which 50% DNA denatured
- **Higher GC content → Higher Tm** (3 bonds vs. 2)
- Formula: Tm ≈ 2(A+T) + 4(G+C) (simple estimate)

**Renaturation (Annealing):**

```
Two single strands ─cool→ Double-stranded DNA
(ssDNA)                   (dsDNA)
```

**Applications:**

- PCR (Polymerase Chain Reaction)
- Hybridization assays
- DNA sequencing
- Assembly algorithms use complementarity

#### Stability Factors

**GC Content:**

- **High GC (>60%):** Very stable, high Tm
- **Low GC (<30%):** Less stable, low Tm
- Affects DNA structure and function

**Length:**

- Longer sequences more stable
- More total hydrogen bonds
- Cooperative melting

**Ionic Strength:**

- Salts stabilize DNA (shield negative charges)
- Low salt → DNA unstable

**pH:**

- Extreme pH denatures DNA
- Protonation/deprotonation of bases

### DNA Replication

#### Semiconservative Replication

**Principle:** Each strand serves as template

```
Original:    ════════════

Replication: ════════════  ════════════
             (Original +   (Original +
              New strand)   New strand)
```

**Proven by:** Meselson-Stahl experiment (1958)

#### The Replication Machinery

**Origin of Replication (oriC):**

- Where replication starts
- AT-rich region (easier to open)
- Bacteria: Single origin
- Eukaryotes: Multiple origins

**DNA Helicase:**

- Unwinds double helix
- Breaks hydrogen bonds
- Creates replication fork

**DNA Polymerase:**

- Synthesizes new DNA
- 5' → 3' direction only
- Requires primer (RNA or DNA)
- Proofreading activity (3' → 5' exonuclease)

**Leading vs. Lagging Strand:**

```
       5' ←───────────── 3'  Template
       3' ─────────────→ 5'  Leading (continuous)
            Replication fork
       5' ←─────────────── 3'  Template
       3' ←──←──←──←──←── 5'  Lagging (discontinuous, Okazaki fragments)
```

**Leading Strand:**

- Synthesized continuously
- Same direction as fork movement

**Lagging Strand:**

- Synthesized discontinuously
- Okazaki fragments (1000-2000 bp bacteria, 100-200 bp eukaryotes)
- DNA ligase joins fragments

#### Replication Errors

**Error Rate:**

- DNA polymerase: 1 error per 10⁷-10⁸ bases (with proofreading)
- Mismatch repair: Reduces to 1 per 10⁹-10¹⁰ bases

**Consequences:**

- Mutations (can be beneficial, neutral, or harmful)
- Evolution substrate
- Disease (if in germ cells)
- Cancer (if in somatic cells)

### DNA Repair Mechanisms

**Types of Damage:**

1. **Spontaneous:** Depurination, deamination, oxidation
2. **Environmental:** UV radiation, chemicals, radiation
3. **Replication errors:** Mismatches

**Repair Systems:**

**Mismatch Repair:**

- Fixes replication errors
- Recognizes distortions in helix

**Base Excision Repair:**

- Removes damaged bases
- DNA glycosylases

**Nucleotide Excision Repair:**

- Removes bulky lesions (thymine dimers)
- Cuts out damaged region

**Double-Strand Break Repair:**

- Homologous recombination
- Non-homologous end joining

**Clinical Relevance:**

- Defects cause cancer predisposition
- Lynch syndrome (mismatch repair deficiency)

---

## Genome Organization

### Prokaryotic Genomes

#### General Features

**Size Range:**

- **Typical:** 0.5-10 Mb
- **Smallest:** _Mycoplasma genitalium_ (0.58 Mb)
- **Largest bacteria:** _Sorangium cellulosum_ (13 Mb)

**Structure:**

- **Circular chromosome** (usually)
- **Linear chromosomes** (some, e.g., _Borrelia_)
- **Single origin of replication**
- **Plasmids:** Extra-chromosomal DNA

**Gene Density:**

- Very high: 85-95% coding
- Little non-coding DNA
- Minimal intergenic space

**Example: _E. coli_**

```
Genome size: 4.6 Mb
Genes: ~4,300
Coding: ~87%
Average gene length: ~1 kb
Intergenic regions: ~130 bp average
```

#### Operons

**Definition:** Multiple genes transcribed as single mRNA

**Classic Example: lac Operon**

```
Promoter─Operator─lacZ─lacY─lacA
         (Control)  (β-gal)(Permease)(Transacetylase)

Single transcript → Polycistronic mRNA → Multiple proteins
```

**Advantages:**

- Coordinate regulation
- Functional gene clustering
- Efficient gene expression

**Common in Bacteria:**

- ~50% of genes in operons
- Functionally related genes grouped

#### Mobile Genetic Elements

**Plasmids:**

- Circular, extrachromosomal DNA
- Size: 1-400 kb typically
- Copy number: 1-100+ per cell
- Traits: Antibiotic resistance, virulence, metabolic capabilities

**Insertion Sequences (IS):**

- Simplest transposons
- 0.7-2.5 kb
- Encode only transposase
- Can disrupt genes

**Transposons:**

- Larger than IS elements
- Carry additional genes (e.g., antibiotic resistance)
- "Jumping genes"

**Bacteriophages:**

- Viral DNA integrated in chromosome
- Can transfer genes between bacteria
- Lysogenic vs. lytic cycle

**Integrons:**

- Capture and express gene cassettes
- Role in antibiotic resistance spread

### Eukaryotic Genomes

#### General Features

**Size Range (enormous variation):**

- **Yeast:** 12 Mb
- **Human:** 3,200 Mb (3.2 Gb)
- **Wheat:** 16,000 Mb (16 Gb)
- **Paris japonica** (plant): 150,000 Mb (150 Gb)!

**C-Value Paradox:**

- No correlation between genome size and complexity
- Much non-coding DNA in large genomes

**Structure:**

- **Multiple linear chromosomes**
- **Centromeres:** Chromosome attachment during division
- **Telomeres:** Chromosome ends, protective caps
- **Multiple origins of replication**

**Gene Density:**

- Lower than prokaryotes: 1-50% coding
- Humans: ~1.5% coding
- Introns, intergenic regions

#### Chromatin Organization

**Nucleosome:**

- DNA wrapped around histone octamer
- ~147 bp per nucleosome
- Basic unit of chromatin

**Structure:**

```
DNA ─┐
     ├─ Wrap around histones
     │  (H2A, H2B, H3, H4) × 2
     └─ ~1.65 turns

Linker DNA (~20-80 bp) → Next nucleosome
```

**Higher-Order Structure:**

1. **Beads-on-a-string:** 11 nm fiber
2. **30 nm fiber:** Solenoid structure
3. **Loops and domains:** 300 nm
4. **Condensed chromosome:** Mitotic/meiotic

**Chromatin States:**

- **Euchromatin:** Open, transcriptionally active
- **Heterochromatin:** Condensed, transcriptionally silent
  - Constitutive: Always condensed (centromeres)
  - Facultative: Conditionally condensed (X-inactivation)

**Epigenetic Modifications:**

- **DNA methylation:** CpG islands, gene silencing
- **Histone modifications:** Acetylation, methylation, phosphorylation
- Heritable without DNA sequence change

#### Chromosome Structure

**Centromere:**

- Attachment point for spindle fibers
- AT-rich regions (humans)
- Satellite DNA repeats
- Essential for segregation

**Telomere:**

- Repetitive sequences at chromosome ends
- Humans: TTAGGG repeats
- Protects from degradation
- Shortens with each division (except in stem cells, cancer)
- Telomerase replenishes

**Chromosome Banding:**

- G-banding (Giemsa staining)
- AT-rich (A-bands) vs. GC-rich (G-bands)
- Karyotyping for abnormalities

### Genome Architecture

#### Gene Distribution

**Gene-Rich Regions:**

- GC-rich
- Open chromatin
- High recombination

**Gene-Poor Regions:**

- AT-rich
- Heterochromatin
- Low recombination

**Isochores:**

- Large regions (>300 kb) of homogeneous GC content
- Correlate with gene density

#### Repetitive DNA

**Types:**

**1. Tandem Repeats:**

**Satellite DNA:**

- Very long arrays (100s kb to Mb)
- Centromeres, heterochromatin

**Minisatellites:**

- 10-100 bp repeat unit
- Arrays 0.5-30 kb
- VNTR (Variable Number Tandem Repeats)
- Forensics (DNA fingerprinting)

**Microsatellites (STR):**

- 1-6 bp repeat unit
- Arrays <150 bp
- Highly polymorphic
- Genetic markers

**Example:**

```
CA repeat: CACACACACACACACA
           (CA)₈
```

**2. Interspersed Repeats:**

**SINEs (Short Interspersed Nuclear Elements):**

- <500 bp
- RNA polymerase III transcribed
- Humans: Alu elements (~1.1 million copies, ~11% of genome)

**LINEs (Long Interspersed Nuclear Elements):**

- > 6 kb
- RNA polymerase II transcribed
- Humans: LINE-1 (~500,000 copies, ~17% of genome)
- Autonomous (encode own transposase)

**DNA Transposons:**

- Move via DNA intermediate
- Mostly inactive in mammals
- Active in some organisms

**Retrotransposons:**

- Move via RNA intermediate
- LINEs and SINEs
- LTR retrotransposons

**Impact:**

- ~45% of human genome from transposable elements
- Mostly inactive ("fossil" transposons)
- Ongoing transposition causes mutations
- Source of genetic diversity

---

## Gene Structure and Expression

### Gene Structure

#### Prokaryotic Genes

**Simple Structure:**

```
Promoter ─ RBS ─ Start (ATG) ─ Coding Sequence ─ Stop ─ Terminator
         (Shine-  (Initiation)  (Open Reading    (TAA/  (Rho-dependent/
          Dalgarno)              Frame, ORF)      TAG/   independent)
                                                  TGA)
```

**Characteristics:**

- No introns (usually)
- Polycistronic (operons)
- Shine-Dalgarno sequence (ribosome binding)
- Translation can begin during transcription

**Promoter Elements:**

- **-10 box:** TATAAT (Pribnow box)
- **-35 box:** TTGACA
- Recognized by sigma factor

#### Eukaryotic Genes

**Complex Structure:**

```
Promoter ─ 5'UTR ─ Exon 1 ─ Intron 1 ─ Exon 2 ─ Intron 2 ─ Exon 3 ─ 3'UTR ─ PolyA signal
(TATA box,                 (Spliced out)       (Spliced out)       (Untranslated)
 CAAT box,
 GC box)
```

**Characteristics:**

- Introns and exons
- Monocistronic (one gene per mRNA)
- 5' cap and 3' poly(A) tail
- RNA processing required

**Exons:**

- Expressed sequences
- Coding and UTRs
- Alternative splicing

**Introns:**

- Intervening sequences
- Spliced out
- Can be huge (>100 kb)
- Contain regulatory elements

**Example: Human Dystrophin Gene**

- Gene: 2.4 Mb (largest human gene)
- mRNA: 14 kb
- Protein: 3685 amino acids
- 99.4% of gene is introns!

### Transcription

#### Overview

**Central Dogma:**

```
DNA ─Transcription→ RNA ─Translation→ Protein
```

**Process:**

1. **Initiation:** RNA polymerase binds promoter
2. **Elongation:** RNA synthesized 5' → 3'
3. **Termination:** Release of RNA and polymerase

#### Prokaryotic Transcription

**RNA Polymerase:**

- Single enzyme for all genes
- Core enzyme + sigma factor
- Sigma factor recognizes promoter

**Regulation:**

- Operons (coordinate control)
- Repressors and activators
- Attenuation
- Riboswitches

#### Eukaryotic Transcription

**Three RNA Polymerases:**

- **Pol I:** rRNA genes (except 5S)
- **Pol II:** mRNA, most ncRNA
- **Pol III:** tRNA, 5S rRNA, U6 snRNA

**Promoter Elements (Pol II):**

- **Core promoter:** TATA box, Initiator, DPE
- **Proximal elements:** CAAT box, GC box
- **Enhancers:** Can be far away (100s kb), orientation-independent

**Transcription Factors:**

- General TFs (TFIIA, TFIIB, etc.)
- Specific TFs (activate or repress)
- Mediator complex

**RNA Processing:**

**5' Capping:**

- 7-methylguanosine cap
- Protection from degradation
- Translation initiation

**Splicing:**

- Remove introns
- Join exons
- Spliceosome (snRNPs)
- Alternative splicing → protein diversity

**3' Polyadenylation:**

- Add poly(A) tail (~200 As)
- Stability
- Translation
- mRNA export

### Translation

#### Genetic Code

**Codon Table:**

- 64 codons (4³ combinations)
- 61 sense codons (code for amino acids)
- 3 stop codons (UAA, UAG, UGA)
- Start codon: AUG (methionine)

**Properties:**

- **Universal:** Nearly identical in all life
- **Degenerate:** Multiple codons per amino acid (usually)
- **Non-overlapping:** Read in triplets
- **No punctuation:** Continuous reading frame

**Wobble Position:**

- Third position of codon
- Allows flexibility (wobble base pairing)
- Reduces effect of mutations

#### The Translation Process

**Initiation:**

```
Ribosome (small subunit) → mRNA
↓
tRNA (with Met) → Start codon (AUG)
↓
Large subunit joins → Translation begins
```

**Elongation:**

1. **Aminoacyl-tRNA binding:** Correct tRNA enters A site
2. **Peptide bond formation:** Transfer amino acid to growing chain
3. **Translocation:** Ribosome moves 3 nucleotides

**Ribosome Sites:**

- **A site:** Aminoacyl-tRNA (incoming)
- **P site:** Peptidyl-tRNA (growing chain)
- **E site:** Exit (empty tRNA leaves)

**Termination:**

- Stop codon (UAA, UAG, UGA) enters A site
- Release factor binds
- Peptide released, ribosome dissociates

**Post-Translational Modifications:**

- Phosphorylation
- Glycosylation
- Ubiquitination
- Proteolytic cleavage
- Folding (chaperones)

### Gene Regulation

#### Levels of Regulation

**1. Transcriptional:**

- Promoter accessibility (chromatin)
- Transcription factor binding
- Enhancers and silencers

**2. Post-Transcriptional:**

- RNA splicing (alternative)
- RNA stability
- RNA localization
- RNA editing

**3. Translational:**

- Ribosome binding
- Upstream ORFs
- RNA secondary structure
- miRNA regulation

**4. Post-Translational:**

- Protein modifications
- Protein stability
- Protein localization
- Protein-protein interactions

#### Regulatory Mechanisms

**Prokaryotic:**

- Negative control: Repressors (lac repressor)
- Positive control: Activators (CAP-cAMP)
- Attenuation: Premature termination (trp operon)
- Two-component systems: Sensor-regulator

**Eukaryotic:**

- Chromatin remodeling
- Transcription factors (combinatorial)
- Enhanceosomes
- miRNAs (microRNAs)
- lncRNAs (long non-coding RNAs)
- RNA interference (RNAi)

---

## Genome Evolution and Variation

### Sources of Variation

#### Mutations

**Point Mutations:**

**1. Substitutions:**

- **Transition:** Purine ↔ Purine (A ↔ G) or Pyrimidine ↔ Pyrimidine (C ↔ T)
- **Transversion:** Purine ↔ Pyrimidine (A/G ↔ C/T)

**Effects:**

- **Silent (Synonymous):** No amino acid change (wobble position)
- **Missense:** Different amino acid (can be conservative or non-conservative)
- **Nonsense:** Create stop codon (truncated protein)

**2. Insertions/Deletions (Indels):**

- Add or remove nucleotides
- **Frameshift:** If not multiple of 3, changes reading frame
- Can be devastating or minimal

**Example:**

```
Original:  ATG CAT GAT TAG
           Met His Asp Stop

+1 Indel:  ATG GCA TGA TTA G
           Met Ala Stop (truncated!)
```

**Rates:**

- Bacteria: ~10⁻⁹ to 10⁻¹⁰ per bp per generation
- Eukaryotes: ~10⁻⁸ to 10⁻⁹ per bp per generation
- RNA viruses: Much higher (~10⁻⁴ to 10⁻⁶)

#### Recombination

**Homologous Recombination:**

- Exchange between similar sequences
- Meiotic recombination (crossover)
- Repair mechanism
- Generates diversity

**Non-Homologous Recombination:**

- Illegitimate recombination
- Transposable elements
- Can cause rearrangements

#### Large-Scale Changes

**Copy Number Variation (CNV):**

- Deletions
- Duplications
- Can be large (kb to Mb)

**Chromosomal Rearrangements:**

- **Deletion:** Loss of segment
- **Duplication:** Extra copy of segment
- **Inversion:** Segment flipped
- **Translocation:** Transfer between chromosomes

**Whole Genome Duplication:**

- Polyploidy
- Common in plants
- Rare in animals (frogs, fish)
- Ancient duplications in vertebrates

#### Horizontal Gene Transfer (HGT)

**Mechanisms in Bacteria:**

**1. Transformation:**

- Uptake of naked DNA from environment
- Natural competence (some bacteria)
- Laboratory transformation

**2. Transduction:**

- Bacteriophage-mediated transfer
- Generalized: Random DNA
- Specialized: Specific genes near prophage

**3. Conjugation:**

- Direct transfer through pilus
- Plasmid transfer
- Chromosomal transfer (Hfr)

**Impact:**

- Major force in bacterial evolution
- Antibiotic resistance spread
- Metabolic innovations
- Adaptation to new environments

**Detection in Genomes:**

- Atypical GC content
- Atypical codon usage
- Phylogenetic incongruence
- Presence in related species

### Molecular Clock

**Principle:** Mutations accumulate over time at relatively constant rate

**Calculation:**

```
Divergence time = Genetic distance / (2 × mutation rate)

Factor of 2: Both lineages accumulating mutations
```

**Applications:**

- Date evolutionary events
- Construct phylogenetic trees
- Estimate population divergence

**Complications:**

- Rates vary across lineages
- Selection affects accumulation
- Generation time effects

### Selection and Adaptation

#### Types of Selection

**Positive (Darwinian) Selection:**

- Favors beneficial mutations
- Drives adaptation
- Signature: dN/dS > 1 (more non-synonymous than synonymous)

**Negative (Purifying) Selection:**

- Removes deleterious mutations
- Maintains function
- Signature: dN/dS < 1 (fewer non-synonymous)

**Neutral Selection:**

- No fitness effect
- Genetic drift
- Signature: dN/dS ≈ 1

**Balancing Selection:**

- Maintains polymorphism
- Heterozygote advantage
- Frequency-dependent selection
- Example: HLA genes, sickle cell

#### Adaptation

**Environmental:**

- Temperature (thermophiles, psychrophiles)
- pH (acidophiles, alkaliphiles)
- Salinity (halophiles)
- Oxygen (aerobes, anaerobes)

**Pathogenicity:**

- Host-pathogen arms race
- Immune evasion
- Virulence factors

**Metabolic:**

- New substrates
- Antibiotic resistance
- Xenobiotic degradation

---

## Microbial Genomics

### Bacterial Genomes

#### Size and Gene Content

**Trends:**

**Obligate Intracellular Parasites:**

- Small genomes (0.16-2 Mb)
- Gene loss (rely on host)
- Examples: _Mycoplasma_, _Chlamydia_, _Rickettsia_

**Free-Living:**

- Larger genomes (2-10 Mb typical)
- Metabolic versatility
- Environmental adaptation

**Soil Bacteria:**

- Often large (8-13 Mb)
- Many secondary metabolite genes
- Complex life cycles

**Examples:**

```
Mycoplasma genitalium: 0.58 Mb, 470 genes (minimal)
Escherichia coli K-12:  4.6 Mb, 4,300 genes (model)
Sorangium cellulosum:   13 Mb, 9,400 genes (giant)
```

#### Core vs. Accessory Genes

**Core Genome:**

- Essential housekeeping genes
- Present in all strains
- Conserved functions
- Examples: DNA replication, transcription, translation

**Accessory Genome:**

- Strain-specific genes
- Adaptations to niches
- Often acquired by HGT
- Examples: Antibiotic resistance, virulence factors

**Pangenome = Core + Accessory**

#### Genomic Islands

**Definition:** Large (10-200 kb) regions of probable foreign origin

**Characteristics:**

- Different GC content
- Flanked by repeats or tRNAs
- Often mobile
- Carry important functions

**Types:**

**Pathogenicity Islands:**

- Virulence factors
- Type III secretion systems
- Toxins
- Example: LEE locus in _E. coli_ O157:H7

**Metabolic Islands:**

- Degradation pathways
- Nutrient utilization
- Example: Nitrogen fixation

**Resistance Islands:**

- Multiple antibiotic resistance genes
- Integrons
- Clinical concern

**Symbiosis Islands:**

- Enable mutualistic relationships
- Example: _Mesorhizobium_ symbiosis island

### Archaeal Genomes

#### Unique Features

**Size:** 0.5-5.7 Mb typically

**Gene Density:** Similar to bacteria (85-95%)

**Characteristics:**

- Circular chromosomes (mostly)
- Operons (like bacteria)
- Transcription/translation machinery (like eukaryotes)
- Unique metabolism genes

**Extremophiles:**

- Thermophiles: Heat-stable proteins
- Halophiles: Salt adaptation
- Methanogens: Unique methanogenesis pathway
- Acidophiles: pH tolerance

**Example: _Methanocaldococcus jannaschii_**

- First archaeon sequenced (1996)
- 1.66 Mb
- Hyperthermophile (85°C)
- Methanogen
- Many unique genes (44% no homologs)

### Viral Genomes

#### Diversity

**Size Range:**

- **Smallest:** Circoviruses (1.7 kb ssDNA)
- **Largest:** Pandoravirus (2.5 Mb dsDNA, larger than some bacteria!)

**Genetic Material:**

- dsDNA (most)
- ssDNA (parvoviruses, geminiviruses)
- dsRNA (reoviruses)
- (+)ssRNA (picornaviruses, coronaviruses)
- (-)ssRNA (influenza, Ebola)
- Retroviruses (RNA → DNA)

#### Bacteriophages

**Importance:**

- Most abundant biological entities (~10³¹ on Earth)
- Major drivers of bacterial evolution
- Gene transfer agents
- Bacterial population control

**Genome Types:**

- Linear or circular
- dsDNA (most)
- ssDNA (φX174, M13)
- RNA (rare)

**Life Cycles:**

**Lytic:**

- Infect → Replicate → Lyse host
- Produce many progeny

**Lysogenic:**

- Integrate into chromosome (prophage)
- Replicate with host
- Can induce to lytic

**Impact on Host:**

- Lysogenic conversion (new traits)
- Immunity to superinfection
- Fitness costs

### Mitochondrial and Chloroplast Genomes

#### Mitochondrial Genomes

**Origin:** Endosymbiotic α-proteobacterium

**Size:**

- Animals: 15-20 kb (reduced)
- Plants: 200-2,400 kb (expanded, lots of repeats)
- Fungi: Variable

**Human mtDNA:**

- 16.6 kb
- 37 genes:
  - 13 protein-coding (respiratory chain)
  - 22 tRNAs
  - 2 rRNAs
- Circular
- Maternal inheritance
- High copy number per cell

**Features:**

- Own genetic code variations
- No introns (animals)
- High mutation rate
- Used for phylogenetics (maternal lineage)

#### Chloroplast Genomes

**Origin:** Endosymbiotic cyanobacterium

**Size:** 120-200 kb typically

**Structure:**

- Circular
- Quadripartite: LSC (large single-copy), SSC (small single-copy), two IRs (inverted repeats)

**Gene Content:**

- ~100-130 genes
- Photosynthesis genes
- rRNAs, tRNAs
- Some plastid-specific genes

**Inheritance:** Usually maternal

---

## Comparative Genomics

### Principles and Methods

#### Sequence Alignment

**Purpose:** Identify homologous regions

**Types:**

**Global Alignment:**

- Align entire sequences
- Needleman-Wunsch algorithm
- For similar-length sequences

**Local Alignment:**

- Find best local matches
- Smith-Waterman algorithm
- BLAST (heuristic, faster)
- For finding conserved domains

**Multiple Alignment:**

- Three or more sequences
- ClustalW, MUSCLE, MAFFT
- For phylogenetics, conservation

#### Homology Types

**Orthology:**

- Genes in different species from common ancestor
- Usually similar function
- Speciation event separated them
- Example: Human vs. mouse hemoglobin

**Paralogy:**

- Genes in same genome from duplication
- May have diverged in function
- Example: α-globin vs. β-globin in humans

**Xenology:**

- Genes related by horizontal transfer
- Common in bacteria
- Complicates phylogenetics

#### Synteny

**Definition:** Conserved gene order on chromosomes

**Types:**

**Conserved Synteny:**

- Same order and orientation
- Indicates recent divergence
- Functional constraints

**Shuffled Synteny:**

- Same genes, different order
- Rearrangements occurred
- More distant relationships

**Uses:**

- Predict gene function
- Identify orthologs
- Understand chromosomal evolution

### Evolutionary Analysis

#### Phylogenetics

**Gene Trees:**

- Based on single gene
- Can differ from species tree (HGT, paralogy)

**Species Trees:**

- Based on many genes (concatenated or consensus)
- Represent organismal relationships

**Methods:**

**Distance-Based:**

- UPGMA, Neighbor-Joining
- Fast, less accurate

**Maximum Parsimony:**

- Fewest evolutionary changes
- Can be misleading (long-branch attraction)

**Maximum Likelihood:**

- Statistical model of evolution
- Computationally intensive
- Accurate

**Bayesian:**

- Posterior probability
- Incorporates prior information
- Very computationally intensive

#### Rates of Evolution

**dN/dS Ratio:**

```
dN = Non-synonymous substitution rate
dS = Synonymous substitution rate

dN/dS < 1: Purifying selection
dN/dS = 1: Neutral evolution
dN/dS > 1: Positive selection
```

**Applications:**

- Identify selected genes
- Detect functional constraints
- Study adaptation

#### Genome Rearrangements

**Events:**

- Inversions
- Translocations
- Fissions/Fusions
- Duplications/Deletions

**Measurement:**

- Breakpoint distance
- Reversal distance
- DCJ (Double-Cut-and-Join) distance

**Implications:**

- Chromosome evolution
- Speciation mechanisms
- Karyotype differences

### Functional Inference

#### Guilt by Association

**Principle:** Genes that function together tend to:

- Cluster in genome (bacteria)
- Co-express
- Co-evolve
- Have similar phylogenetic profiles

**Applications:**

- Predict function of unknown genes
- Identify pathway members
- Discover protein interactions

#### Comparative Annotation Transfer

**Principle:** Transfer annotations from well-studied organisms

**Steps:**

1. Identify orthologs
2. Transfer function annotation
3. Validate with additional evidence

**Caution:**

- Functions can diverge
- Experimental validation best
- Computational prediction starting point

---

## Functional Genomics

### Gene Expression Analysis

#### Transcriptomics

**RNA-Seq:**

- Sequence all RNA in sample
- Quantify expression levels
- Discover novel transcripts
- Alternative splicing analysis

**Advantages over Microarrays:**

- No prior knowledge needed
- Single-nucleotide resolution
- Large dynamic range
- Detect isoforms

**Applications:**

- Differential expression
- Condition responses
- Time-series
- Single-cell transcriptomics

#### Proteomics

**Mass Spectrometry:**

- Identify proteins present
- Quantify abundance
- Post-translational modifications
- Protein interactions

**2D Gel Electrophoresis:**

- Separate by charge and size
- Visualize protein spot patterns
- Older method, less common now

**Proteome Complexity:**

- Alternative splicing
- Post-translational modifications
- Protein degradation
- Many proteins per gene possible

### Metabolomics

**Definition:** Study of all small molecules (metabolites) in system

**Methods:**

- Mass spectrometry
- NMR spectroscopy

**Metabolites:**

- Intermediates
- Signaling molecules
- Waste products
- Secondary metabolites

**Integration:**

```
Genome → Transcriptome → Proteome → Metabolome → Phenotype
(Potential)  (Expression)  (Function)  (Activity)   (Outcome)
```

### Functional Genomics Approaches

#### Gene Knockout/Knockdown

**Purpose:** Determine gene function by loss-of-function

**Methods:**

**Knockout (Gene Deletion):**

- Homologous recombination
- CRISPR-Cas9 (modern, efficient)
- Complete loss of function

**Knockdown (Reduced Expression):**

- RNA interference (RNAi)
- Antisense oligonucleotides
- Partial loss of function

**Phenotype Analysis:**

- What happens without the gene?
- Infer normal function
- Genetic interactions

#### Overexpression

**Purpose:** Study gain-of-function

**Methods:**

- Strong promoters
- Multiple copies
- Inducible systems

**Applications:**

- Protein production
- Study dominant effects
- Rescue experiments

#### Genome-Wide Association Studies (GWAS)

**Principle:** Link genetic variants to phenotypes

**Approach:**

```
Many individuals → Genotype + Phenotype data
                 ↓
        Statistical association
                 ↓
   Variants associated with traits
```

**Applications:**

- Disease susceptibility
- Quantitative traits
- Drug response
- Agricultural traits

**Limitations:**

- Association ≠ causation
- Complex traits multigenic
- Missing heritability
- Population structure confounds

---

## Pangenomics

### Concept and Definitions

**Pangenome:** All genes found in a species or clade

**Components:**

**Core Genome:**

- Present in all strains
- Essential functions
- Usually 70-90% of average genome (bacteria)

**Accessory (Dispensable) Genome:**

- Present in some strains
- Adaptation to niches
- Can be larger than core

**Unique Genes:**

- Strain-specific
- Recent acquisitions or losses

### Pangenome Types

**Closed Pangenome:**

```
Adding more genomes adds few new genes
Curve plateaus

Example: Bacillus anthracis (monomorphic pathogen)
```

**Open Pangenome:**

```
Adding more genomes keeps adding new genes
Curve increases

Example: Streptococcus pneumoniae (highly recombinogenic)
```

**Mathematical Model:**

```
Pangenome = Core + Accessory

Heaps' Law: P(n) = k × n^γ

Where:
P(n) = pangenome size with n genomes
k = constant
γ < 1: closed
γ ≈ 1: open
```

### Applications

#### Clinical Microbiology

**Pathogen Variability:**

- Not all _E. coli_ are pathogens
- Identify virulence genes
- Track outbreaks (strain-specific genes)

**Vaccine Design:**

- Target core antigens (all strains)
- Or variable antigens (population coverage)

**Diagnostics:**

- Species identification (core genes)
- Strain typing (variable genes)

#### Agriculture

**Crop Pathogens:**

- Understand host range
- Identify effectors
- Resistance genes

**Beneficial Microbes:**

- Nitrogen fixation genes
- Plant growth promotion
- Biocontrol agents

### Pangenome Analysis Tools

**Roary:**

- Rapid pangenome analysis
- Uses pre-annotated genomes
- Clustering-based

**Panaroo:**

- Improved accuracy
- Handles fragmented assemblies
- Graph-based approach

**PGAP:**

- Pangenome Analysis Pipeline
- Multiple algorithms

**GET_HOMOLOGUES:**

- Flexible clustering
- Multiple methods

---

## Genomic Technologies

### Genome Sequencing

**Covered in detail in other documents, but key points:**

**Generations:**

1. **First:** Sanger sequencing (still used for validation)
2. **Second:** Illumina, Ion Torrent (short reads, high accuracy)
3. **Third:** PacBio, Nanopore (long reads, lower accuracy but improving)

**Applications:**

- Whole genome sequencing (WGS)
- Targeted sequencing
- Exome sequencing
- RNA-seq (transcriptomes)

### Genome Annotation

#### Structural Annotation

**Gene Prediction:**

**Ab Initio:**

- Based on sequence signals alone
- ORFs, start/stop codons
- Promoters, splice sites
- Tools: GeneMark, Glimmer, Augustus

**Evidence-Based:**

- Use RNA-seq data
- Protein homology
- EST (Expressed Sequence Tag) mapping
- More accurate

**Combined:**

- Integrate multiple sources
- Tools: MAKER, BRAKER

**Non-Coding RNAs:**

- tRNAscan-SE (tRNAs)
- RNAmmer (rRNAs)
- Rfam (other ncRNAs)

#### Functional Annotation

**Homology-Based:**

- BLAST against databases
- Transfer annotations

**Domain-Based:**

- Pfam (protein families)
- InterPro (integrated domains)
- Functional domains indicate function

**Pathway-Based:**

- KEGG (Kyoto Encyclopedia of Genes and Genomes)
- MetaCyc
- Complete pathways

**Ontology-Based:**

- Gene Ontology (GO)
  - Biological Process
  - Molecular Function
  - Cellular Component

**Example Annotation:**

```
Gene: yhjX
Coordinates: 3,524,100..3,525,200
Strand: +
Product: Pyruvate formate-lyase
EC: 2.3.1.54
KEGG: Mixed acid fermentation
GO: GO:0006006 (glucose metabolic process)
Pfam: PF00456 (Formate-C-acetyltransferase)
```

### Genome Editing

#### CRISPR-Cas9

**Mechanism:**

1. Design guide RNA (gRNA) to target sequence
2. Cas9 nuclease cuts DNA at target
3. Cell repairs break:
   - NHEJ (Non-Homologous End Joining): Often introduces frameshift
   - HDR (Homology-Directed Repair): Can insert specific sequence

**Applications:**

- Gene knockout
- Gene editing (change sequences)
- Gene regulation (dCas9 - dead Cas9)
- Base editing (change single nucleotides)
- Prime editing (insert, delete, replace)

**Advantages:**

- Precise
- Efficient
- Versatile
- Relatively easy

**Limitations:**

- Off-target effects
- Delivery challenges
- Mosaicism
- Ethical concerns (germline editing)

#### Other Methods

**TALENs:**

- Transcription Activator-Like Effector Nucleases
- Protein-based targeting
- More specific, more complex

**Zinc Finger Nucleases:**

- Earlier method
- Difficult to design

**Homologous Recombination:**

- Traditional method
- Lower efficiency
- Still used in some systems

### Synthetic Biology

**Definition:** Design and construction of new biological parts, devices, systems

**Applications:**

**Metabolic Engineering:**

- Biofuel production
- Drug synthesis
- Biomaterials

**Minimal Genomes:**

- JCVI-syn3.0: 531 kb, 473 genes
- Understand essential functions

**Genetic Circuits:**

- Toggle switches
- Oscillators
- Sensors

**Xenobiology:**

- Expanded genetic code
- Synthetic base pairs
- Orthogonal systems

---

## Bioinformatics Foundations

### Sequence Analysis

#### BLAST (Basic Local Alignment Search Tool)

**Algorithm:**

1. **Seed:** Find short exact matches (words)
2. **Extend:** Extend matches in both directions
3. **Score:** Evaluate significance

**Variants:**

- **blastn:** Nucleotide vs. nucleotide
- **blastp:** Protein vs. protein
- **blastx:** Translated nucleotide vs. protein
- **tblastn:** Protein vs. translated nucleotide
- **tblastx:** Translated nucleotide vs. translated nucleotide

**E-value:**

- Expect value: Number of hits expected by chance
- Lower = more significant
- E < 10⁻⁵ typically significant

**Applications:**

- Find homologs
- Annotate sequences
- Identify organisms

#### Sequence Motifs and Patterns

**Consensus Sequences:**

- Conserved patterns
- Regulatory elements
- Binding sites

**Representation:**

- IUPAC codes (R = A or G, Y = C or T, etc.)
- Position Weight Matrices (PWM)
- Logos (visual)

**Tools:**

- MEME (Motif discovery)
- FIMO (Motif scanning)
- TRANSFAC (Transcription factors)

### Phylogenetic Analysis

**Distance Methods:**

- Calculate pairwise distances
- Build tree (UPGMA, NJ)
- Fast

**Character-Based:**

- Maximum Parsimony
- Maximum Likelihood
- Bayesian
- More accurate

**Software:**

- MEGA: User-friendly
- RAxML: Maximum Likelihood
- MrBayes: Bayesian
- IQ-TREE: Fast ML

**Tree Evaluation:**

- Bootstrap (confidence)
- Branch support
- Model selection

### Genomic Databases

#### Sequence Databases

**GenBank (NCBI):**

- Comprehensive
- Public submissions
- Annotated sequences

**RefSeq:**

- Curated reference sequences
- Non-redundant
- High quality

**Ensembl:**

- Eukaryotic genomes
- Comprehensive annotation
- Comparative genomics

**UniProt:**

- Protein sequences
- Functional annotation
- Swiss-Prot (curated) + TrEMBL (automatic)

#### Specialized Databases

**KEGG:**

- Pathways
- Enzymes
- Reactions

**Pfam:**

- Protein families
- Domains
- Hidden Markov Models

**GO (Gene Ontology):**

- Standardized functional terms
- Three ontologies
- Cross-species

**GTDB:**

- Genome Taxonomy Database
- Microbial genomes
- Phylogeny-based

---

## Applications in Metagenomics

### From Genomics to Metagenomics

**Single Genome:**

- One organism
- Reference available (or assembled)
- Clear gene assignments
- Direct functional inference

**Metagenome:**

- Mixed community
- Multiple genomes simultaneously
- Fragmented assemblies
- Indirect functional inference

**Challenges:**

- Separate organisms (binning!)
- Assign functions to bins
- Understand interactions
- Incomplete genomes

### Genomic Context in Binning

#### Compositional Features

**TNF Reflects:**

- Codon usage
- DNA structural properties
- Replication/repair mechanisms
- All genomic properties

**Why It Works:**

- Organisms have distinct "genomic signatures"
- Reflects evolutionary history
- Relatively stable across genome

**Limitations:**

- Similar organisms have similar TNF
- HGT regions differ
- Short sequences unreliable

#### Functional Gene Markers

**Single-Copy Core Genes:**

- Essential housekeeping genes
- Present once per genome
- Used by CheckM for completeness/contamination

**Why Useful:**

1. **Binning:** Should all be in same bin
2. **Linking:** If on different contigs, contigs likely same organism
3. **Coverage:** Should have similar coverage
4. **Validation:** Check bin quality

**Examples:**

- RecA (recombinase)
- DnaK (chaperone)
- RpoB (RNA polymerase)
- Ribosomal proteins

### GraphBin's Genomic Foundation

#### Why Connectivity Matters

**Genomic Reality:**

```
Real genome: ─A─B─C─D─E─F─G─ (continuous DNA molecule)
                ↓
Assembly:    ─A─B─  ─C─D─  ─E─F─G─ (contigs)
             Contig1  Contig2  Contig3
```

**Assembly Graph:**

```
Contig1 ──edge→ Contig2 ──edge→ Contig3
(Edges from assembly process = biological connection!)
```

**Binning Logic:**

```
If Contig1 is from Organism X
AND Contig1 connects to Contig2
THEN Contig2 likely also from Organism X
```

#### Complementary to Traditional Features

**Traditional Binning:**

- Composition (TNF): Averaged over contig
- Coverage: Single value per contig
- No spatial information

**GraphBin Adds:**

- Local connectivity
- Biological adjacency
- Resolves ambiguous cases

**Synergy:**

```
TNF says: "Could be Organism A or B"
Coverage says: "Could be Organism A or B"
Graph says: "Connected to confirmed Organism A contigs"
→ High confidence: Organism A
```

#### Handling Genomic Complexity

**Repeats:**

- Create branches in graph
- GraphBin careful about ambiguous connections
- Removes ambiguous labels

**HGT:**

- Foreign genes have different TNF
- But still connected in genome
- Graph helps keep together

**Strains:**

- Similar genomes
- Similar TNF
- Graph can help distinguish (different connection patterns)

### Integration with Genomic Knowledge

**Gene Context:**

- Genes in operons functionally related
- Conservation of gene order
- GraphBin preserves genomic context

**Pathway Completeness:**

- After binning, can check pathways
- Complete pathways → functional organism
- Partial pathways → incomplete bin or actual

**Phylogenetic Validation:**

- After binning, classify bins
- Check for consistency
- Hybrid bins phylogenetically mixed

---

## Summary

Genomics provides the foundational knowledge for understanding GraphBin's approach to metagenomic binning:

**Key Concepts:**

1. **DNA Structure and Function:**

   - Double helix, base pairing
   - Replication, transcription, translation
   - Genetic code universal
   - Foundation for all genomic analysis

2. **Genome Organization:**

   - Prokaryotes: Compact, circular, operons
   - Eukaryotes: Complex, linear, introns
   - Mobile elements common
   - Understanding structure aids analysis

3. **Gene Structure:**

   - Coding sequences
   - Regulatory elements
   - Expression mechanisms
   - Context for annotation

4. **Evolution:**

   - Mutations, selection, drift
   - Horizontal gene transfer
   - Comparative genomics
   - Explains diversity in metagenomes

5. **Microbial Genomes:**

   - Small, gene-dense
   - Accessory genes variable
   - Pangenome concept
   - Relevant for metagenomic binning

6. **Functional Genomics:**

   - Gene expression
   - Metabolic pathways
   - Phenotype connections
   - Validates binning results

7. **Technologies:**
   - Sequencing advances enable metagenomics
   - Assembly creates data for binning
   - Annotation interprets bins
   - CRISPR and synthetic biology applications

**For GraphBin:**

- **Genome continuity:** Justifies using assembly graphs
- **Compositional signatures:** Why TNF works
- **Gene clustering:** Why connected contigs likely from same genome
- **HGT and repeats:** Challenges GraphBin addresses
- **Validation:** Genomic markers assess bin quality

**Integration:**

```
Genomics ──────────┐
                   ├─→ Metagenomics
Assembly ──────────┤
                   ├─→ GraphBin
Binning ───────────┘
```

Understanding genomics is essential for:

- Interpreting metagenomic data
- Developing binning algorithms
- Validating results
- Making biological sense of bins

---

## Additional Resources

### Textbooks

**General Genomics:**

- Brown, "Genomes" (comprehensive introduction)
- Lesk, "Introduction to Genomics"
- Primrose & Twyman, "Genomics: Applications in Human Biology"

**Microbial Genomics:**

- Fraser et al., "Microbial Genomes"
- Bacterial Genomes volumes in Methods in Molecular Biology

**Molecular Biology:**

- Alberts et al., "Molecular Biology of the Cell"
- Lodish et al., "Molecular Cell Biology"

### Online Resources

**Databases:**

- NCBI: https://www.ncbi.nlm.nih.gov/
- Ensembl: https://www.ensembl.org/
- UniProt: https://www.uniprot.org/
- KEGG: https://www.genome.jp/kegg/

**Tools:**

- BLAST: https://blast.ncbi.nlm.nih.gov/
- Galaxy: https://usegalaxy.org/ (web-based analysis)
- Bioconda: https://bioconda.github.io/ (bioinformatics software)

**Learning:**

- Rosalind: http://rosalind.info/ (bioinformatics problems)
- NCBI Education: https://www.ncbi.nlm.nih.gov/home/learn/
- EBI Training: https://www.ebi.ac.uk/training

### Key Papers

**Foundational:**

- Watson & Crick (1953): DNA structure
- Sanger et al. (1977): First genome sequence
- Fleischmann et al. (1995): First free-living organism genome

**Human Genome:**

- Lander et al. (2001): Initial sequencing and analysis
- Venter et al. (2001): The sequence of the human genome

**Metagenomics:**

- Handelsman et al. (1998): Coined "metagenomics"
- Tyson et al. (2004): Community structure from genomic data
- Venter et al. (2004): Sargasso Sea metagenome

**Pangenomics:**

- Tettelin et al. (2005): Introduced pangenome concept
- Medini et al. (2005): The microbial pan-genome

### Courses

**Coursera:**

- "Genomic Data Science Specialization"
- "Bioinformatics Specialization"

**edX:**

- "Introduction to Genomics"
- "Principles of Synthetic Biology"

**MIT OpenCourseWare:**

- 7.012 Introduction to Biology
- 7.91J Foundations of Computational and Systems Biology

---

**Last Updated:** December 2025  
**Part of:** GraphBin Developer Documentation  
**Complements:** domain_metagenomics.md, domain_assembly.md, domain_binning.md
