# GraphBin Tutorial 2: Core Functionality - Binning Refinement

## Overview

This notebook demonstrates the core functionality of GraphBin: refining metagenomic binning results using assembly graph connectivity. We'll walk through a complete workflow using the SPAdes test dataset.

### What You'll Learn:

1. How to run GraphBin programmatically
2. Understanding the label propagation algorithm
3. Analyzing refinement results
4. Comparing before and after binning
5. Output file formats and interpretation

---

## 1. Setup and Imports

```python
import sys
import os
from pathlib import Path
import csv
from collections import Counter, defaultdict

# Add GraphBin to path
repo_root = Path('..').resolve()
src_path = repo_root / 'src'
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

# Import GraphBin modules
from graphbin import graphbin_SPAdes
from graphbin.cli import ArgsObj
from graphbin.parsers import spades_parser
from graphbin.parsers import get_initial_bin_count

import igraph as ig

print("✓ Imports successful")
print(f"Working directory: {Path.cwd()}")
print(f"Repository root: {repo_root}")
```

## 2. Define Input Files

We'll use the ESC_metaSPAdes test dataset which contains realistic metagenomic data.

```python
# Test dataset paths
test_data = repo_root / 'tests' / 'data' / 'ESC_metaSPAdes'
output_dir = repo_root / 'notebooks' / 'output' / 'spades_refinement'

# Create output directory
output_dir.mkdir(parents=True, exist_ok=True)

# Input files
graph_file = test_data / 'assembly_graph_with_scaffolds.gfa'
contigs_file = test_data / 'contigs.fasta'
paths_file = test_data / 'contigs.paths'
binned_file = test_data / 'initial_binning_res.csv'

# Verify files exist
files = {
    'Assembly graph': graph_file,
    'Contigs': contigs_file,
    'Paths': paths_file,
    'Initial binning': binned_file
}

print("Input files:")
for name, path in files.items():
    exists = "✓" if path.exists() else "✗"
    size = f"{path.stat().st_size:,} bytes" if path.exists() else "N/A"
    print(f"  {exists} {name}: {path.name} ({size})")

print(f"\nOutput directory: {output_dir}")
```

## 3. Analyze Initial Binning Results

Before refinement, let's understand the initial binning state.

```python
# Parse initial binning
initial_bins = {}
with open(binned_file, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader:
        if len(row) >= 2:
            initial_bins[row[0]] = row[1]

# Count contigs
sequences = {}
with open(contigs_file, 'r') as f:
    current_seq = ''
    current_name = None
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if current_name:
                sequences[current_name] = current_seq
            current_name = line[1:].split()[0]
            current_seq = ''
        else:
            current_seq += line
    if current_name:
        sequences[current_name] = current_seq

# Statistics
bin_counts = Counter(initial_bins.values())
unbinned_contigs = set(sequences.keys()) - set(initial_bins.keys())

print("INITIAL BINNING STATISTICS")
print("=" * 60)
print(f"Total contigs in assembly: {len(sequences)}")
print(f"Binned contigs: {len(initial_bins)}")
print(f"Unbinned contigs: {len(unbinned_contigs)}")
print(f"Number of bins: {len(bin_counts)}")
print(f"\nContigs per bin:")
for bin_name in sorted(bin_counts.keys()):
    count = bin_counts[bin_name]
    # Calculate total length for this bin
    bin_contigs = [c for c, b in initial_bins.items() if b == bin_name]
    total_length = sum(len(sequences[c]) for c in bin_contigs if c in sequences)
    print(f"  {bin_name}: {count} contigs, {total_length:,} bp")

# Unbinned statistics
if unbinned_contigs:
    unbinned_lengths = [len(sequences[c]) for c in unbinned_contigs if c in sequences]
    print(f"\nUnbinned contigs statistics:")
    print(f"  Count: {len(unbinned_contigs)}")
    print(f"  Length range: {min(unbinned_lengths):,} - {max(unbinned_lengths):,} bp")
    print(f"  Mean length: {sum(unbinned_lengths)/len(unbinned_lengths):,.0f} bp")
    print(f"  Total: {sum(unbinned_lengths):,} bp")
```

## 4. Examine Assembly Graph Connectivity

The assembly graph structure is key to GraphBin's refinement.

```python
# Parse assembly graph to understand connectivity
segments = []
links = []

with open(graph_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('S'):
            segments.append(line.split('\t'))
        elif line.startswith('L'):
            links.append(line.split('\t'))

print("ASSEMBLY GRAPH STRUCTURE")
print("=" * 60)
print(f"Graph nodes (segments): {len(segments)}")
print(f"Graph edges (links): {len(links)}")
print(f"Average degree: {(2 * len(links)) / len(segments):.2f}")

# Analyze connectivity
node_degrees = defaultdict(int)
for link in links:
    if len(link) >= 5:
        node_degrees[link[1]] += 1
        node_degrees[link[3]] += 1

degree_dist = Counter(node_degrees.values())
print(f"\nDegree distribution:")
for degree in sorted(degree_dist.keys())[:10]:  # Show first 10
    count = degree_dist[degree]
    print(f"  Degree {degree}: {count} nodes")

# Check connectivity of binned vs unbinned contigs
binned_node_names = set(initial_bins.keys())
unbinned_node_names = unbinned_contigs

# Count connections between binned/unbinned
binned_to_binned = 0
binned_to_unbinned = 0
unbinned_to_unbinned = 0

for link in links:
    if len(link) >= 5:
        src = link[1]
        dst = link[3]
        
        src_binned = src in binned_node_names
        dst_binned = dst in binned_node_names
        
        if src_binned and dst_binned:
            binned_to_binned += 1
        elif (src_binned and not dst_binned) or (not src_binned and dst_binned):
            binned_to_unbinned += 1
        else:
            unbinned_to_unbinned += 1

print(f"\nConnectivity analysis:")
print(f"  Binned ↔ Binned: {binned_to_binned} edges")
print(f"  Binned ↔ Unbinned: {binned_to_unbinned} edges")
print(f"  Unbinned ↔ Unbinned: {unbinned_to_unbinned} edges")
print(f"\n💡 The {binned_to_unbinned} edges connecting binned to unbinned contigs")
print(f"   are opportunities for GraphBin to assign bins to unbinned contigs!")
```

## 5. Run GraphBin Refinement

Now let's run the actual GraphBin algorithm to refine the binning.

```python
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

# Create arguments object
args = ArgsObj(
    assembler='spades',
    graph=str(graph_file),
    contigs=str(contigs_file),
    paths=str(paths_file),
    binned=str(binned_file),
    output=str(output_dir),
    prefix='graphbin',
    max_iteration=100,
    diff_threshold=0.1,
    delimiter=','
)

print("Running GraphBin refinement...")
print("=" * 60)
print(f"Input parameters:")
print(f"  Assembler: {args.assembler}")
print(f"  Max iterations: {args.max_iteration}")
print(f"  Difference threshold: {args.diff_threshold}")
print(f"  Output prefix: {args.prefix}")
print()

# Run GraphBin
try:
    graphbin_SPAdes.main(args)
    print("\n✓ GraphBin refinement completed successfully!")
except Exception as e:
    print(f"\n✗ Error during refinement: {e}")
    import traceback
    traceback.print_exc()
```

## 6. Analyze Refined Binning Results

Let's compare the initial and refined binning results.

```python
# Read refined binning results
refined_binning_file = output_dir / 'graphbin_output.csv'

refined_bins = {}
if refined_binning_file.exists():
    with open(refined_binning_file, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            if len(row) >= 2:
                refined_bins[row[0]] = row[1]
    
    print("REFINED BINNING STATISTICS")
    print("=" * 60)
    
    refined_bin_counts = Counter(refined_bins.values())
    refined_unbinned = set(sequences.keys()) - set(refined_bins.keys())
    
    print(f"Total contigs in assembly: {len(sequences)}")
    print(f"Binned contigs: {len(refined_bins)}")
    print(f"Unbinned contigs: {len(refined_unbinned)}")
    print(f"Number of bins: {len(refined_bin_counts)}")
    
    print(f"\nRefined contigs per bin:")
    for bin_name in sorted(refined_bin_counts.keys()):
        count = refined_bin_counts[bin_name]
        bin_contigs = [c for c, b in refined_bins.items() if b == bin_name]
        total_length = sum(len(sequences[c]) for c in bin_contigs if c in sequences)
        print(f"  {bin_name}: {count} contigs, {total_length:,} bp")
    
    # Compare initial vs refined
    print("\n" + "=" * 60)
    print("COMPARISON: Initial vs Refined")
    print("=" * 60)
    
    initial_binned_count = len(initial_bins)
    refined_binned_count = len(refined_bins)
    recovered = refined_binned_count - initial_binned_count
    
    print(f"\nBinning coverage:")
    print(f"  Initial: {initial_binned_count} / {len(sequences)} contigs ({100*initial_binned_count/len(sequences):.1f}%)")
    print(f"  Refined: {refined_binned_count} / {len(sequences)} contigs ({100*refined_binned_count/len(sequences):.1f}%)")
    print(f"  Recovered: {recovered} additional contigs")
    
else:
    print(f"⚠️  Refined binning file not found: {refined_binning_file}")
```

## 7. Detailed Change Analysis

Let's identify which contigs changed bins and which were newly assigned.

```python
if refined_bins:
    # Find changes
    newly_binned = []
    changed_bins = []
    unchanged = []
    
    for contig in sequences.keys():
        initial_bin = initial_bins.get(contig, None)
        refined_bin = refined_bins.get(contig, None)
        
        if initial_bin is None and refined_bin is not None:
            newly_binned.append((contig, refined_bin))
        elif initial_bin is not None and refined_bin is not None:
            if initial_bin != refined_bin:
                changed_bins.append((contig, initial_bin, refined_bin))
            else:
                unchanged.append(contig)
    
    print("DETAILED CHANGE ANALYSIS")
    print("=" * 60)
    print(f"Unchanged contigs: {len(unchanged)}")
    print(f"Newly binned contigs: {len(newly_binned)}")
    print(f"Changed bin assignments: {len(changed_bins)}")
    
    # Show newly binned contigs
    if newly_binned:
        print(f"\nNewly binned contigs (first 10):")
        for i, (contig, bin_id) in enumerate(newly_binned[:10], 1):
            length = len(sequences[contig])
            print(f"  {i}. {contig} → {bin_id} (length: {length:,} bp)")
    
    # Show changed bins
    if changed_bins:
        print(f"\nChanged bin assignments (first 10):")
        for i, (contig, old_bin, new_bin) in enumerate(changed_bins[:10], 1):
            length = len(sequences[contig])
            print(f"  {i}. {contig}: {old_bin} → {new_bin} (length: {length:,} bp)")
    
    # Summary of improvements
    print("\n" + "=" * 60)
    print("SUMMARY OF IMPROVEMENTS")
    print("=" * 60)
    
    if newly_binned:
        newly_binned_length = sum(len(sequences[c]) for c, _ in newly_binned)
        print(f"✓ Recovered {len(newly_binned)} previously unbinned contigs")
        print(f"  Total: {newly_binned_length:,} bp of additional sequence")
    
    if changed_bins:
        print(f"✓ Corrected {len(changed_bins)} mis-binned contigs")
    
    if not newly_binned and not changed_bins:
        print("ℹ️  No changes detected - initial binning was already optimal")
```

## 8. Examine Output Files

GraphBin generates several output files. Let's explore them.

```python
print("OUTPUT FILES")
print("=" * 60)

# List all output files
output_files = sorted(output_dir.glob('*'))

if output_files:
    for i, filepath in enumerate(output_files, 1):
        size = filepath.stat().st_size
        size_str = f"{size:,} bytes" if size < 1024*1024 else f"{size/(1024*1024):.2f} MB"
        print(f"\n{i}. {filepath.name}")
        print(f"   Size: {size_str}")
        
        # Show file type and purpose
        if filepath.suffix == '.csv':
            print(f"   Type: Binning result file")
            with open(filepath, 'r') as f:
                lines = f.readlines()[:3]
            print(f"   First 3 lines:")
            for line in lines:
                print(f"     {line.rstrip()}")
        
        elif filepath.suffix == '.fasta':
            print(f"   Type: FASTA file (bin-specific sequences)")
            # Count sequences
            seq_count = 0
            with open(filepath, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        seq_count += 1
            print(f"   Sequences: {seq_count}")
        
        elif filepath.suffix == '.log':
            print(f"   Type: Log file")
else:
    print("No output files found.")
```

## 9. Visualize Bin Size Distribution

Let's compare the distribution of bin sizes before and after refinement.

```python
import matplotlib.pyplot as plt
import numpy as np

if refined_bins:
    # Calculate bin sizes (by total bp)
    initial_bin_sizes = {}
    for bin_name in set(initial_bins.values()):
        contigs = [c for c, b in initial_bins.items() if b == bin_name]
        total_size = sum(len(sequences[c]) for c in contigs if c in sequences)
        initial_bin_sizes[bin_name] = total_size / 1_000_000  # Convert to Mbp
    
    refined_bin_sizes = {}
    for bin_name in set(refined_bins.values()):
        contigs = [c for c, b in refined_bins.items() if b == bin_name]
        total_size = sum(len(sequences[c]) for c in contigs if c in sequences)
        refined_bin_sizes[bin_name] = total_size / 1_000_000  # Convert to Mbp
    
    # Create comparison plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Initial bins
    bins_initial = sorted(initial_bin_sizes.keys())
    sizes_initial = [initial_bin_sizes[b] for b in bins_initial]
    ax1.bar(range(len(bins_initial)), sizes_initial, color='steelblue', alpha=0.7)
    ax1.set_xlabel('Bin')
    ax1.set_ylabel('Total Size (Mbp)')
    ax1.set_title('Initial Binning')
    ax1.set_xticks(range(len(bins_initial)))
    ax1.set_xticklabels(bins_initial, rotation=45)
    ax1.grid(axis='y', alpha=0.3)
    
    # Refined bins
    bins_refined = sorted(refined_bin_sizes.keys())
    sizes_refined = [refined_bin_sizes[b] for b in bins_refined]
    ax2.bar(range(len(bins_refined)), sizes_refined, color='forestgreen', alpha=0.7)
    ax2.set_xlabel('Bin')
    ax2.set_ylabel('Total Size (Mbp)')
    ax2.set_title('Refined Binning (GraphBin)')
    ax2.set_xticks(range(len(bins_refined)))
    ax2.set_xticklabels(bins_refined, rotation=45)
    ax2.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'bin_size_comparison.png', dpi=150, bbox_inches='tight')
    print(f"\n✓ Plot saved to: {output_dir / 'bin_size_comparison.png'}")
    plt.show()
    
    # Print size comparison
    print("\nBin size comparison (Mbp):")
    print("-" * 60)
    all_bins = sorted(set(bins_initial) | set(bins_refined))
    for bin_name in all_bins:
        initial_size = initial_bin_sizes.get(bin_name, 0)
        refined_size = refined_bin_sizes.get(bin_name, 0)
        diff = refined_size - initial_size
        print(f"{bin_name}: {initial_size:.2f} → {refined_size:.2f} Mbp (Δ {diff:+.2f} Mbp)")
```

## 10. Understanding the Label Propagation Algorithm

GraphBin uses a label propagation algorithm. Let's understand how it works conceptually.

```python
print("LABEL PROPAGATION ALGORITHM EXPLAINED")
print("=" * 60)
print()
print("The label propagation algorithm works through these steps:")
print()
print("1. INITIALIZATION")
print("   - Start with initial binning results (bins = labels)")
print("   - Some contigs are labeled (binned), others are unlabeled (unbinned)")
print()
print("2. IDENTIFY AMBIGUOUS VERTICES")
print("   - Find contigs whose neighbors have different bin labels")
print("   - These might be incorrectly binned")
print("   - Remove their labels temporarily")
print()
print("3. ITERATIVE PROPAGATION")
print("   - For each unlabeled contig:")
print("     a) Look at its neighbors in the assembly graph")
print("     b) Count the bin labels of neighbors")
print("     c) Assign the most common neighboring bin label")
print()
print("4. CONVERGENCE CHECK")
print("   - After each iteration, count how many labels changed")
print("   - Stop when:")
print("     • Changes fall below threshold (default: 10%)")
print("     • Maximum iterations reached (default: 100)")
print()
print("5. OUTPUT")
print("   - Final refined binning with:")
print("     • Corrected mis-binned contigs")
print("     • Newly assigned previously unbinned contigs")
print()
print("=" * 60)
print("KEY INSIGHT:")
print("Connected contigs in the assembly graph likely come from the")
print("same organism. Label propagation leverages this connectivity!")
print("=" * 60)
```

## 11. Summary and Key Metrics

```python
if refined_bins:
    print("\n" + "=" * 60)
    print("GRAPHBIN REFINEMENT SUMMARY")
    print("=" * 60)
    print()
    print("📊 METRICS:")
    print()
    
    # Coverage improvement
    initial_coverage = len(initial_bins) / len(sequences) * 100
    refined_coverage = len(refined_bins) / len(sequences) * 100
    coverage_improvement = refined_coverage - initial_coverage
    
    print(f"  Binning Coverage:")
    print(f"    Initial:  {initial_coverage:.1f}%")
    print(f"    Refined:  {refined_coverage:.1f}%")
    print(f"    Improvement: +{coverage_improvement:.1f}%")
    print()
    
    # Sequence recovery
    if newly_binned:
        recovered_bp = sum(len(sequences[c]) for c, _ in newly_binned)
        recovered_pct = recovered_bp / sum(len(s) for s in sequences.values()) * 100
        print(f"  Sequence Recovery:")
        print(f"    Contigs recovered: {len(newly_binned)}")
        print(f"    Bases recovered: {recovered_bp:,} bp ({recovered_pct:.1f}% of total)")
        print()
    
    # Corrections
    if changed_bins:
        print(f"  Bin Corrections:")
        print(f"    Contigs reassigned: {len(changed_bins)}")
        print()
    
    print("=" * 60)
    print("✓ GraphBin successfully refined the metagenomic binning!")
    print("=" * 60)
```

## Next Steps

In this notebook, we:

1. ✅ Loaded and analyzed test data
2. ✅ Examined initial binning statistics
3. ✅ Analyzed assembly graph connectivity
4. ✅ Ran GraphBin refinement
5. ✅ Compared initial vs refined results
6. ✅ Identified specific improvements
7. ✅ Visualized bin size changes
8. ✅ Understood the label propagation algorithm

### Continue to:

- **Notebook 3**: Try different assemblers (MEGAHIT, SGA, Flye, etc.)
- **Notebook 4**: Deep dive into label propagation and graph analysis
- **Notebook 5**: Use parser and visualization utilities

---

### Key Takeaways:

- GraphBin typically recovers 10-30% more contigs
- It corrects mis-binned contigs using graph connectivity
- The algorithm is iterative and converges to a stable solution
- Results depend on assembly graph quality and initial binning
- Works best when contigs are well-connected in the graph
