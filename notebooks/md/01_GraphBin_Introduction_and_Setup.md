# GraphBin Tutorial 1: Introduction and Setup

## Overview

**GraphBin** is a metagenomic contig bin refinement tool that uses assembly graph connectivity information to improve binning results. This notebook introduces GraphBin and demonstrates basic setup and verification.

### What GraphBin Does:

1. **Refines existing binning results** - Takes binning output from tools like MaxBin2 or MetaBAT2 and improves it
2. **Uses assembly graph connectivity** - Leverages the fact that connected contigs likely come from the same organism
3. **Recovers unbinned sequences** - Assigns bins to short contigs that other tools couldn't classify
4. **Fixes mis-binned contigs** - Corrects contigs that were assigned to the wrong bin

### Key Concepts:

- **Contig**: A contiguous DNA sequence assembled from shorter reads
- **Bin**: A group of contigs believed to come from the same organism
- **Assembly Graph**: A graph showing how contigs are connected during assembly
- **Label Propagation**: Algorithm that spreads bin labels along connected contigs

---

## 1. Installation Verification

Let's verify that GraphBin is properly installed and check the version.

```python
# Check GraphBin installation
import sys
import os

# Add the src directory to the path if running from the repository
repo_root = os.path.abspath('..')
src_path = os.path.join(repo_root, 'src')
if src_path not in sys.path:
    sys.path.insert(0, src_path)

# Import GraphBin modules
import graphbin
from graphbin import __version__

print(f"GraphBin version: {__version__}")
print(f"GraphBin location: {graphbin.__file__}")
```

## 2. Check Dependencies

GraphBin requires several Python packages. Let's verify they're installed.

```python
# Check required dependencies
dependencies = {
    'igraph': None,
    'cogent3': None,
    'cairocffi': None,
    'click': None
}

for package in dependencies.keys():
    try:
        mod = __import__(package)
        version = getattr(mod, '__version__', 'unknown')
        dependencies[package] = version
        print(f"✓ {package}: {version}")
    except ImportError:
        print(f"✗ {package}: NOT INSTALLED")
        dependencies[package] = None
```

## 3. Explore Test Data Structure

GraphBin comes with test datasets. Let's explore what files are available.

```python
import os
from pathlib import Path

# Define test data directory
test_data_dir = Path(repo_root) / 'tests' / 'data'

print(f"Test data directory: {test_data_dir}")
print(f"\nAvailable test datasets:\n")

# List all test datasets
if test_data_dir.exists():
    for dataset in sorted(test_data_dir.iterdir()):
        if dataset.is_dir():
            print(f"📁 {dataset.name}")
            # List files in each dataset
            for file in sorted(dataset.iterdir()):
                size = file.stat().st_size if file.is_file() else 0
                size_str = f"{size:,} bytes" if size < 1024*1024 else f"{size/(1024*1024):.2f} MB"
                print(f"  └─ {file.name} ({size_str})")
            print()
else:
    print(f"⚠️  Test data directory not found: {test_data_dir}")
```

## 4. Understanding Input File Formats

GraphBin requires specific input files. Let's examine their formats using the test data.

```python
# Use ESC_metaSPAdes test dataset as example
test_dataset = test_data_dir / 'ESC_metaSPAdes'

print("=" * 80)
print("REQUIRED INPUT FILES FOR GRAPHBIN")
print("=" * 80)

files_info = {
    'contigs.fasta': 'FASTA file containing assembled contigs',
    'assembly_graph_with_scaffolds.gfa': 'Assembly graph in GFA format',
    'contigs.paths': 'Paths file mapping contigs to graph (SPAdes only)',
    'initial_binning_res.csv': 'Initial binning results from another tool'
}

for filename, description in files_info.items():
    filepath = test_dataset / filename
    print(f"\n{'='*80}")
    print(f"File: {filename}")
    print(f"Description: {description}")
    print(f"{'='*80}")
    
    if filepath.exists():
        print(f"✓ Found: {filepath}")
        print(f"  Size: {filepath.stat().st_size:,} bytes")
        
        # Show first few lines
        with open(filepath, 'r') as f:
            lines = [next(f) for _ in range(min(5, sum(1 for _ in open(filepath))))]  
        print(f"\n  First {len(lines)} lines:")
        for i, line in enumerate(lines, 1):
            print(f"  {i}: {line.rstrip()[:100]}{'...' if len(line.rstrip()) > 100 else ''}")
    else:
        print(f"✗ Not found: {filepath}")
```

## 5. Parse and Understand Contig FASTA Format

Let's examine the contigs file in detail to understand the sequence data.

```python
from cogent3 import load_seq

contigs_file = test_dataset / 'contigs.fasta'

# Read FASTA file and get basic statistics
sequences = {}
current_seq = None
current_name = None

with open(contigs_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if current_name:
                sequences[current_name] = current_seq
            current_name = line[1:].split()[0]  # Get contig name
            current_seq = ''
        else:
            current_seq += line
    if current_name:
        sequences[current_name] = current_seq

print(f"Total contigs: {len(sequences)}")
print(f"\nContig statistics:")

lengths = [len(seq) for seq in sequences.values()]
print(f"  Min length: {min(lengths):,} bp")
print(f"  Max length: {max(lengths):,} bp")
print(f"  Mean length: {sum(lengths)/len(lengths):,.0f} bp")
print(f"  Total bases: {sum(lengths):,} bp")

# Show first 3 contigs
print(f"\nFirst 3 contigs:")
for i, (name, seq) in enumerate(list(sequences.items())[:3], 1):
    print(f"\n{i}. {name}")
    print(f"   Length: {len(seq):,} bp")
    print(f"   Sequence: {seq[:60]}...")
```

## 6. Parse Initial Binning Results

The initial binning CSV file shows how contigs were grouped before GraphBin refinement.

```python
import csv
from collections import Counter

binning_file = test_dataset / 'initial_binning_res.csv'

# Parse binning results
contig_bins = {}
with open(binning_file, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader:
        if len(row) >= 2:
            contig_bins[row[0]] = row[1]

print(f"Total binned contigs: {len(contig_bins)}")
print(f"Total contigs in assembly: {len(sequences)}")
print(f"Unbinned contigs: {len(sequences) - len(contig_bins)}")

# Count bins
bin_counts = Counter(contig_bins.values())
print(f"\nNumber of bins: {len(bin_counts)}")
print(f"\nContigs per bin:")
for bin_name, count in sorted(bin_counts.items()):
    print(f"  {bin_name}: {count} contigs")

# Show some examples
print(f"\nExample binning assignments:")
for i, (contig, bin_id) in enumerate(list(contig_bins.items())[:5], 1):
    length = len(sequences.get(contig, ''))
    print(f"  {i}. {contig} → {bin_id} (length: {length:,} bp)")
```

## 7. Understanding Assembly Graph Structure

The assembly graph shows how contigs connect to each other. Let's examine the GFA format.

```python
graph_file = test_dataset / 'assembly_graph_with_scaffolds.gfa'

# Parse GFA file
segments = []  # S lines (contigs/nodes)
links = []     # L lines (connections/edges)

with open(graph_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('S'):
            segments.append(line)
        elif line.startswith('L'):
            links.append(line)

print("GFA File Structure:")
print(f"  Segments (S): {len(segments)} nodes (contigs)")
print(f"  Links (L): {len(links)} edges (connections)")

print(f"\nExample segment lines (first 3):")
for i, seg in enumerate(segments[:3], 1):
    parts = seg.split('\t')
    print(f"  {i}. Node: {parts[1]}, Length: {len(parts[2]) if len(parts) > 2 else 'N/A'}")
    print(f"     {seg[:100]}...")

print(f"\nExample link lines (first 5):")
for i, link in enumerate(links[:5], 1):
    parts = link.split('\t')
    if len(parts) >= 5:
        print(f"  {i}. {parts[1]} ({parts[2]}) → {parts[3]} ({parts[4]})")
    else:
        print(f"  {i}. {link}")
```

## 8. Load Assembly Graph with igraph

GraphBin uses the igraph library to represent and analyze the assembly graph.

```python
import igraph as ig

# Create a simple example graph to demonstrate
print("Creating example graph structure...\n")

# Simple example with 5 nodes
g = ig.Graph()
g.add_vertices(5)
g.add_edges([(0, 1), (1, 2), (2, 3), (3, 4), (1, 3)])

# Add some attributes
g.vs['name'] = ['NODE_1', 'NODE_2', 'NODE_3', 'NODE_4', 'NODE_5']
g.vs['length'] = [1000, 1500, 2000, 1200, 800]

print(f"Graph properties:")
print(f"  Vertices: {g.vcount()}")
print(f"  Edges: {g.ecount()}")
print(f"  Average degree: {sum(g.degree()) / g.vcount():.2f}")

print(f"\nVertex information:")
for v in g.vs:
    neighbors = g.neighbors(v.index)
    print(f"  {v['name']}: degree={g.degree(v.index)}, neighbors={[g.vs[n]['name'] for n in neighbors]}")

print(f"\nEdge list:")
for e in g.es:
    src = g.vs[e.source]['name']
    dst = g.vs[e.target]['name']
    print(f"  {src} ↔ {dst}")
```

## 9. GraphBin Module Structure

Let's explore the main modules and their purposes.

```python
import graphbin
from graphbin import (
    graphbin_SPAdes,
    graphbin_SGA,
    graphbin_MEGAHIT,
    graphbin_Flye,
    graphbin_Canu,
    graphbin_Miniasm,
    graphbin_Func
)

modules_info = {
    'graphbin_SPAdes': 'Handles SPAdes/metaSPAdes assemblies (GFA format with paths)',
    'graphbin_SGA': 'Handles SGA assemblies (ASQG format)',
    'graphbin_MEGAHIT': 'Handles MEGAHIT assemblies (GFA format)',
    'graphbin_Flye': 'Handles Flye long-read assemblies (GFA format)',
    'graphbin_Canu': 'Handles Canu long-read assemblies (GFA format)',
    'graphbin_Miniasm': 'Handles Miniasm long-read assemblies (GFA format)',
    'graphbin_Func': 'Core algorithm functions (label propagation)'
}

print("GraphBin Module Structure:\n")
print("=" * 80)
for module_name, description in modules_info.items():
    module = eval(module_name)
    print(f"\n📦 {module_name}")
    print(f"   {description}")
    
    # List main functions
    functions = [name for name in dir(module) if callable(getattr(module, name)) and not name.startswith('_')]
    if functions:
        print(f"   Functions: {', '.join(functions[:5])}{'...' if len(functions) > 5 else ''}")
```

## 10. Command-Line Interface

GraphBin provides a CLI through the `click` library. Let's examine the available options.

```python
from graphbin.cli import main as graphbin_main
from click.testing import CliRunner

# Use CliRunner to capture help output
runner = CliRunner()
result = runner.invoke(graphbin_main, ['--help'])

print("GraphBin Command-Line Interface:\n")
print(result.output)
```

## 11. Summary and Next Steps

In this notebook, we:

1. ✅ Verified GraphBin installation
2. ✅ Checked dependencies
3. ✅ Explored test data structure
4. ✅ Understood input file formats:
   - Contigs (FASTA)
   - Assembly graph (GFA)
   - Initial binning (CSV)
5. ✅ Examined assembly graph structure
6. ✅ Learned about GraphBin modules

### Next Steps:

- **Notebook 2**: Core functionality - Run GraphBin refinement
- **Notebook 3**: Different assemblers - Compare SPAdes, MEGAHIT, SGA, etc.
- **Notebook 4**: Advanced features - Label propagation algorithm details
- **Notebook 5**: Utilities - Parsers, visualization, and output analysis

---

### Key Takeaways:

- GraphBin refines existing binning results using assembly graph connectivity
- Requires 3-4 input files: contigs, graph, binning results, (optional) paths
- Supports 6 different assemblers
- Uses label propagation algorithm on the assembly graph
- Built with Python, igraph, and cogent3
