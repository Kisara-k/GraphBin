# GraphBin Tutorial 5: Parsers, Utilities, and Visualization

## Overview

GraphBin includes several utility modules for:
- Parsing different assembly graph formats
- Converting between file formats
- Visualizing binning results
- Processing and analyzing outputs

This notebook demonstrates how to use these utilities programmatically.

---

## 1. Setup and Imports

```python
import sys
import os
from pathlib import Path
import csv
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import numpy as np

# Add GraphBin to path
repo_root = Path('..').resolve()
src_path = repo_root / 'src'
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

# Import GraphBin parsers
from graphbin.parsers import (
    spades_parser,
    megahit_parser,
    sga_parser,
    flye_parser,
    canu_parser,
    miniasm_parser,
    get_initial_bin_count
)

# Import support utilities
from graphbin.support import gfa2fasta, prep_result

import igraph as ig

print("✓ Imports successful")
print(f"Repository root: {repo_root}")
```

## 2. Understanding Parser Modules

GraphBin has specialized parsers for each assembler type.

```python
print("GRAPHBIN PARSER MODULES")
print("=" * 80)
print()

parsers_info = [
    {
        'Module': 'spades_parser',
        'Assembler': 'SPAdes/metaSPAdes',
        'Graph Format': 'GFA',
        'Key Functions': [
            'parse_graph(graph_file, contig_paths)',
            'get_initial_binning_result(binned_file, delimiter)',
            'write_output(output_path, contig_bins, ...)',
        ]
    },
    {
        'Module': 'megahit_parser',
        'Assembler': 'MEGAHIT',
        'Graph Format': 'GFA',
        'Key Functions': [
            'parse_graph(graph_file, original_contigs)',
            'get_contig_descriptors(contigs_file)',
            'write_output(output_path, contig_bins, ...)',
        ]
    },
    {
        'Module': 'sga_parser',
        'Assembler': 'SGA',
        'Graph Format': 'ASQG',
        'Key Functions': [
            'parse_graph(graph_file)',
            'get_contig_descriptions(contigs_file)',
            'write_output(output_path, contig_bins, ...)',
        ]
    },
    {
        'Module': 'flye_parser',
        'Assembler': 'Flye',
        'Graph Format': 'GFA',
        'Key Functions': [
            'parse_graph(graph_file, contig_paths)',
            'get_initial_binning_result(binned_file, delimiter)',
            'write_output(output_path, contig_bins, ...)',
        ]
    },
    {
        'Module': 'canu_parser',
        'Assembler': 'Canu',
        'Graph Format': 'GFA',
        'Key Functions': [
            'parse_graph(graph_file)',
            'get_initial_binning_result(binned_file, delimiter)',
            'write_output(output_path, contig_bins, ...)',
        ]
    },
    {
        'Module': 'miniasm_parser',
        'Assembler': 'Miniasm',
        'Graph Format': 'GFA',
        'Key Functions': [
            'parse_graph(graph_file)',
            'get_initial_binning_result(binned_file, delimiter)',
            'write_output(output_path, contig_bins, ...)',
        ]
    }
]

for i, parser in enumerate(parsers_info, 1):
    print(f"{i}. {parser['Module']}")
    print(f"   Assembler: {parser['Assembler']}")
    print(f"   Graph Format: {parser['Graph Format']}")
    print(f"   Key Functions:")
    for func in parser['Key Functions']:
        print(f"     • {func}")
    print()

print("=" * 80)
print("COMMON PARSER FUNCTIONS (all modules):")
print("=" * 80)
print()
print("parse_graph()")
print("  • Reads assembly graph file (GFA or ASQG)")
print("  • Creates igraph Graph object")
print("  • Returns: (graph, contig_names, node_count)")
print()
print("get_initial_binning_result()")
print("  • Reads initial binning CSV file")
print("  • Returns: (contig_bins, bins, bins_list)")
print()
print("write_output()")
print("  • Writes refined binning results")
print("  • Creates per-bin FASTA files")
print("  • Generates output CSV")
print()
print("=" * 80)
```

## 3. Parsing Assembly Graphs

Let's use the parsers to load and analyze different graph formats.

```python
# Example: Parse SPAdes assembly graph
test_data = repo_root / 'tests' / 'data' / 'ESC_metaSPAdes'
graph_file = test_data / 'assembly_graph_with_scaffolds.gfa'
paths_file = test_data / 'contigs.paths'

print("PARSING SPADES ASSEMBLY GRAPH")
print("=" * 80)

# Read contig paths
contig_paths = {}
if paths_file.exists():
    with open(paths_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split()
                if len(parts) >= 2:
                    contig_paths[parts[0]] = parts[1]

print(f"Loaded {len(contig_paths)} contig paths")

# Parse graph
assembly_graph, contig_names, node_count = spades_parser.parse_graph(
    str(graph_file),
    contig_paths
)

print(f"\nGraph statistics:")
print(f"  Vertices: {assembly_graph.vcount()}")
print(f"  Edges: {assembly_graph.ecount()}")
print(f"  Contigs: {len(contig_names)}")
print(f"  Node count: {node_count}")

# Show some contig mappings
print(f"\nExample contig mappings (first 5):")
for i, contig in enumerate(list(contig_names.keys())[:5]):
    node_id = contig_names[contig]
    print(f"  {i+1}. {contig} → Node {node_id}")

print("\n" + "=" * 80)
```

## 4. Parsing Initial Binning Results

All parsers provide functions to read initial binning files.

```python
binned_file = test_data / 'initial_binning_res.csv'

print("PARSING INITIAL BINNING RESULTS")
print("=" * 80)

# Parse using SPAdes parser (same function for all parsers)
contig_bins, bins, bins_list = spades_parser.get_initial_binning_result(
    str(binned_file),
    delimiter=','
)

print(f"Parsed binning results:")
print(f"  Total binned contigs: {len(contig_bins)}")
print(f"  Number of bins: {len(bins)}")
print()

# Show bin contents
print(f"Bin details:")
for bin_id in sorted(bins.keys()):
    contigs_in_bin = bins[bin_id]
    print(f"  Bin {bin_id}: {len(contigs_in_bin)} contigs")

# Show bins_list structure
print(f"\nBins list structure:")
for i, bin_name in enumerate(bins_list[:5]):
    print(f"  {i}. {bin_name}")

print("\n" + "=" * 80)

# Alternative: Use common function
from graphbin.parsers import get_initial_bin_count

n_bins = get_initial_bin_count(str(binned_file), ',')
print(f"\nUsing get_initial_bin_count():")
print(f"  Number of bins: {n_bins}")
print("\n" + "=" * 80)
```

## 5. Using the GFA to FASTA Converter

The `gfa2fasta` utility extracts sequences from GFA files.

```python
print("GFA TO FASTA CONVERSION")
print("=" * 80)
print()
print("The gfa2fasta utility extracts node sequences from GFA files.")
print()
print("Use case: Some assemblers only output GFA files")
print("          gfa2fasta extracts sequences into FASTA format")
print()
print("Usage:")
print("  python -m graphbin.support.gfa2fasta --graph input.gfa --output output.fasta")
print()

# Demonstrate reading GFA segments
segments_with_seq = []
with open(graph_file, 'r') as f:
    for line in f:
        if line.startswith('S'):
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                node_name = parts[1]
                sequence = parts[2]
                if sequence != '*':  # Some GFA files use * for no sequence
                    segments_with_seq.append((node_name, sequence))

print(f"Found {len(segments_with_seq)} segments with sequences in GFA file")
print(f"\nExample segments (first 3):")
for i, (name, seq) in enumerate(segments_with_seq[:3], 1):
    print(f"  {i}. {name}: {len(seq)} bp")
    print(f"     Sequence: {seq[:60]}...")

print("\n" + "=" * 80)
```

## 6. Output Writing Functions

Parsers include functions to write refined binning results.

```python
print("OUTPUT WRITING FUNCTIONS")
print("=" * 80)
print()
print("All parsers include a write_output() function that:")
print()
print("1. Creates refined binning CSV file")
print("   Format: contig_name,bin_id")
print("   Example: NODE_1_length_12345,bin_1")
print()
print("2. Generates per-bin FASTA files")
print("   One file per bin containing all sequences")
print("   Example: graphbin_bin_1.fasta")
print()
print("3. Maintains contig order and metadata")
print()

# Show signature
print("Function signature (SPAdes example):")
print("""
spades_parser.write_output(
    output_path,        # Output directory path
    contig_bins,        # Dict: {contig_name: bin_id}
    contigs_file,       # Path to contigs FASTA
    binned_length,      # Total bp in bins
    prefix="graphbin",  # Output file prefix
    delimiter=","       # CSV delimiter
)
""")

print("=" * 80)

# Demonstrate creating a sample output
output_dir = repo_root / 'notebooks' / 'output' / 'demo_output'
output_dir.mkdir(parents=True, exist_ok=True)

# Create a simple example CSV
example_output = output_dir / 'example_binning.csv'
with open(example_output, 'w', newline='') as f:
    writer = csv.writer(f)
    for contig, bin_id in list(contig_bins.items())[:10]:
        writer.writerow([contig, bin_id])

print(f"\nCreated example output: {example_output}")
print(f"\nFirst 5 lines:")
with open(example_output, 'r') as f:
    for i, line in enumerate(f):
        if i >= 5:
            break
        print(f"  {line.rstrip()}")

print("\n" + "=" * 80)
```

## 7. Visualizing Binning Results

GraphBin includes visualization utilities for different assemblers.

```python
print("VISUALIZATION UTILITIES")
print("=" * 80)
print()
print("GraphBin includes visualization scripts for each assembler:")
print()
print("1. visualise_result_SPAdes.py")
print("   • Plots SPAdes assembly graph with bin colors")
print("   • Shows connectivity between bins")
print()
print("2. visualise_result_MEGAHIT.py")
print("   • Visualizes MEGAHIT assembly graph")
print("   • Color-codes contigs by bin")
print()
print("3. visualise_result_SGA.py")
print("   • Plots SGA string graph")
print("   • Displays overlap relationships")
print()
print("4. visualise_result_Flye_Canu_Miniasm.py")
print("   • Unified visualization for long-read assemblers")
print("   • Shows repeat graph structure")
print()
print("=" * 80)
print("USAGE:")
print("=" * 80)
print()
print("Command-line:")
print("  python visualise_result_SPAdes.py \\")
print("    --graph assembly_graph.gfa \\")
print("    --binned graphbin_output.csv \\")
print("    --output visualization.png")
print()
print("=" * 80)
```

## 8. Custom Visualization Example

Let's create a custom visualization of the binning results.

```python
# Create a visualization of binning results on the assembly graph

# Sample a subgraph for visualization (full graph too large)
sample_size = min(50, assembly_graph.vcount())
sample_indices = np.random.choice(assembly_graph.vcount(), sample_size, replace=False)
subgraph = assembly_graph.subgraph(sample_indices)

# Map nodes to bins
node_colors = []
bin_color_map = {
    'bin.1': 'red',
    'bin.2': 'blue',
    'bin.3': 'green',
    'bin.4': 'orange',
    'unbinned': 'lightgray'
}

for idx in sample_indices:
    # Find contig name for this node
    node_contig = None
    for contig, node_id in contig_names.items():
        if node_id == idx:
            node_contig = contig
            break
    
    if node_contig and node_contig in contig_bins:
        bin_id = contig_bins[node_contig]
        color = bin_color_map.get(bin_id, 'lightgray')
    else:
        color = 'lightgray'
    
    node_colors.append(color)

# Create visualization
fig, ax = plt.subplots(1, 1, figsize=(12, 10))

layout = subgraph.layout('fr')  # Fruchterman-Reingold layout

ig.plot(
    subgraph,
    target=ax,
    layout=layout,
    vertex_size=20,
    vertex_color=node_colors,
    edge_width=0.5,
    edge_color='gray'
)

ax.set_title('Assembly Graph Colored by Bin Assignment\n(Sample of 50 nodes)', 
             fontsize=14, fontweight='bold')
ax.axis('off')

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='red', label='Bin 1'),
    Patch(facecolor='blue', label='Bin 2'),
    Patch(facecolor='green', label='Bin 3'),
    Patch(facecolor='orange', label='Bin 4'),
    Patch(facecolor='lightgray', label='Unbinned')
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

plt.tight_layout()
plt.savefig(repo_root / 'notebooks' / 'output' / 'assembly_graph_binning.png',
            dpi=150, bbox_inches='tight')
print("✓ Assembly graph visualization saved")
plt.show()

print("\n💡 Interpretation:")
print("  • Nodes of same color belong to same bin")
print("  • Connected nodes often share bin labels")
print("  • Gray nodes are unbinned (targets for refinement)")
```

## 9. Comparing Different Graph Formats

Let's compare GFA vs ASQG format parsing.

```python
print("ASSEMBLY GRAPH FORMAT COMPARISON")
print("=" * 80)
print()
print("GFA (Graphical Fragment Assembly) Format:")
print("  Used by: SPAdes, MEGAHIT, Flye, Canu, Miniasm")
print()
print("  Line types:")
print("    S (Segment)  - Graph node with sequence")
print("                   Example: S  NODE_1  ATCG...  LN:i:1234")
print()
print("    L (Link)     - Edge connecting two nodes")
print("                   Example: L  NODE_1  +  NODE_2  +  0M")
print()
print("    P (Path)     - Ordered sequence of nodes")
print("                   Example: P  path1  NODE_1+,NODE_2+  *")
print()
print("  Advantages:")
print("    • Widely supported standard")
print("    • Compact representation")
print("    • Includes optional metadata")
print()
print("=" * 80)
print()
print("ASQG (Assembly String Graph) Format:")
print("  Used by: SGA")
print()
print("  Line types:")
print("    VT (Vertex) - Graph node")
print("                  Example: VT  NODE_1  ATCG...")
print()
print("    ED (Edge)   - Overlap between nodes")
print("                  Example: ED  NODE_1  NODE_2  0  100  50  150  0  100")
print()
print("  Advantages:")
print("    • Explicit overlap information")
print("    • Detailed alignment data")
print("    • String graph representation")
print()
print("=" * 80)
print()
print("GraphBin handles both formats transparently!")
print("Parser automatically detects format and processes accordingly.")
print()
print("=" * 80)
```

## 10. Result Preparation Utility

The `prep_result` module helps prepare and validate outputs.

```python
print("RESULT PREPARATION UTILITIES")
print("=" * 80)
print()
print("The prep_result module provides utilities for:")
print()
print("1. VALIDATING OUTPUT FILES")
print("   • Check CSV format correctness")
print("   • Verify FASTA file integrity")
print("   • Ensure all contigs are accounted for")
print()
print("2. CONVERTING FORMATS")
print("   • Convert between different binning formats")
print("   • Standardize contig names")
print("   • Merge multiple binning results")
print()
print("3. GENERATING STATISTICS")
print("   • Count contigs per bin")
print("   • Calculate total sequence lengths")
print("   • Identify unbinned contigs")
print()
print("=" * 80)

# Demonstrate validation
print("\nExample: Validate binning results")
print("-" * 80)

def validate_binning_result(binning_file, contigs_file):
    """Simple validation of binning results."""
    
    # Read binning
    binned_contigs = set()
    bins_found = set()
    
    with open(binning_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 2:
                binned_contigs.add(row[0])
                bins_found.add(row[1])
    
    # Read contigs
    all_contigs = set()
    with open(contigs_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                contig_name = line[1:].split()[0]
                all_contigs.add(contig_name)
    
    # Validate
    print(f"Total contigs: {len(all_contigs)}")
    print(f"Binned contigs: {len(binned_contigs)}")
    print(f"Unbinned contigs: {len(all_contigs - binned_contigs)}")
    print(f"Number of bins: {len(bins_found)}")
    print(f"Bins: {sorted(bins_found)}")
    
    # Check for contigs in binning but not in FASTA
    extra = binned_contigs - all_contigs
    if extra:
        print(f"\n⚠️  Warning: {len(extra)} contigs in binning but not in FASTA")
    else:
        print(f"\n✓ All binned contigs found in FASTA file")
    
    return len(all_contigs), len(binned_contigs), len(bins_found)

# Run validation
contigs_file = test_data / 'contigs.fasta'
validate_binning_result(binned_file, contigs_file)

print("\n" + "=" * 80)
```

## 11. Working with BidirectionalMap

GraphBin uses a bidirectional map for efficient node-contig lookups.

```python
from graphbin.bidirectionalmap.bidirectionalmap import BidirectionalMap

print("BIDIRECTIONAL MAP DATA STRUCTURE")
print("=" * 80)
print()
print("The BidirectionalMap allows O(1) lookups in both directions:")
print()

# Create example bidirectional map
bimap = BidirectionalMap()

# Add mappings
bimap['NODE_1'] = 0
bimap['NODE_2'] = 1
bimap['NODE_3'] = 2

print("Forward lookups (name → id):")
print(f"  bimap['NODE_1'] = {bimap['NODE_1']}")
print(f"  bimap['NODE_2'] = {bimap['NODE_2']}")
print()

print("Reverse lookups (id → name):")
print(f"  bimap.inverse[0] = {bimap.inverse[0]}")
print(f"  bimap.inverse[1] = {bimap.inverse[1]}")
print()

print("Use case in GraphBin:")
print("  • Map contig names to graph node indices")
print("  • Quickly find contig name from node ID")
print("  • Efficiently process large assemblies")
print()

# Try to add duplicate (will raise error)
print("Uniqueness enforcement:")
try:
    bimap['NODE_4'] = 0  # 0 is already mapped to NODE_1
except Exception as e:
    print(f"  ✓ Prevented duplicate: {type(e).__name__}")

print("\n" + "=" * 80)
```

## 12. Best Practices for Using GraphBin Utilities

```python
print("BEST PRACTICES FOR GRAPHBIN UTILITIES")
print("=" * 80)
print()
print("1. ALWAYS VALIDATE INPUT FILES")
print("   • Check file formats before running GraphBin")
print("   • Ensure graph file matches assembler type")
print("   • Verify contig names are consistent")
print()
print("2. USE APPROPRIATE PARSER")
print("   • SPAdes → spades_parser")
print("   • MEGAHIT → megahit_parser")
print("   • Don't mix parsers and assemblers!")
print()
print("3. CHECK OUTPUT COMPLETENESS")
print("   • Verify all expected files are created")
print("   • Check for empty or corrupted FASTA files")
print("   • Validate bin assignments make sense")
print()
print("4. VISUALIZE RESULTS")
print("   • Use provided visualization scripts")
print("   • Look for suspicious binning patterns")
print("   • Check if graph structure supports binning")
print()
print("5. COMPARE WITH INITIAL BINNING")
print("   • Count improvements (new bins, corrections)")
print("   • Identify which contigs were affected")
print("   • Validate with external tools (CheckM)")
print()
print("6. HANDLE EDGE CASES")
print("   • Empty bins (remove or keep?)")
print("   • Very short contigs (may remain unbinned)")
print("   • Disconnected graph components")
print()
print("=" * 80)
print()
print("COMMON ISSUES AND SOLUTIONS:")
print("=" * 80)
print()
print("Issue: Parser fails on graph file")
print("  → Check assembler type matches parser")
print("  → Verify graph file format is correct")
print()
print("Issue: Contig names don't match")
print("  → Ensure same contig file used for binning and GraphBin")
print("  → Check for extra spaces or special characters")
print()
print("Issue: No improvement in binning")
print("  → Graph may be poorly connected")
print("  → Initial binning may be already optimal")
print("  → Try different initial binning tool")
print()
print("Issue: Output files are empty")
print("  → Check for errors in log file")
print("  → Verify input files are not corrupted")
print("  → Ensure sufficient disk space")
print()
print("=" * 80)
```

## Summary

In this notebook, we:

1. ✅ Explored all parser modules
2. ✅ Demonstrated graph parsing for different formats
3. ✅ Used binning result parsers
4. ✅ Learned about format conversion utilities
5. ✅ Created custom visualizations
6. ✅ Compared GFA vs ASQG formats
7. ✅ Used result validation utilities
8. ✅ Understood bidirectional map structure
9. ✅ Reviewed best practices

---

### Complete Tutorial Series:

- **Notebook 1**: ✅ Introduction and Setup
- **Notebook 2**: ✅ Core Functionality
- **Notebook 3**: ✅ Multiple Assemblers
- **Notebook 4**: ✅ Advanced Features
- **Notebook 5**: ✅ Parsers and Utilities (You are here!)

---

### Key Takeaways:

- GraphBin provides specialized parsers for 6 assemblers
- All parsers follow consistent interface patterns
- Utilities handle format conversion and validation
- Visualization tools help interpret results
- BidirectionalMap enables efficient lookups
- Always validate inputs and outputs
- Understanding formats helps debug issues

---

## Further Resources

- **GraphBin Documentation**: https://graphbin.readthedocs.io/
- **GraphBin Paper**: https://doi.org/10.1093/bioinformatics/btaa180
- **GitHub Repository**: https://github.com/metagentools/GraphBin
- **Report Issues**: https://github.com/metagentools/GraphBin/issues

Happy binning! 🧬
