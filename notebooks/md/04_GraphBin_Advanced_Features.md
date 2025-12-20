# GraphBin Tutorial 4: Advanced Features - Label Propagation and Graph Analysis

## Overview

This notebook provides a deep dive into GraphBin's core algorithm and advanced features:
- Label Propagation algorithm implementation
- Graph analysis and connectivity
- Algorithm visualization
- Parameter tuning
- Performance optimization

---

## 1. Setup and Imports

```python
import sys
import os
from pathlib import Path
import csv
from collections import Counter, defaultdict
import numpy as np
import matplotlib.pyplot as plt

# Add GraphBin to path
repo_root = Path('..').resolve()
src_path = repo_root / 'src'
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

# Import GraphBin modules
from graphbin.labelpropagation.labelprop import LabelProp, Edge
from graphbin.parsers.spades_parser import parse_graph, get_initial_binning_result
from graphbin import graphbin_Func
import igraph as ig

print("✓ Imports successful")
print(f"Repository root: {repo_root}")
```

## 2. Understanding Label Propagation Algorithm

Label propagation is a semi-supervised learning algorithm that spreads labels through a graph based on connectivity.

```python
print("LABEL PROPAGATION ALGORITHM - CONCEPTUAL OVERVIEW")
print("=" * 80)
print()
print("Think of it like spreading paint on a connected network:")
print()
print("1. INITIAL STATE")
print("   • Some nodes are colored (labeled/binned contigs)")
print("   • Other nodes are white (unlabeled/unbinned contigs)")
print("   • Nodes are connected by edges (assembly graph links)")
print()
print("2. PROPAGATION STEP (repeated iteratively)")
print("   For each unlabeled node:")
print("     a) Look at all its neighbors")
print("     b) Calculate weighted average of neighbor labels")
print("     c) Update this node's label probabilities")
print()
print("3. CONVERGENCE")
print("   • Repeat until changes become very small")
print("   • Labeled nodes keep their original labels")
print("   • Unlabeled nodes adopt most likely label from neighbors")
print()
print("4. FINAL ASSIGNMENT")
print("   • Each unlabeled node gets the label with highest probability")
print()
print("=" * 80)
print("KEY EQUATIONS:")
print("=" * 80)
print()
print("For unlabeled vertex v and label l:")
print()
print("  F(v, l) = Σ [F(u, l) × w(u,v) / deg(v)]")
print("            u∈neighbors(v)")
print()
print("Where:")
print("  • F(v, l) = probability that vertex v has label l")
print("  • w(u,v) = edge weight between vertices u and v")
print("  • deg(v) = degree (sum of edge weights) of vertex v")
print()
print("=" * 80)
```

## 3. Simple Label Propagation Example

Let's create a simple example to see label propagation in action.

```python
# Create a simple example graph
# Graph structure:
#   [1] --- [2] --- [3]
#    |       |       |
#   [4] --- [5] --- [6]
#
# Initial labels: 1, 3 = Label A; 4, 6 = Label B; 2, 5 = unlabeled

example_data = [
    [1, 'A', [[2, 1.0], [4, 1.0]]],           # Node 1: Label A, connected to 2, 4
    [2, 0, [[1, 1.0], [3, 1.0], [5, 1.0]]],   # Node 2: Unlabeled, connected to 1, 3, 5
    [3, 'A', [[2, 1.0], [6, 1.0]]],           # Node 3: Label A, connected to 2, 6
    [4, 'B', [[1, 1.0], [5, 1.0]]],           # Node 4: Label B, connected to 1, 5
    [5, 0, [[2, 1.0], [4, 1.0], [6, 1.0]]],   # Node 5: Unlabeled, connected to 2, 4, 6
    [6, 'B', [[3, 1.0], [5, 1.0]]]            # Node 6: Label B, connected to 3, 5
]

# Create and run label propagation
lp = LabelProp()
lp.load_data_from_mem(example_data)

print("SIMPLE LABEL PROPAGATION EXAMPLE")
print("=" * 80)
print()
print("Initial state:")
print("  Node 1: Label A (fixed)")
print("  Node 2: Unlabeled")
print("  Node 3: Label A (fixed)")
print("  Node 4: Label B (fixed)")
print("  Node 5: Unlabeled")
print("  Node 6: Label B (fixed)")
print()
print("Running label propagation...")
print()

# Run with detailed output
result = lp.run(eps=1.0, max_iter=10, show_log=True, clean_result=False)

print("\nFinal labels:")
for item in result:
    node_id = item[0]
    final_label = item[1]
    label_probs = item[2:]
    print(f"  Node {node_id}: {final_label} (probabilities: {label_probs})")

print("\n" + "=" * 80)
print("INTERPRETATION:")
print("=" * 80)
print("Node 2 is connected to nodes 1 (A), 3 (A), and 5 (unlabeled)")
print("→ Strong influence from Label A neighbors → likely Label A")
print()
print("Node 5 is connected to nodes 2 (unlabeled), 4 (B), and 6 (B)")
print("→ Strong influence from Label B neighbors → likely Label B")
print("=" * 80)
```

## 4. Visualize Simple Example

```python
# Create visualization of the example graph
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Graph structure
G = ig.Graph()
G.add_vertices(6)
G.add_edges([(0,1), (0,3), (1,2), (1,4), (2,5), (3,4), (4,5)])
G.vs['label'] = [1, 2, 3, 4, 5, 6]

# Initial state
initial_colors = ['red', 'white', 'red', 'blue', 'white', 'blue']
G.vs['color'] = initial_colors
layout = G.layout('fr')

ig.plot(G, target=ax1, layout=layout, 
        vertex_size=50, 
        vertex_label=G.vs['label'],
        edge_width=2,
        vertex_label_size=12)
ax1.set_title('Initial State\n(Red=A, Blue=B, White=Unlabeled)', fontsize=12)
ax1.axis('off')

# After label propagation
# Node 2 likely becomes A, Node 5 likely becomes B
final_colors = ['red', 'red', 'red', 'blue', 'blue', 'blue']
G.vs['color'] = final_colors

ig.plot(G, target=ax2, layout=layout,
        vertex_size=50,
        vertex_label=G.vs['label'],
        edge_width=2,
        vertex_label_size=12)
ax2.set_title('After Label Propagation\n(All nodes labeled)', fontsize=12)
ax2.axis('off')

plt.tight_layout()
plt.savefig(repo_root / 'notebooks' / 'output' / 'label_propagation_example.png', 
            dpi=150, bbox_inches='tight')
print(f"\n✓ Visualization saved")
plt.show()
```

## 5. Load Real Assembly Graph

Now let's analyze a real assembly graph from the test data.

```python
# Load test data
test_data = repo_root / 'tests' / 'data' / 'ESC_metaSPAdes'
graph_file = test_data / 'assembly_graph_with_scaffolds.gfa'
contigs_file = test_data / 'contigs.fasta'
paths_file = test_data / 'contigs.paths'
binned_file = test_data / 'initial_binning_res.csv'

print("Loading assembly graph...")

# Parse graph using GraphBin parser
from graphbin.parsers.spades_parser import parse_graph

# Read contig paths
contig_paths = {}
with open(paths_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line:
            parts = line.split()
            if len(parts) >= 2:
                contig_name = parts[0]
                path = parts[1]
                contig_paths[contig_name] = path

print(f"✓ Loaded {len(contig_paths)} contig paths")

# Parse the assembly graph
assembly_graph, contig_names, node_count = parse_graph(
    str(graph_file), 
    contig_paths
)

print(f"✓ Parsed assembly graph")
print(f"  Vertices: {assembly_graph.vcount()}")
print(f"  Edges: {assembly_graph.ecount()}")
print(f"  Contigs mapped: {len(contig_names)}")
```

## 6. Analyze Graph Properties

```python
print("ASSEMBLY GRAPH ANALYSIS")
print("=" * 80)

# Basic statistics
degrees = assembly_graph.degree()
print(f"\nDegree statistics:")
print(f"  Min degree: {min(degrees)}")
print(f"  Max degree: {max(degrees)}")
print(f"  Mean degree: {np.mean(degrees):.2f}")
print(f"  Median degree: {np.median(degrees):.0f}")

# Connected components
components = assembly_graph.connected_components()
print(f"\nConnected components: {len(components)}")
component_sizes = [len(c) for c in components]
print(f"  Largest component: {max(component_sizes)} nodes")
print(f"  Smallest component: {min(component_sizes)} nodes")
print(f"  Mean component size: {np.mean(component_sizes):.1f} nodes")

# Clustering coefficient
try:
    clustering = assembly_graph.transitivity_undirected()
    print(f"\nGlobal clustering coefficient: {clustering:.4f}")
except:
    print(f"\nClustering coefficient: N/A")

# Diameter (for largest component if graph is disconnected)
if len(components) > 1:
    largest_component = components[component_sizes.index(max(component_sizes))]
    subgraph = assembly_graph.subgraph(largest_component)
    try:
        diameter = subgraph.diameter()
        print(f"\nDiameter (largest component): {diameter}")
    except:
        print(f"\nDiameter: N/A (graph too large)")
else:
    try:
        diameter = assembly_graph.diameter()
        print(f"\nDiameter: {diameter}")
    except:
        print(f"\nDiameter: N/A (graph too large)")

print("\n" + "=" * 80)
```

## 7. Visualize Degree Distribution

```python
# Plot degree distribution
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Histogram
ax1.hist(degrees, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
ax1.set_xlabel('Degree', fontsize=12)
ax1.set_ylabel('Frequency', fontsize=12)
ax1.set_title('Degree Distribution (Linear)', fontsize=14)
ax1.grid(alpha=0.3)

# Log-log plot
degree_counts = Counter(degrees)
degrees_sorted = sorted(degree_counts.keys())
counts = [degree_counts[d] for d in degrees_sorted]

ax2.loglog(degrees_sorted, counts, 'o', color='forestgreen', alpha=0.6, markersize=8)
ax2.set_xlabel('Degree (log scale)', fontsize=12)
ax2.set_ylabel('Count (log scale)', fontsize=12)
ax2.set_title('Degree Distribution (Log-Log)', fontsize=14)
ax2.grid(alpha=0.3, which='both')

plt.tight_layout()
plt.savefig(repo_root / 'notebooks' / 'output' / 'degree_distribution.png', 
            dpi=150, bbox_inches='tight')
print("✓ Degree distribution plot saved")
plt.show()

print("\n💡 Interpretation:")
print("  • Higher degree nodes are more influential in label propagation")
print("  • Well-connected contigs propagate labels more effectively")
print("  • Isolated contigs (degree 0-1) are hard to bin")
```

## 8. Analyze Binning Impact on Graph Structure

```python
# Load initial binning
initial_bins = {}
with open(binned_file, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader:
        if len(row) >= 2:
            initial_bins[row[0]] = row[1]

# Map contigs to node indices
contig_to_node = {}
for node_idx, contig_name in enumerate(contig_names):
    contig_to_node[contig_name] = node_idx

# Analyze binned vs unbinned connectivity
binned_nodes = set()
unbinned_nodes = set()

for contig_name in contig_names:
    if contig_name in initial_bins:
        if contig_name in contig_to_node:
            binned_nodes.add(contig_to_node[contig_name])
    else:
        if contig_name in contig_to_node:
            unbinned_nodes.add(contig_to_node[contig_name])

print("BINNING AND GRAPH CONNECTIVITY")
print("=" * 80)
print(f"\nBinned nodes: {len(binned_nodes)}")
print(f"Unbinned nodes: {len(unbinned_nodes)}")

# Calculate average degree for each group
binned_degrees = [degrees[n] for n in binned_nodes]
unbinned_degrees = [degrees[n] for n in unbinned_nodes]

print(f"\nAverage degree:")
print(f"  Binned nodes: {np.mean(binned_degrees):.2f}")
print(f"  Unbinned nodes: {np.mean(unbinned_degrees):.2f}")

# Count edges between groups
binned_to_binned = 0
binned_to_unbinned = 0
unbinned_to_unbinned = 0

for edge in assembly_graph.es:
    src = edge.source
    dst = edge.target
    
    src_binned = src in binned_nodes
    dst_binned = dst in binned_nodes
    
    if src_binned and dst_binned:
        binned_to_binned += 1
    elif (src_binned and not dst_binned) or (not src_binned and dst_binned):
        binned_to_unbinned += 1
    else:
        unbinned_to_unbinned += 1

print(f"\nEdge connectivity:")
print(f"  Binned ↔ Binned: {binned_to_binned} edges")
print(f"  Binned ↔ Unbinned: {binned_to_unbinned} edges ← KEY for refinement!")
print(f"  Unbinned ↔ Unbinned: {unbinned_to_unbinned} edges")

print("\n" + "=" * 80)
print("💡 INSIGHT:")
print(f"  {binned_to_unbinned} edges connect binned to unbinned contigs.")
print(f"  Label propagation uses these connections to assign bins to unbinned contigs!")
print("=" * 80)
```

## 9. Parameter Sensitivity Analysis

GraphBin has two key parameters: `max_iteration` and `diff_threshold`. Let's see how they affect convergence.

```python
print("PARAMETER SENSITIVITY")
print("=" * 80)
print()
print("GraphBin Parameters:")
print()
print("1. max_iteration (default: 100)")
print("   • Maximum number of label propagation iterations")
print("   • Higher = more refinement but slower")
print("   • Typical convergence: 10-50 iterations")
print()
print("2. diff_threshold (default: 0.1)")
print("   • Convergence threshold (percentage of vertex_size)")
print("   • Stop when change < threshold × vertex_size")
print("   • Lower = stricter convergence criteria")
print()
print("3. delimiter (default: ',')")
print("   • CSV delimiter for input/output files")
print()
print("=" * 80)
print("RECOMMENDATIONS:")
print("=" * 80)
print()
print("Standard datasets:")
print("  max_iteration = 100, diff_threshold = 0.1")
print()
print("Large/complex datasets:")
print("  max_iteration = 200, diff_threshold = 0.05")
print()
print("Quick/draft refinement:")
print("  max_iteration = 50, diff_threshold = 0.2")
print()
print("High-quality refinement:")
print("  max_iteration = 500, diff_threshold = 0.01")
print()
print("=" * 80)
```

## 10. Understanding the GraphBin Workflow

Let's trace through the actual GraphBin algorithm step-by-step.

```python
from graphbin import graphbin_Func

print("GRAPHBIN ALGORITHM WORKFLOW")
print("=" * 80)
print()
print("Step 1: LOAD DATA")
print("  • Parse assembly graph (GFA/ASQG format)")
print("  • Read contigs and map to graph nodes")
print("  • Load initial binning results")
print()
print("Step 2: IDENTIFY AMBIGUOUS VERTICES")
print("  Function: graphbin_Func.getClosestLabelledVertices()")
print("  • For each binned contig, check its neighbors")
print("  • If neighbors have different bin labels → AMBIGUOUS")
print("  • Remove ambiguous labels (will be re-propagated)")
print()
print("Step 3: PREPARE LABEL PROPAGATION DATA")
print("  • Convert igraph to label propagation format")
print("  • Each node: [node_id, label, [(neighbor_id, weight), ...]]")
print("  • Labels: bin names for binned, 0 for unbinned")
print()
print("Step 4: RUN LABEL PROPAGATION")
print("  Function: LabelProp.run()")
print("  • Initialize probability distributions")
print("  • Iterate: update probabilities based on neighbors")
print("  • Converge when changes < threshold")
print()
print("Step 5: ASSIGN FINAL LABELS")
print("  • Each unlabeled node gets label with max probability")
print("  • Previously labeled nodes keep original labels")
print()
print("Step 6: WRITE OUTPUTS")
print("  • graphbin_output.csv: refined binning results")
print("  • graphbin_bin_X.fasta: sequences for each bin")
print("  • graphbin.log: detailed log file")
print()
print("=" * 80)

# Show the key functions
print("\nKEY FUNCTIONS:")
print("=" * 80)
print()
print("graphbin_Func.getClosestLabelledVertices(graph, node, binned_contigs)")
print("  • Finds nearest labeled neighbors using BFS")
print("  • Returns list of labeled vertices at closest distance")
print()
print("graphbin_Func.graphbin_main(...)")
print("  • Main refinement algorithm")
print("  • Identifies ambiguous vertices")
print("  • Runs label propagation")
print("  • Returns refined binning")
print()
print("LabelProp.iterate()")
print("  • Single label propagation iteration")
print("  • Updates probability for each unlabeled vertex")
print("  • Returns total change (for convergence check)")
print()
print("=" * 80)
```

## 11. Edge Cases and Limitations

```python
print("EDGE CASES AND LIMITATIONS")
print("=" * 80)
print()
print("⚠️  When GraphBin may struggle:")
print()
print("1. DISCONNECTED GRAPH")
print("   • Isolated contigs have no neighbors")
print("   • Cannot propagate labels to disconnected components")
print("   • Solution: These contigs remain unbinned")
print()
print("2. POOR INITIAL BINNING")
print("   • If initial bins are mostly wrong")
print("   • Label propagation amplifies errors")
print("   • Solution: Use better initial binning tool")
print()
print("3. CHIMERIC CONTIGS")
print("   • Contigs incorrectly joining two organisms")
print("   • May connect unrelated bins in graph")
print("   • Solution: Ambiguous vertex detection helps")
print()
print("4. LOW COVERAGE REGIONS")
print("   • Sparse graph connectivity")
print("   • Few edges to propagate labels")
print("   • Solution: Limited improvement possible")
print()
print("5. HIGHLY SIMILAR ORGANISMS")
print("   • Contigs may be misassembled across strains")
print("   • Graph connects related organisms")
print("   • Solution: May merge similar bins (expected)")
print()
print("=" * 80)
print("✓ BEST PRACTICES:")
print("=" * 80)
print()
print("1. Use high-quality initial binning (MaxBin2, MetaBAT2, CONCOCT)")
print("2. Ensure good assembly quality (SPAdes metagenomic mode)")
print("3. Check assembly graph connectivity before refinement")
print("4. Validate results with CheckM or similar tools")
print("5. Compare multiple binning tools + GraphBin")
print()
print("=" * 80)
```

## 12. Performance Optimization Tips

```python
print("PERFORMANCE OPTIMIZATION")
print("=" * 80)
print()
print("🚀 Speed Improvements:")
print()
print("1. USE MEMORY-EFFICIENT ASSEMBLERS")
print("   • MEGAHIT instead of SPAdes for large datasets")
print("   • Reduces graph parsing time")
print()
print("2. FILTER SHORT CONTIGS")
print("   • Pre-filter contigs < 1000 bp before binning")
print("   • Reduces number of nodes in graph")
print("   • Can recover these later with GraphBin")
print()
print("3. LOWER MAX_ITERATION")
print("   • Most convergence happens in first 20-30 iterations")
print("   • Set max_iteration = 50 for faster results")
print()
print("4. INCREASE DIFF_THRESHOLD")
print("   • Higher threshold = earlier stopping")
print("   • diff_threshold = 0.2 for quick drafts")
print()
print("=" * 80)
print("📊 Typical Performance:")
print("=" * 80)
print()
print("Small dataset (100 contigs):")
print("  • Runtime: < 1 minute")
print("  • Memory: < 100 MB")
print()
print("Medium dataset (1,000 contigs):")
print("  • Runtime: 2-5 minutes")
print("  • Memory: 200-500 MB")
print()
print("Large dataset (10,000 contigs):")
print("  • Runtime: 10-30 minutes")
print("  • Memory: 1-3 GB")
print()
print("Very large dataset (100,000 contigs):")
print("  • Runtime: 1-3 hours")
print("  • Memory: 10-20 GB")
print()
print("=" * 80)
```

## Summary

In this notebook, we:

1. ✅ Explained the label propagation algorithm conceptually
2. ✅ Demonstrated label propagation with a simple example
3. ✅ Analyzed real assembly graph properties
4. ✅ Examined graph connectivity and binning
5. ✅ Discussed parameter tuning
6. ✅ Traced through the GraphBin workflow
7. ✅ Identified edge cases and limitations
8. ✅ Provided optimization tips

### Continue to:

- **Notebook 5**: Parsers, utilities, and visualization tools

---

### Key Takeaways:

- Label propagation spreads labels through graph connectivity
- Well-connected graphs produce better refinement
- Algorithm typically converges in 10-50 iterations
- Parameters can be tuned for speed vs accuracy
- GraphBin works best with good initial binning and assembly
- Understanding graph structure helps predict refinement success
