# GraphBin Tutorial 3: Multiple Assemblers Support

## Overview

GraphBin supports 6 different assembly tools:
- **SPAdes/metaSPAdes** - Short-read assembler (de Bruijn graph)
- **MEGAHIT** - Memory-efficient short-read assembler
- **SGA** - String Graph Assembler
- **Flye** - Long-read assembler (PacBio/Nanopore)
- **Canu** - Long-read assembler
- **Miniasm** - Ultra-fast long-read assembler

This notebook demonstrates how to use GraphBin with different assemblers and compares their characteristics.

---

## 1. Setup and Imports

```python
import sys
import os
from pathlib import Path
import csv
from collections import Counter

# Add GraphBin to path
repo_root = Path('..').resolve()
src_path = repo_root / 'src'
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

# Import all assembler-specific modules
from graphbin import (
    graphbin_SPAdes,
    graphbin_SGA,
    graphbin_MEGAHIT,
    graphbin_Flye,
    graphbin_Canu,
    graphbin_Miniasm
)
from graphbin.cli import ArgsObj

print("✓ Imports successful")
print(f"Repository root: {repo_root}")
```

## 2. Assembler Characteristics Comparison

```python
import pandas as pd

# Assembler comparison data
assemblers_info = {
    'Assembler': ['SPAdes', 'MEGAHIT', 'SGA', 'Flye', 'Canu', 'Miniasm'],
    'Type': ['Short-read', 'Short-read', 'Short-read', 'Long-read', 'Long-read', 'Long-read'],
    'Graph Format': ['GFA', 'GFA', 'ASQG', 'GFA', 'GFA', 'GFA'],
    'Paths File': ['Yes', 'No', 'No', 'Yes (assembly.info)', 'No', 'No'],
    'Algorithm': ['de Bruijn', 'de Bruijn', 'String graph', 'OLC', 'OLC', 'OLC'],
    'Memory': ['High', 'Low', 'Medium', 'Medium', 'High', 'Low'],
    'Speed': ['Medium', 'Fast', 'Slow', 'Medium', 'Slow', 'Very Fast']
}

df = pd.DataFrame(assemblers_info)
print("ASSEMBLER COMPARISON")
print("=" * 80)
print(df.to_string(index=False))
print()
print("Notes:")
print("  OLC = Overlap-Layout-Consensus")
print("  GFA = Graphical Fragment Assembly format")
print("  ASQG = Assembly String Graph format")
```

## 3. Explore Available Test Datasets

```python
test_data_root = repo_root / 'tests' / 'data'

print("AVAILABLE TEST DATASETS")
print("=" * 80)

# Map dataset directories to assemblers
dataset_map = {
    'ESC_metaSPAdes': 'SPAdes',
    'ESC_MEGAHIT': 'MEGAHIT',
    'ESC_SGA': 'SGA',
    '1Y3B_Flye': 'Flye',
    '1Y3B_Canu': 'Canu',
    '1Y3B_Miniasm': 'Miniasm'
}

available_datasets = {}
for dataset_dir, assembler in dataset_map.items():
    dataset_path = test_data_root / dataset_dir
    if dataset_path.exists():
        files = list(dataset_path.glob('*'))
        available_datasets[assembler] = {
            'path': dataset_path,
            'files': [f.name for f in files if f.is_file()]
        }
        print(f"\n✓ {assembler} ({dataset_dir})")
        print(f"  Files:")
        for fname in sorted(available_datasets[assembler]['files']):
            fpath = dataset_path / fname
            size = fpath.stat().st_size
            size_str = f"{size:,} bytes" if size < 1024*1024 else f"{size/(1024*1024):.2f} MB"
            print(f"    • {fname} ({size_str})")
    else:
        print(f"\n✗ {assembler} ({dataset_dir}) - NOT FOUND")

print(f"\n\nTotal datasets available: {len(available_datasets)}")
```

## 4. Function to Run GraphBin with Any Assembler

Let's create a unified function to run GraphBin with any assembler.

```python
import logging

def run_graphbin_for_assembler(assembler, dataset_path, output_subdir):
    """
    Run GraphBin for a specific assembler.
    
    Args:
        assembler: Assembler name (spades, megahit, sga, flye, canu, miniasm)
        dataset_path: Path to the test dataset
        output_subdir: Subdirectory name for output
    
    Returns:
        Dictionary with results or None if failed
    """
    
    # Define expected file patterns for each assembler
    file_patterns = {
        'spades': {
            'graph': 'assembly_graph_with_scaffolds.gfa',
            'contigs': 'contigs.fasta',
            'paths': 'contigs.paths',
            'binned': 'initial_binning_res.csv'
        },
        'megahit': {
            'graph': 'final.gfa',
            'contigs': 'final.contigs.fa',
            'paths': None,
            'binned': 'initial_binning_res.csv'
        },
        'sga': {
            'graph': 'final.asqg',
            'contigs': 'final.fa',
            'paths': None,
            'binned': 'initial_binning_res.csv'
        },
        'flye': {
            'graph': 'assembly_graph.gfa',
            'contigs': 'assembly.fasta',
            'paths': 'assembly.info',
            'binned': 'initial_binning_res.csv'
        },
        'canu': {
            'graph': 'assembly.gfa',
            'contigs': 'assembly.fasta',
            'paths': None,
            'binned': 'initial_binning_res.csv'
        },
        'miniasm': {
            'graph': 'assembly.gfa',
            'contigs': 'assembly.fasta',
            'paths': None,
            'binned': 'initial_binning_res.csv'
        }
    }
    
    if assembler.lower() not in file_patterns:
        print(f"❌ Unknown assembler: {assembler}")
        return None
    
    patterns = file_patterns[assembler.lower()]
    
    # Build file paths
    graph_file = dataset_path / patterns['graph']
    contigs_file = dataset_path / patterns['contigs']
    paths_file = dataset_path / patterns['paths'] if patterns['paths'] else None
    binned_file = dataset_path / patterns['binned']
    
    # Check files exist
    if not graph_file.exists():
        print(f"❌ Graph file not found: {graph_file}")
        return None
    if not contigs_file.exists():
        print(f"❌ Contigs file not found: {contigs_file}")
        return None
    if not binned_file.exists():
        print(f"❌ Binned file not found: {binned_file}")
        return None
    
    # Create output directory
    output_dir = repo_root / 'notebooks' / 'output' / output_subdir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create arguments
    args = ArgsObj(
        assembler=assembler.lower(),
        graph=str(graph_file),
        contigs=str(contigs_file),
        paths=str(paths_file) if paths_file and paths_file.exists() else None,
        binned=str(binned_file),
        output=str(output_dir),
        prefix='graphbin',
        max_iteration=100,
        diff_threshold=0.1,
        delimiter=','
    )
    
    print(f"\n{'='*80}")
    print(f"Running GraphBin with {assembler.upper()}")
    print(f"{'='*80}")
    print(f"Graph: {graph_file.name}")
    print(f"Contigs: {contigs_file.name}")
    print(f"Paths: {paths_file.name if paths_file and paths_file.exists() else 'N/A'}")
    print(f"Output: {output_dir}")
    
    # Select correct module
    modules = {
        'spades': graphbin_SPAdes,
        'megahit': graphbin_MEGAHIT,
        'sga': graphbin_SGA,
        'flye': graphbin_Flye,
        'canu': graphbin_Canu,
        'miniasm': graphbin_Miniasm
    }
    
    module = modules[assembler.lower()]
    
    try:
        # Run GraphBin
        module.main(args)
        print(f"\n✓ GraphBin completed successfully for {assembler}!")
        
        # Return results
        return {
            'assembler': assembler,
            'output_dir': output_dir,
            'success': True
        }
    except Exception as e:
        print(f"\n❌ Error running GraphBin for {assembler}: {e}")
        import traceback
        traceback.print_exc()
        return None

print("✓ Function defined: run_graphbin_for_assembler()")
```

## 5. Run GraphBin with SPAdes

```python
if 'SPAdes' in available_datasets:
    spades_result = run_graphbin_for_assembler(
        'spades',
        available_datasets['SPAdes']['path'],
        'spades_run'
    )
else:
    print("⚠️  SPAdes dataset not available")
    spades_result = None
```

## 6. Run GraphBin with MEGAHIT

```python
if 'MEGAHIT' in available_datasets:
    megahit_result = run_graphbin_for_assembler(
        'megahit',
        available_datasets['MEGAHIT']['path'],
        'megahit_run'
    )
else:
    print("⚠️  MEGAHIT dataset not available")
    megahit_result = None
```

## 7. Run GraphBin with SGA

```python
if 'SGA' in available_datasets:
    sga_result = run_graphbin_for_assembler(
        'sga',
        available_datasets['SGA']['path'],
        'sga_run'
    )
else:
    print("⚠️  SGA dataset not available")
    sga_result = None
```

## 8. Run GraphBin with Flye (Long-read)

```python
if 'Flye' in available_datasets:
    flye_result = run_graphbin_for_assembler(
        'flye',
        available_datasets['Flye']['path'],
        'flye_run'
    )
else:
    print("⚠️  Flye dataset not available")
    flye_result = None
```

## 9. Run GraphBin with Canu (Long-read)

```python
if 'Canu' in available_datasets:
    canu_result = run_graphbin_for_assembler(
        'canu',
        available_datasets['Canu']['path'],
        'canu_run'
    )
else:
    print("⚠️  Canu dataset not available")
    canu_result = None
```

## 10. Run GraphBin with Miniasm (Long-read)

```python
if 'Miniasm' in available_datasets:
    miniasm_result = run_graphbin_for_assembler(
        'miniasm',
        available_datasets['Miniasm']['path'],
        'miniasm_run'
    )
else:
    print("⚠️  Miniasm dataset not available")
    miniasm_result = None
```

## 11. Compare Results Across Assemblers

Let's analyze and compare the refinement results from different assemblers.

```python
def analyze_results(result_info, dataset_path):
    """Analyze GraphBin results for an assembler."""
    
    if not result_info or not result_info.get('success'):
        return None
    
    output_dir = result_info['output_dir']
    assembler = result_info['assembler']
    
    # Read initial binning
    initial_file = dataset_path / 'initial_binning_res.csv'
    initial_bins = {}
    if initial_file.exists():
        with open(initial_file, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                if len(row) >= 2:
                    initial_bins[row[0]] = row[1]
    
    # Read refined binning
    refined_file = output_dir / 'graphbin_output.csv'
    refined_bins = {}
    if refined_file.exists():
        with open(refined_file, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            for row in reader:
                if len(row) >= 2:
                    refined_bins[row[0]] = row[1]
    
    # Get contig count from contigs file
    contigs_patterns = ['contigs.fasta', 'final.contigs.fa', 'final.fa', 'assembly.fasta']
    total_contigs = 0
    for pattern in contigs_patterns:
        contigs_file = dataset_path / pattern
        if contigs_file.exists():
            with open(contigs_file, 'r') as f:
                total_contigs = sum(1 for line in f if line.startswith('>'))
            break
    
    # Calculate statistics
    initial_count = len(initial_bins)
    refined_count = len(refined_bins)
    recovered = refined_count - initial_count
    
    # Count changes
    changed = 0
    for contig in initial_bins:
        if contig in refined_bins and initial_bins[contig] != refined_bins[contig]:
            changed += 1
    
    return {
        'Assembler': assembler,
        'Total Contigs': total_contigs,
        'Initial Binned': initial_count,
        'Refined Binned': refined_count,
        'Recovered': recovered,
        'Changed': changed,
        'Initial Coverage %': f"{100*initial_count/total_contigs:.1f}" if total_contigs > 0 else 'N/A',
        'Refined Coverage %': f"{100*refined_count/total_contigs:.1f}" if total_contigs > 0 else 'N/A'
    }

# Collect all results
all_results = []
results_map = {
    'SPAdes': (spades_result, 'SPAdes'),
    'MEGAHIT': (megahit_result, 'MEGAHIT'),
    'SGA': (sga_result, 'SGA'),
    'Flye': (flye_result, 'Flye'),
    'Canu': (canu_result, 'Canu'),
    'Miniasm': (miniasm_result, 'Miniasm')
}

for assembler_name, (result, key) in results_map.items():
    if result and key in available_datasets:
        stats = analyze_results(result, available_datasets[key]['path'])
        if stats:
            all_results.append(stats)

if all_results:
    comparison_df = pd.DataFrame(all_results)
    print("\n" + "="*100)
    print("GRAPHBIN RESULTS COMPARISON ACROSS ASSEMBLERS")
    print("="*100)
    print(comparison_df.to_string(index=False))
    print("\n" + "="*100)
else:
    print("\n⚠️  No results to compare")
```

## 12. Visualize Assembler Comparison

```python
import matplotlib.pyplot as plt
import numpy as np

if all_results:
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    
    assemblers = [r['Assembler'] for r in all_results]
    initial_binned = [r['Initial Binned'] for r in all_results]
    refined_binned = [r['Refined Binned'] for r in all_results]
    recovered = [r['Recovered'] for r in all_results]
    changed = [r['Changed'] for r in all_results]
    
    x = np.arange(len(assemblers))
    width = 0.35
    
    # Plot 1: Initial vs Refined binned contigs
    axes[0, 0].bar(x - width/2, initial_binned, width, label='Initial', color='steelblue', alpha=0.7)
    axes[0, 0].bar(x + width/2, refined_binned, width, label='Refined', color='forestgreen', alpha=0.7)
    axes[0, 0].set_xlabel('Assembler')
    axes[0, 0].set_ylabel('Binned Contigs')
    axes[0, 0].set_title('Initial vs Refined Binning')
    axes[0, 0].set_xticks(x)
    axes[0, 0].set_xticklabels(assemblers, rotation=45)
    axes[0, 0].legend()
    axes[0, 0].grid(axis='y', alpha=0.3)
    
    # Plot 2: Recovered contigs
    colors = ['forestgreen' if r > 0 else 'gray' for r in recovered]
    axes[0, 1].bar(assemblers, recovered, color=colors, alpha=0.7)
    axes[0, 1].set_xlabel('Assembler')
    axes[0, 1].set_ylabel('Recovered Contigs')
    axes[0, 1].set_title('Newly Binned Contigs')
    axes[0, 1].tick_params(axis='x', rotation=45)
    axes[0, 1].grid(axis='y', alpha=0.3)
    
    # Plot 3: Changed bin assignments
    axes[1, 0].bar(assemblers, changed, color='orange', alpha=0.7)
    axes[1, 0].set_xlabel('Assembler')
    axes[1, 0].set_ylabel('Changed Assignments')
    axes[1, 0].set_title('Corrected Mis-binned Contigs')
    axes[1, 0].tick_params(axis='x', rotation=45)
    axes[1, 0].grid(axis='y', alpha=0.3)
    
    # Plot 4: Coverage percentages
    initial_pct = [float(r['Initial Coverage %']) for r in all_results]
    refined_pct = [float(r['Refined Coverage %']) for r in all_results]
    
    axes[1, 1].bar(x - width/2, initial_pct, width, label='Initial', color='steelblue', alpha=0.7)
    axes[1, 1].bar(x + width/2, refined_pct, width, label='Refined', color='forestgreen', alpha=0.7)
    axes[1, 1].set_xlabel('Assembler')
    axes[1, 1].set_ylabel('Coverage %')
    axes[1, 1].set_title('Binning Coverage Percentage')
    axes[1, 1].set_xticks(x)
    axes[1, 1].set_xticklabels(assemblers, rotation=45)
    axes[1, 1].legend()
    axes[1, 1].grid(axis='y', alpha=0.3)
    axes[1, 1].set_ylim([0, 100])
    
    plt.tight_layout()
    
    # Save plot
    output_plot = repo_root / 'notebooks' / 'output' / 'assembler_comparison.png'
    plt.savefig(output_plot, dpi=150, bbox_inches='tight')
    print(f"\n✓ Comparison plot saved to: {output_plot}")
    plt.show()
else:
    print("\n⚠️  No results available for visualization")
```

## 13. Key Differences Between Assemblers

Let's understand the key differences in how GraphBin handles each assembler.

```python
print("KEY DIFFERENCES IN GRAPHBIN PROCESSING")
print("=" * 80)
print()
print("📊 SHORT-READ ASSEMBLERS:")
print()
print("SPAdes/metaSPAdes:")
print("  • Uses contigs.paths file to map contigs to graph paths")
print("  • GFA format with multiple node types (edges, scaffolds)")
print("  • Highest complexity but most information")
print("  • Best for complex metagenomes")
print()
print("MEGAHIT:")
print("  • No paths file - direct contig-to-node mapping")
print("  • Simplified GFA format")
print("  • Memory efficient, faster processing")
print("  • Good for large datasets")
print()
print("SGA (String Graph Assembler):")
print("  • Uses ASQG format instead of GFA")
print("  • String graph representation (overlaps, not k-mers)")
print("  • Better for repetitive sequences")
print("  • Different graph structure than de Bruijn graphs")
print()
print("=" * 80)
print("📊 LONG-READ ASSEMBLERS:")
print()
print("Flye:")
print("  • Uses assembly.info file (like SPAdes paths)")
print("  • Repeat graph structure")
print("  • Designed for PacBio/Nanopore reads")
print("  • Handles long read errors well")
print()
print("Canu:")
print("  • No paths/info file")
print("  • Overlap-Layout-Consensus approach")
print("  • Most accurate but slowest")
print("  • Large genome support")
print()
print("Miniasm:")
print("  • Simplest, fastest long-read assembler")
print("  • No consensus/polishing step")
print("  • Higher error rate but very fast")
print("  • Good for quick assembly drafts")
print()
print("=" * 80)
print("💡 KEY INSIGHT:")
print()
print("GraphBin's label propagation algorithm works the SAME way for all")
print("assemblers, but the input parsing differs based on graph format!")
print()
print("The quality of refinement depends on:")
print("  1. Assembly graph connectivity (more edges = better)")
print("  2. Initial binning quality")
print("  3. Contig length distribution")
print("  4. Complexity of the microbial community")
print("=" * 80)
```

## 14. Summary and Recommendations

```python
print("\n" + "=" * 80)
print("ASSEMBLER SELECTION GUIDE")
print("=" * 80)
print()
print("Choose your assembler based on:")
print()
print("🧬 DATA TYPE:")
print("  Short reads (Illumina)     → SPAdes, MEGAHIT, or SGA")
print("  Long reads (PacBio/Oxford) → Flye, Canu, or Miniasm")
print()
print("💾 MEMORY CONSTRAINTS:")
print("  Limited memory  → MEGAHIT (short) or Miniasm (long)")
print("  Plenty of RAM   → SPAdes (short) or Canu (long)")
print()
print("⚡ SPEED REQUIREMENTS:")
print("  Quick results   → MEGAHIT (short) or Miniasm (long)")
print("  Quality focus   → SPAdes (short) or Flye/Canu (long)")
print()
print("🎯 ACCURACY NEEDS:")
print("  Highest quality → SPAdes or Canu")
print("  Quick draft     → MEGAHIT or Miniasm")
print()
print("=" * 80)
print("GRAPHBIN PERFORMANCE WITH DIFFERENT ASSEMBLERS:")
print("=" * 80)
print()
print("Best results typically with:")
print("  ✓ Well-connected assembly graphs")
print("  ✓ Higher-quality initial assemblies (SPAdes, Flye, Canu)")
print("  ✓ Moderate to high coverage datasets")
print()
print("GraphBin adds most value when:")
print("  • Initial binning has many short unbinned contigs")
print("  • Assembly graph has good connectivity")
print("  • Community is moderately complex (3-20 species)")
print()
print("=" * 80)
```

## Next Steps

In this notebook, we:

1. ✅ Compared 6 different assemblers supported by GraphBin
2. ✅ Ran GraphBin with multiple assemblers
3. ✅ Analyzed and compared refinement results
4. ✅ Visualized performance differences
5. ✅ Understood key differences in processing
6. ✅ Learned when to use each assembler

### Continue to:

- **Notebook 4**: Deep dive into label propagation algorithm and graph analysis
- **Notebook 5**: Parser functions, visualization, and utilities

---

### Key Takeaways:

- GraphBin supports both short-read and long-read assemblers
- The core algorithm is the same; only parsing differs
- Better assembly quality → better refinement results
- Choose assembler based on data type, resources, and goals
- GraphBin can work with any assembler that produces a graph
