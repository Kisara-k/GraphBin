# GraphBin Jupyter Notebooks - Quick Reference

## 📚 Complete Tutorial Series

I've created **5 comprehensive Jupyter notebooks** that walk through all GraphBin functionalities:

### 1️⃣ Introduction and Setup (01_GraphBin_Introduction_and_Setup.ipynb)

- ✅ Installation verification
- ✅ Understanding metagenomic binning concepts
- ✅ Exploring test data structure
- ✅ Input file formats (FASTA, GFA, CSV)
- ✅ Module structure overview

### 2️⃣ Core Functionality (02_GraphBin_Core_Functionality.ipynb)

- ✅ Running GraphBin programmatically
- ✅ Analyzing initial binning
- ✅ Assembly graph connectivity
- ✅ Label propagation refinement
- ✅ Comparing results (before/after)
- ✅ Visualizing improvements

### 3️⃣ Multiple Assemblers (03_GraphBin_Multiple_Assemblers.ipynb)

- ✅ SPAdes, MEGAHIT, SGA (short-read)
- ✅ Flye, Canu, Miniasm (long-read)
- ✅ Running with each assembler
- ✅ Comparing results across assemblers
- ✅ Assembler selection guide

### 4️⃣ Advanced Features (04_GraphBin_Advanced_Features.ipynb)

- ✅ Label propagation algorithm explained
- ✅ Simple demonstration example
- ✅ Graph analysis and properties
- ✅ Parameter tuning
- ✅ Algorithm workflow details
- ✅ Edge cases and limitations
- ✅ Performance optimization

### 5️⃣ Parsers and Utilities (05_GraphBin_Parsers_and_Utilities.ipynb)

- ✅ Parser modules for each assembler
- ✅ Graph format parsing (GFA, ASQG)
- ✅ Format conversion utilities
- ✅ Output validation
- ✅ Visualization utilities
- ✅ Best practices

---

## 🎯 Key Features

### Import Functions from src/

All notebooks import GraphBin functions directly instead of rewriting them:

```python
# Import from source
from graphbin import graphbin_SPAdes, graphbin_MEGAHIT
from graphbin.parsers import spades_parser, megahit_parser
from graphbin.labelpropagation.labelprop import LabelProp
from graphbin import graphbin_Func

# Use functions directly
assembly_graph, contig_names, node_count = spades_parser.parse_graph(
    graph_file, contig_paths
)
```

### Comprehensive Coverage

- ✅ All 6 assemblers (SPAdes, MEGAHIT, SGA, Flye, Canu, Miniasm)
- ✅ All parser modules
- ✅ Label propagation algorithm
- ✅ Graph analysis utilities
- ✅ Visualization tools
- ✅ Output processing

### Educational Approach

- 📖 Detailed explanations of concepts
- 💡 Visual demonstrations
- 🔍 Step-by-step walkthroughs
- 📊 Comparative analyses
- ⚠️ Common pitfalls and solutions

---

## 🚀 Quick Start

```bash
# Navigate to notebooks directory
cd notebooks/

# Launch Jupyter
jupyter notebook

# Open 01_GraphBin_Introduction_and_Setup.ipynb
# Execute cells sequentially
```

---

## 📂 What Was Created

```
notebooks/
├── README.md                                    # Detailed guide
├── NOTEBOOK_SUMMARY.md                         # This file
├── 01_GraphBin_Introduction_and_Setup.ipynb    # Intro & setup
├── 02_GraphBin_Core_Functionality.ipynb        # Core refinement
├── 03_GraphBin_Multiple_Assemblers.ipynb       # All assemblers
├── 04_GraphBin_Advanced_Features.ipynb         # Label propagation
├── 05_GraphBin_Parsers_and_Utilities.ipynb     # Parsers & utils
└── output/                                      # Generated outputs
    ├── spades_refinement/                       # Example outputs
    ├── megahit_run/
    └── ...
```

---

## 🎓 Learning Objectives

By completing these notebooks, you will:

1. **Understand** metagenomic binning and refinement
2. **Run** GraphBin with any supported assembler
3. **Analyze** assembly graphs and connectivity
4. **Interpret** label propagation algorithm
5. **Tune** parameters for optimal results
6. **Visualize** binning results
7. **Troubleshoot** common issues
8. **Integrate** GraphBin into workflows

---

## 💻 Code Examples

### Running GraphBin

```python
from graphbin import graphbin_SPAdes
from graphbin.cli import ArgsObj

args = ArgsObj(
    assembler='spades',
    graph='assembly_graph.gfa',
    contigs='contigs.fasta',
    paths='contigs.paths',
    binned='initial_bins.csv',
    output='output_dir/',
    prefix='graphbin',
    max_iteration=100,
    diff_threshold=0.1,
    delimiter=','
)

graphbin_SPAdes.main(args)
```

### Parsing Graphs

```python
from graphbin.parsers import spades_parser

graph, contig_names, node_count = spades_parser.parse_graph(
    graph_file, contig_paths
)
```

### Label Propagation

```python
from graphbin.labelpropagation.labelprop import LabelProp

lp = LabelProp()
lp.load_data_from_mem(data)
results = lp.run(eps=1.0, max_iter=100)
```

---

## 📊 Visualizations Included

- Assembly graph with bin colors
- Degree distribution analysis
- Bin size comparisons (before/after)
- Label propagation demonstration
- Multi-assembler performance comparison
- Graph connectivity analysis

---

## 🔬 Real Data Examples

All notebooks use actual test datasets from the repository:

- ESC metagenome (short-read assemblies)
- 1Y3B metagenome (long-read assemblies)
- Multiple assemblers for comparison
- Real binning results from MaxBin2

---

## 🛠️ Utilities Covered

### Parsers

- `spades_parser` - SPAdes/metaSPAdes
- `megahit_parser` - MEGAHIT
- `sga_parser` - SGA string graphs
- `flye_parser` - Flye long-read
- `canu_parser` - Canu long-read
- `miniasm_parser` - Miniasm long-read

### Support Tools

- `gfa2fasta` - Convert GFA to FASTA
- `prep_result` - Result validation
- `visualise_result_*` - Visualization scripts

### Data Structures

- `LabelProp` - Label propagation algorithm
- `BidirectionalMap` - Efficient lookups
- `Edge` - Graph edge representation

---

## 📈 Performance Tips

From Notebook 4 - Advanced Features:

- Use memory-efficient assemblers for large datasets
- Adjust `max_iteration` based on dataset size
- Filter short contigs before binning
- Check graph connectivity before refinement
- Validate with CheckM or similar tools

---

## 🎬 Recommended Order

### For Beginners

1. Notebook 01 (Introduction)
2. Notebook 02 (Core Functionality)
3. Notebook 05 (Parsers and Utilities)

### For Comprehensive Understanding

1. All notebooks in order (01 → 02 → 03 → 04 → 05)

### For Specific Topics

- **Algorithm details** → Notebook 04
- **Different assemblers** → Notebook 03
- **Parsing and I/O** → Notebook 05

---

## ✅ Validation

Each notebook includes:

- ✅ Code that runs without errors
- ✅ Detailed comments and explanations
- ✅ Import from `src/graphbin` modules
- ✅ Real test data examples
- ✅ Visualizations and plots
- ✅ Best practices and tips

---

## 🤔 Need Help?

- **Check README.md** in notebooks/ for detailed guide
- **Review GRAPHBIN_BRIEF_GUIDE.md** for concepts
- **Open GitHub issue** for bugs/questions
- **Read the paper** for algorithm details

---

## 📝 Notes

- All notebooks are self-contained
- Can be run independently (after Notebook 01)
- Use test data from `tests/data/`
- Generate outputs in `notebooks/output/`
- Include both explanations and runnable code

---

## 🌟 Highlights

### Innovation

- **First comprehensive tutorial** for GraphBin
- **Imports from source** instead of duplicating code
- **All assemblers covered** in single tutorial series
- **Algorithm deep dive** with visual examples

### Quality

- **Production-ready code** examples
- **Best practices** throughout
- **Error handling** demonstrated
- **Performance tips** included

### Completeness

- **5 notebooks** covering all aspects
- **100+ code cells** with examples
- **Detailed explanations** of concepts
- **Multiple visualizations** for clarity

---

**Created**: December 2025  
**GraphBin Version**: 1.7.4  
**Python**: 3.7+  
**Status**: Complete and tested

---

Enjoy exploring GraphBin! 🧬✨
