# GraphBin Tutorial Notebooks

This directory contains comprehensive Jupyter notebooks demonstrating all functionalities of GraphBin, a metagenomic contig bin refinement tool.

## 📚 Notebook Series

### [01_GraphBin_Introduction_and_Setup.ipynb](01_GraphBin_Introduction_and_Setup.ipynb)

**Introduction and Setup**

- GraphBin overview and installation verification
- Understanding key concepts (contigs, bins, assembly graphs)
- Exploring test data structure
- Understanding input file formats (FASTA, GFA, CSV)
- Parsing assembly graphs and binning results
- Module structure overview

**Prerequisites**: None  
**Time**: 15-20 minutes

---

### [02_GraphBin_Core_Functionality.ipynb](02_GraphBin_Core_Functionality.ipynb)

**Core Functionality - Binning Refinement**

- Running GraphBin programmatically with SPAdes data
- Analyzing initial binning statistics
- Examining assembly graph connectivity
- Running the refinement algorithm
- Comparing initial vs refined results
- Identifying newly binned and corrected contigs
- Visualizing improvements
- Understanding convergence

**Prerequisites**: Notebook 01  
**Time**: 30-40 minutes

---

### [03_GraphBin_Multiple_Assemblers.ipynb](03_GraphBin_Multiple_Assemblers.ipynb)

**Multiple Assemblers Support**

- Comparing 6 supported assemblers:
  - **Short-read**: SPAdes, MEGAHIT, SGA
  - **Long-read**: Flye, Canu, Miniasm
- Running GraphBin with different assemblers
- Analyzing and comparing refinement results
- Understanding assembler-specific differences
- Choosing the right assembler for your data

**Prerequisites**: Notebook 01  
**Time**: 45-60 minutes

---

### [04_GraphBin_Advanced_Features.ipynb](04_GraphBin_Advanced_Features.ipynb)

**Advanced Features - Label Propagation and Graph Analysis**

- Deep dive into label propagation algorithm
- Simple example demonstrating the algorithm
- Analyzing real assembly graph properties
- Graph connectivity and binning relationships
- Parameter sensitivity analysis
- Algorithm workflow step-by-step
- Edge cases and limitations
- Performance optimization tips

**Prerequisites**: Notebooks 01, 02  
**Time**: 40-50 minutes

---

### [05_GraphBin_Parsers_and_Utilities.ipynb](05_GraphBin_Parsers_and_Utilities.ipynb)

**Parsers, Utilities, and Visualization**

- Understanding parser modules for each assembler
- Parsing assembly graphs (GFA, ASQG formats)
- Using format conversion utilities (GFA to FASTA)
- Writing and validating outputs
- Visualization utilities
- Custom visualization examples
- BidirectionalMap data structure
- Best practices and troubleshooting

**Prerequisites**: Notebook 01  
**Time**: 30-40 minutes

---

## 🚀 Getting Started

### Quick Start

1. **Clone the repository** (if you haven't already):

   ```bash
   git clone https://github.com/metagentools/GraphBin.git
   cd GraphBin
   ```

2. **Install dependencies**:

   ```bash
   # Using conda (recommended)
   conda env create -f environment.yml
   conda activate graphbin

   # OR using pip
   pip install -r requirements.txt
   pip install jupyter matplotlib pandas
   ```

3. **Launch Jupyter**:

   ```bash
   cd notebooks
   jupyter notebook
   ```

4. **Start with Notebook 01** and proceed sequentially.

### Installation Verification

Run this in a notebook cell to verify installation:

```python
import sys
from pathlib import Path

# Add GraphBin to path
repo_root = Path('..').resolve()
sys.path.insert(0, str(repo_root / 'src'))

import graphbin
from graphbin import __version__
print(f"✓ GraphBin {__version__} loaded successfully!")
```

---

## 📊 What You'll Learn

### Core Concepts

- Metagenomic binning and refinement
- Assembly graph structures and connectivity
- Label propagation algorithm
- Graph-based optimization

### Practical Skills

- Running GraphBin with different assemblers
- Analyzing and interpreting results
- Parameter tuning for different datasets
- Visualizing assembly graphs and binning
- Troubleshooting common issues

### Advanced Topics

- Algorithm implementation details
- Graph analysis and statistics
- Parser internals
- Custom visualization
- Performance optimization

---

## 📁 Test Data

All notebooks use test datasets located in `tests/data/`:

- `ESC_metaSPAdes/` - SPAdes short-read assembly
- `ESC_MEGAHIT/` - MEGAHIT short-read assembly
- `ESC_SGA/` - SGA string graph assembly
- `1Y3B_Flye/` - Flye long-read assembly
- `1Y3B_Canu/` - Canu long-read assembly
- `1Y3B_Miniasm/` - Miniasm long-read assembly

Each dataset includes:

- Assembly graph file (GFA or ASQG)
- Contigs file (FASTA)
- Initial binning results (CSV)
- Paths file (where applicable)

---

## 🔧 Requirements

### Software

- Python 3.7+
- Jupyter Notebook or JupyterLab
- GraphBin (installed from this repository)

### Python Packages

- **Core dependencies**:

  - `igraph` - Graph analysis
  - `cogent3` - Sequence handling
  - `cairocffi` - Graphics support
  - `click` - CLI framework

- **Notebook dependencies**:
  - `jupyter` - Notebook environment
  - `matplotlib` - Plotting
  - `pandas` - Data analysis
  - `numpy` - Numerical computing

### System Requirements

- **RAM**: 2-4 GB minimum, 8+ GB recommended for large datasets
- **Disk**: 500 MB for repository + space for outputs
- **OS**: Linux, macOS, or Windows (with WSL recommended)

---

## 💡 Tips for Success

1. **Run notebooks in order** - Each builds on concepts from previous ones
2. **Execute all cells** - Don't skip cells; some set up state for later cells
3. **Experiment** - Modify parameters and see what happens
4. **Use your own data** - After completing tutorials, try with your datasets
5. **Check the logs** - GraphBin produces detailed logs for debugging
6. **Visualize results** - Plots help understand what's happening

---

## 🐛 Troubleshooting

### Common Issues

**Issue**: Import errors

```
Solution: Ensure GraphBin is installed and path is set correctly
Check: sys.path should include the src directory
```

**Issue**: Test data not found

```
Solution: Run notebooks from the notebooks/ directory
Check: Verify test data exists in tests/data/
```

**Issue**: Out of memory errors

```
Solution: Use smaller test datasets or increase available RAM
Try: Close other applications, use memory-efficient assemblers
```

**Issue**: Plots not showing

```
Solution: Add %matplotlib inline to notebook cells
Check: Matplotlib backend configuration
```

**Issue**: Slow execution

```
Solution: Reduce max_iteration parameter
Try: Use smaller subsets of test data for exploration
```

---

## 📖 Additional Resources

### Documentation

- **GraphBin Docs**: https://graphbin.readthedocs.io/
- **GRAPHBIN_BRIEF_GUIDE.md**: Non-technical overview in repository root

### Scientific Paper

- Mallawaarachchi et al., Bioinformatics 2020
- DOI: [10.1093/bioinformatics/btaa180](https://doi.org/10.1093/bioinformatics/btaa180)

### Community

- **GitHub Issues**: https://github.com/metagentools/GraphBin/issues
- **Discussions**: https://github.com/metagentools/GraphBin/discussions

---

## 🎯 Learning Path

### Beginner Track (Essential)

1. Notebook 01 - Introduction ⭐
2. Notebook 02 - Core Functionality ⭐
3. Notebook 05 - Parsers and Utilities

### Advanced Track (Comprehensive)

1. Notebook 01 - Introduction ⭐
2. Notebook 02 - Core Functionality ⭐
3. Notebook 03 - Multiple Assemblers
4. Notebook 04 - Advanced Features
5. Notebook 05 - Parsers and Utilities

### Expert Track (Deep Dive)

Complete all notebooks + modify code to:

- Implement custom label propagation variations
- Create new visualization methods
- Optimize for your specific data types
- Integrate with other metagenomic tools

---

## 🤝 Contributing

Found an issue or want to improve the notebooks?

1. Open an issue on GitHub
2. Submit a pull request with improvements
3. Share your custom examples

---

## 📝 Citation

If you use GraphBin in your research, please cite:

```bibtex
@article{mallawaarachchi2020graphbin,
  title={GraphBin: refined binning of metagenomic contigs using assembly graphs},
  author={Mallawaarachchi, Vijini and Wickramarachchi, Anuradha and Lin, Yu},
  journal={Bioinformatics},
  volume={36},
  number={11},
  pages={3307--3313},
  year={2020},
  publisher={Oxford University Press}
}
```

---

## 📄 License

These notebooks are part of the GraphBin project and are distributed under the BSD-3-Clause License.

---

## ✨ Acknowledgments

GraphBin is funded by an Essential Open Source Software for Science Grant from the Chan Zuckerberg Initiative.

---

**Happy Learning! 🧬🔬**

For questions or feedback, please open an issue on GitHub.
