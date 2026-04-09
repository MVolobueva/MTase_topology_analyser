# MTase Topology Analyzer

A comprehensive tool for analyzing methyltransferase (MTase) protein topologies from PDB/mmCIF files. This application provides 2D and 3D visualizations, topology classification, and batch processing capabilities.

## Features

- **Topology Analysis**: Classifies MTase protein structures and identifies topological features
- **2D Visualization**: Generates publication-ready topology diagrams
- **3D Visualization**: Interactive 3D structure viewer for detailed inspection
- **Batch Processing**: Analyze multiple structures simultaneously
- **Web Interface**: User-friendly Streamlit-based interface
- **DSSP Integration**: Secondary structure assignment using embedded mkdssp tool

## Project Structure

```
MTase_topology_analyser/
├── app.py                 # Main Streamlit application
├── batch_analyze.py       # Batch processing module
├── classifier.py          # Topology classification logic
├── analyzer/              # Core analysis modules
│   ├── coordinates.py     # Coordinate handling
│   ├── core.py           # Core analysis engine
│   ├── topology.py       # Topology calculations
│   ├── visualization_2d.py # 2D plotting
│   └── visualization_3d.py # 3D visualization
├── components/           # UI components
│   ├── results_table.py
│   ├── sidebar.py
│   └── visualizations.py
├── pages/               # Streamlit multi-pages
│   ├── analysis_page.py
│   └── documentation_page.py
├── utils/               # Utility functions
│   └── helpers.py
├── input.csv            # Input data template
├── output.csv           # Analysis results
├── mkdssp              # DSSP executable (embedded)
└── requirements.txt     # Python dependencies
```

## Installation

### Prerequisites

- Python 3.11 or higher
- pip package manager

### Setup

1. Clone the repository:
```bash
git clone https://github.com/yourusername/MTase_topology_analyser.git
cd MTase_topology_analyser
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Verify the embedded DSSP tool:
```bash
chmod +x mkdssp
```

## Usage

### Web Interface (Recommended)

Run the Streamlit app:
```bash
streamlit run app.py
```

Then open your browser to `http://localhost:8501`

### Batch Processing

Run batch analysis from the command line:
```bash
python batch_analyze.py
```

To specify custom input/output files:
```bash
python batch_analyze.py input.csv output.csv
```

## Input Format

The batch analyzer accepts a CSV file with **two columns**: `ID` and `Type`

### Supported Types:

| Type | Description | Example ID |
|------|-------------|------------|
| `pdb` | PDB database entry | `3S1S` |
| `alphafold` | AlphaFold model | `A0A7R8ZSU6` |
| `file` | Local PDB file path | `/path/to/structure.pdb` |

### Example `input.csv`:

```csv
ID,Type
3S1S,pdb
A0A7R8ZSU6,alphafold
examples/1xyz.pdb,file
```

## Output

The analyzer generates `output.csv` with the following columns:
- `source_id` — Original ID from input
- `source_type` — Type (pdb/alphafold/file)
- `chain` — Chain identifier
- `found_motif` — Detected MTase motif
- `found_motif_position` — Residue position of the motif
- `full_secondary_elements` — Complete topology with helices and strands
- `strands_secondary_elements` — Topology with strands only
- `strand_directions` — Direction pattern (↑/↓)
- `class` — Topology class (A-F)
- `class_confidence` — Confidence score
- `class_reasons` — Reasoning for classification
- `strand_order` — Order of strands
- `has_s0`, `has_s-1` — Presence of special strands
- `gap_s6_s7` — Gap between S6 and S7
- `has_n_helix`, `has_c_helix` — Presence of terminal helices

## Dependencies

Key Python packages (see `requirements.txt`):
- streamlit
- pandas
- numpy
- plotly
- matplotlib
- biopython
- MDAnalysis

### Bundled Libraries

The repository includes pre-compiled libraries for DSSP functionality:
- `mkdssp` - DSSP executable
- `libcifpp.so.5` - CIF parsing library
- ICU libraries (for Unicode support)

## License

MIT License

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## Support

For issues or questions, please open a GitHub issue.

## Acknowledgments

- DSSP authors for secondary structure assignment
- The structural bioinformatics community

---

**Note**: This tool is designed for research purposes. Always validate results with experimental data.
