```markdown
# Advanced Genomic Analysis Pipeline Documentation

## Overview
A comprehensive genomic analysis pipeline performing:
- BLAST sequence searches (remote/local)
- GenBank sequence retrieval
- Multiple sequence alignment (ClustalW)
- Phylogenetic tree construction
- Primer design (primer3)
- Interactive HTML reporting

## Installation

### Requirements
- Python 3.8+
- System dependencies:
  - ClustalW2 (`sudo apt-get install clustalw`)
  - BLAST+ (for local searches - `sudo apt-get install ncbi-blast+`)

### Python Packages
```bash
pip install biopython pyyaml primer3 logomaker matplotlib jinja2
```

## Usage

### Basic Execution
```bash
python script.py input.fasta --config config.yaml --output results
```

### Key Arguments
| Argument          | Description                                  | Default       |
|-------------------|----------------------------------------------|---------------|
| `query`           | Input FASTA file or raw sequence            | Required      |
| `--output`        | Output directory                            | `results`     |
| `--config`        | Configuration file (YAML/JSON)              | Optional      |
| `--local_blast`   | Path to local BLAST database                | Remote search |
| `--max_hits`      | Maximum BLAST hits to process               | 50            |
| `--threads`       | Parallel processing threads                 | 4             |
| `--force`         | Overwrite existing results                  | False         |
| `--resume`        | Resume from last checkpoint                 | False         |

### Example Config File (`config.yaml`)
```yaml
program: blastn
database: nt
e_value: 0.001
max_hits: 100
email: user@example.com
primer_params:
  size_range: [18, 25]
  tm_temp: 60.0
```

## Pipeline Components

### 1. BLAST Search Module
- Remote (NCBI) or local execution
- XML result caching with SHA256 hashing
- Configurable parameters:
  - Program (`blastn`, `blastp`, etc.)
  - Database (`nt`, `nr`, custom)
  - E-value threshold
  - Percent identity

### 2. Sequence Retrieval
- Parallel fetching from GenBank
- FASTA format output
- Entrez email validation

### 3. Multiple Sequence Alignment
- ClustalW2 implementation
- Output formats: CLUSTAL, FASTA
- Temp file cleanup option

### 4. Phylogenetic Analysis
- Distance matrix calculation
- UPGMA/NJ tree construction
- Newick format output

### 5. Primer Design
- primer3 binding parameters:
  ```python
  {
      'PRIMER_OPT_SIZE': 20,
      'PRIMER_MIN_SIZE': 18,
      'PRIMER_MAX_SIZE': 25,
      'PRIMER_OPT_TM': 60.0,
      'PRIMER_MIN_TM': 58.0,
      'PRIMER_MAX_TM': 62.0
  }
  ```
- JSON output with primer details

### 6. Reporting System
- Interactive HTML report with:
  - Sequence conservation logo
  - Phylogenetic tree visualization
  - Primer table
  - BLAST summary statistics
- Plotly integration for interactive elements

## Input/Output Specifications

### Input
- FASTA file or raw nucleotide sequence
- YAML/JSON configuration file (optional)

### Output Structure
```
results/
├── alignment.aln       # Multiple sequence alignment
├── primers.json        # Designed primers in JSON format
├── report.html         # Interactive analysis report
├── local_blast.xml     # BLAST results (if local)
├── blast_cache/        # Cached BLAST results
└── pipeline.log        # Detailed execution log
```

## Checkpoint System
- Automatic state preservation at key stages:
  1. Initialization
  2. BLAST completion
  3. Sequence retrieval
  4. Alignment
  5. Tree construction
  6. Primer design
- Resume with `--resume` flag
- JSON state storage in `checkpoint.json`

## Advanced Configuration

### Environment Variables
```bash
export NCBI_EMAIL="user@institute.edu"  # Required for Entrez access
```

### Customizing Visualizations
Modify the Jinja2 template in `VisualizationEngine.HTML_TEMPLATE` to:
- Add custom CSS/JavaScript
- Modify report layout
- Integrate additional visualization libraries

## Troubleshooting

### Common Issues
1. **Missing Dependencies**:
   - Ensure ClustalW2 and BLAST+ are in PATH
   - Verify Python package versions

2. **NCBI Access Errors**:
   - Validate Entrez email configuration
   - Check API rate limits (3 requests/sec)

3. **Alignment Failures**:
   - Verify input sequence homogeneity
   - Check ClustalW installation

4. **Primer Design Warnings**:
   - Adjust melting temperature parameters
   - Modify product size ranges

### Debugging Tips
```bash
--log_level DEBUG     # Enable verbose logging
--keep_temp           # Preserve intermediate files
```

## License
MIT License - See included LICENSE file

---

This documentation covers key aspects of the genomic pipeline. For implementation details, refer to inline code comments and class docstrings.
```
