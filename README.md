#### Trso Advanced Genomic Analysis Pipeline

**Description:**
This script performs comprehensive genomic analysis, including BLAST queries, sequence alignment, phylogenetic analysis, primer design, variant calling, motif analysis, and sequence validation. It includes enhanced configuration management, error handling, and checkpointing to ensure robust and efficient execution.

**Usage:**
```bash
$ python script.py input.fasta --config config.yaml --output results
```

**Expected Outputs:**
- `results/report.html`: Interactive analysis report
- `results/alignment.aln`: Multiple sequence alignment in ALN format
- `results/alignment.fasta`: Multiple sequence alignment in FASTA format
- `results/tree.nwk`: Phylogenetic tree in NEWICK format
- `results/primers.json`: Designed primer sequences in JSON format
- `results/variants.vcf`: SNP and indel calls in VCF format
- `results/motifs.txt`: Detected motifs in TXT format
- `results/fastqc_output/`: FastQC validation reports

**Command-Line Arguments:**

| Argument               | Description                                                                                       | Default             |
|------------------------|---------------------------------------------------------------------------------------------------|---------------------|
| `query`                | Input sequence or FASTA file path                                                                 |                     |
| `-o`, `--output`       | Output directory                                                                                  | `results`           |
| `--config`             | Path to YAML/JSON config file                                                                     |                     |
| `--local_blast`        | Path to local BLAST database                                                                      |                     |
| `--program`            | BLAST program                                                                                   | `blastn`            |
| `--database`           | BLAST database                                                                                    | `nt`                |
| `--max_hits`           | Maximum hits to process                                                                         | `50`                |
| `--threads`            | Processing threads                                                                              | `4`                 |
| `--force`              | Overwrite existing results                                                                      |                     |
| `--e_value`            | E-value threshold for BLAST                                                                     | `1e-5`              |
| `--score_threshold`    | Score threshold for BLAST                                                                       | `50`                |
| `--keep_temp`          | Keep temporary files                                                                            |                     |
| `--resume`             | Resume from checkpoint                                                                          |                     |
| `--design_primers`     | Design primers for the query sequence                                                           |                     |
| `--log_level`          | Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)                                           | `INFO`              |
| `--email`              | NCBI email address                                                                              |                     |
| `--call_variants`      | Call variants using bcftools                                                                    |                     |
| `--reference`          | Reference genome for variant calling                                                            |                     |
| `--find_motifs`        | Find motifs using MEME                                                                          |                     |
| `--num_motifs`         | Number of motifs to find                                                                        | `5`                 |
| `--min_motif_width`    | Minimum motif width                                                                             | `6`                 |
| `--max_motif_width`    | Maximum motif width                                                                             | `12`                |
| `--advanced_tree`      | Use IQ-TREE for advanced phylogenetic analysis                                                  |                     |
| `--validate_sequence`  | Validate sequence using FastQC                                                                  |                     |

**Configuration File Example (`config.yaml`):**
```yaml
output: results
local_blast: /path/to/local/blastdb
program: blastn
database: nt
max_hits: 50
threads: 4
force: false
e_value: 1e-5
score_threshold: 50
keep_temp: false
resume: false
design_primers: true
log_level: INFO
email: your.email@example.com
call_variants: true
reference: /path/to/reference.fasta
find_motifs: true
num_motifs: 5
min_motif_width: 6
max_motif_width: 12
advanced_tree: true
validate_sequence: true
```

**Dependencies:**
Ensure the following dependencies are installed in your environment. You can install them using the `requirements.txt` file.

**Requirements File (`requirements.txt`):**
```
biopython
numpy
matplotlib
pandas
scipy
jinja2
logomaker
primer3-py
pyyaml
plotly
subprocess32
fastqc
meme
iqtree
bcftools
```

**Installation:**
1. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```

2. Ensure external tools are installed and available in your PATH:
   - `clustalw2`
   - `blastn`
   - `bcftools`
   - `meme`
   - `iqtree`
   - `fastqc`

**Example Workflow:**
1. Prepare your input FASTA file (`input.fasta`).
2. Create a configuration file (`config.yaml`) with desired parameters.
3. Run the script:
   ```bash
   python script.py input.fasta --config config.yaml --output results
   ```

**Logging:**
Logs are written to `results/pipeline.log` and also printed to the console. The log level can be adjusted using the `--log_level` argument.

**Checkpointing:**
The script uses checkpoints to resume from previous stages if interrupted. Use the `--resume` flag to continue from the last saved checkpoint.

**Interactive Report:**
The final report (`results/report.html`) includes interactive visualizations such as phylogenetic trees and sequence conservation logos.

```

### Notes:
- Ensure that the external tools (`clustalw2`, `blastn`, `bcftools`, `meme`, `iqtree`, `fastqc`) are installed and available in your system's PATH.
- The script assumes that the input FASTA file is correctly formatted and contains valid nucleotide sequences.
- The `bcftools` and `fastqc` tools require their respective dependencies and configurations.
