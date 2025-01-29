#### Trso Advanced Genomic Analysis Pipeline

The script performs comprehensive genomic analysis, including BLAST queries, sequence alignment, phylogenetic analysis, primer design, variant calling, motif analysis, and sequence validation. It supports various BLAST programs (`blastn`, `blastx`, `blastp`) and handles both nucleotide and protein sequences.

#### **Features:**

1. **BLAST Searches:**
   - **Program Support:** Handles `blastn`, `blastx`, and `blastp`.
   - **Local and Remote Execution:** Supports both local and remote BLAST searches.
   - **Caching Mechanism:** Caches BLAST results to avoid redundant computations.
   - **Enhanced Parsing:** Parses BLAST results with enhanced filtering based on `e_value` and `score_threshold`.

2. **Sequence Fetching:**
   - **Parallel Requests:** Fetches sequences from GenBank in parallel using `ThreadPoolExecutor`, improving efficiency.
   - **Validation:** Validates fetched sequences and logs any failures.

3. **Sequence Alignment:**
   - **Multiple Sequence Alignment:** Uses ClustalW for aligning sequences.
   - **Type Handling:** Automatically handles DNA or protein sequences based on the BLAST program used (`dna` for `blastn`, `protein` for `blastx` and `blastp`).

4. **Phylogenetic Analysis:**
   - **Distance-Based Trees:** Constructs phylogenetic trees using `DistanceTreeConstructor` for simple distance-based methods.
   - **Advanced Trees:** Uses IQ-TREE for advanced phylogenetic analysis, supporting both DNA and protein sequences.

5. **Variant Calling:**
   - **Using bcftools:** Calls variants from aligned sequences using `bcftools`.
   - **Reference Genome Required:** Requires a reference genome for variant calling.

6. **Motif Analysis:**
   - **Using MEME:** Detects motifs in sequences using the MEME suite.
   - **Customizable Parameters:** Allows customization of the number of motifs and motif width.

7. **Primer Design:**
   - **Using primer3:** Designs primers for the query sequence using the primer3 library.
   - **Parameters:** Configurable primer design parameters such as size, melting temperature, and GC content.

8. **Sequence Validation:**
   - **Using FastQC:** Validates input sequences for contamination or low-quality regions using FastQC.
   - **Report Generation:** Generates FastQC reports for quality assessment.

9. **Visualization:**
   - **Interactive Phylogenetic Trees:** Generates interactive phylogenetic tree visualizations using Plotly.
   - **Sequence Conservation Logos:** Creates sequence conservation logos using logomaker.

10. **Configuration Management:**
    - **Command-Line Arguments:** Accepts various command-line arguments for configuration.
    - **Config File Support:** Loads configuration from a YAML/JSON file, allowing for flexible settings.

11. **Checkpointing and State Recovery:**
    - **Resumable Runs:** Saves checkpoints to resume interrupted runs without repeating completed steps.
    - **Force Rerun:** Option to force rerunning the pipeline even if a complete run exists.

12. **Logging and Error Handling:**
    - **Detailed Logging:** Comprehensive logging for each step of the pipeline.
    - **Error Handling:** Robust error handling to manage exceptions and ensure graceful failure.

13. **Output Generation:**
    - **HTML Report:** Compiles a comprehensive HTML report with visualizations and results.
    - **Alignment Files:** Saves multiple sequence alignments in ALN and FASTA formats.
    - **Tree Files:** Saves phylogenetic trees in NEWICK format.
    - **Variant Files:** Saves variant calls in VCF format.
    - **Motif Files:** Saves detected motifs in TXT format.
    - **Primer Files:** Saves designed primers in JSON format.
    - **FastQC Reports:** Includes FastQC validation reports.

#### **Command-Line Arguments:**

| Argument               | Description                                                                                       | Default             |
|------------------------|---------------------------------------------------------------------------------------------------|---------------------|
| `query`                | Input sequence or FASTA file path                                                                 |                     |
| `-o`, `--output`       | Output directory                                                                                  | `results`           |
| `--config`             | Path to YAML/JSON config file                                                                     |                     |
| `--local_blast`        | Path to local BLAST database                                                                      |                     |
| `--program`            | BLAST program (`blastn`, `blastx`, `blastp`)                                                    | `blastn`            |
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

#### **Modules:**

1. **PipelineState:**
   - Manages checkpoints and state recovery to ensure resumable runs.

2. **ConfigManager:**
   - Handles configuration loading from command-line arguments and YAML/JSON config files.

3. **NcbiHandler:**
   - Manages NCBI Entrez interactions and sets the NCBI email address.

4. **BlastRunner:**
   - Executes BLAST searches with caching and resume support.
   - Supports both local and remote BLAST searches.

5. **SequenceAnalyzer:**
   - Fetches sequences from GenBank with parallel requests.
   - Aligns sequences using ClustalW.
   - Translates nucleotide sequences to protein sequences if needed.

6. **VariantCaller:**
   - Calls variants using `bcftools`.

7. **MotifAnalyzer:**
   - Finds motifs using the MEME suite.

8. **PhylogeneticAnalyzer:**
   - Constructs phylogenetic trees using IQ-TREE for advanced analysis.

9. **SequenceValidator:**
   - Validates sequences using FastQC.

10. **TreeVisualizer:**
    - Generates interactive phylogenetic tree visualizations using Plotly.

11. **LogoGenerator:**
    - Generates sequence conservation logos using logomaker.

12. **VisualizationEngine:**
    - Compiles a comprehensive HTML report with visualizations and results.

13. **AnalysisPipeline:**
    - Orchestrates the entire analysis pipeline, managing stages and state transitions.

#### **Requirements:**

- **Python Packages:**
  ```txt
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

- **External Tools:**
  - `clustalw2`
  - `blastn`, `blastx`, `blastp`
  - `bcftools`
  - `meme`
  - `iqtree`
  - `fastqc`

#### **Usage Example:**

```bash
$ python script.py input.fasta --config config.yaml --output results
```

#### **Expected Outputs:**

- `results/report.html`: Interactive analysis report
- `results/alignment.aln`: Multiple sequence alignment in ALN format
- `results/alignment.fasta`: Multiple sequence alignment in FASTA format
- `results/tree.nwk`: Phylogenetic tree in NEWICK format
- `results/primers.json`: Designed primer sequences in JSON format
- `results/variants.vcf`: SNP and indel calls in VCF format
- `results/motifs.txt`: Detected motifs in TXT format
- `results/fastqc_output/`: FastQC validation reports
