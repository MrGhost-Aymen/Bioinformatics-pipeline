from io import StringIO
import argparse
import base64
import json
import logging
import os
import sys
import hashlib
import shutil
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import yaml  # Requires PyYAML package
from Bio import SeqIO, Entrez, Phylo, AlignIO
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastxCommandline, NcbiblastpCommandline
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from jinja2 import Template
import primer3
import logomaker
from Bio.Seq import Seq
from collections import defaultdict
import plotly.graph_objects as go
import subprocess
import time

# Configure logging
logger = logging.getLogger(__name__)

__version__ = "1.0.4"  # Updated version

class PipelineState:
    """Manages pipeline checkpoints and state recovery."""
    CHECKPOINT_STAGES = [
        'init', 'blast', 'fetch', 'align', 'tree', 'primers', 'variants',
        'motifs', 'validation', 'complete'
    ]

    def __init__(self, args):
        self.args = args
        self.checkpoint_file = Path(args.output) / "checkpoint.json"
        self.state = self._load_checkpoint() or {
            "stage": "init",
            "data": {},
            "timestamp": datetime.now().isoformat()
        }

    def _load_checkpoint(self):
        """Load existing checkpoint if available and resume flag is set."""
        if self.args.resume and self.checkpoint_file.exists():
            try:
                with open(self.checkpoint_file) as f:
                    state = json.load(f)
                    if state['stage'] in self.CHECKPOINT_STAGES:
                        logger.info(f"Resuming from {state['stage']} stage")
                        return state
            except Exception as e:
                logger.warning(f"Checkpoint load failed: {e}")
        return None

    def save(self, stage, data):
        """Persist current state to checkpoint file."""
        if stage not in self.CHECKPOINT_STAGES:
            raise ValueError(f"Invalid stage: {stage}")
        self.state.update({
            "stage": stage,
            "data": data,
            "timestamp": datetime.now().isoformat()
        })
        with open(self.checkpoint_file, "w") as f:
            json.dump(self.state, f, indent=2)

class ConfigManager:
    """Handles configuration loading from file and command-line."""
    @staticmethod
    def load(config_path, args):
        """Merge config file settings with command-line arguments."""
        config_path = Path(config_path)
        if not config_path.exists():
            raise FileNotFoundError(f"Config file {config_path} not found")
        with open(config_path) as f:
            config = yaml.safe_load(f) if config_path.suffix in ['.yaml', '.yml'] else json.load(f)
        for key, value in config.items():
            if hasattr(args, key):
                setattr(args, key, value)
            else:
                logger.warning(f"Ignoring unknown config key: {key}")
        return args

class NcbiHandler:
    """Manages NCBI Entrez interactions and credentials."""
    @staticmethod
    def set_email(email):
        """Validate and set NCBI email address."""
        Entrez.email = email or os.getenv("NCBI_EMAIL")
        if not Entrez.email:
            raise ValueError("NCBI email required. Set via --email or NCBI_EMAIL env var")
        logger.info(f"Using NCBI email: {Entrez.email}")

class BlastRunner:
    """Handles BLAST searches with local/remote execution and caching."""
    CACHE_DIR = "blast_cache"

    @classmethod
    def run(cls, query, args):
        """Execute BLAST with caching and resume support."""
        cache_key = cls._cache_key(query, args)
        cached = cls._check_cache(cache_key, args)
        if cached and not args.force:
            logger.info("Using cached BLAST results")
            with open(cached, "r") as f:
                return StringIO(f.read())  # Return cached XML content as a file-like object
        result = cls._run_remote(query, args) if not args.local_blast else cls._run_local(query, args)
        cls._save_cache(cache_key, result, args)
        return result

    @classmethod
    def _cache_key(cls, query, args):
        """Generate SHA256 hash for BLAST query caching."""
        key_data = {
            'query': query,
            'program': args.program,
            'database': args.database,
            'params': {
                'e_value': args.e_value,
                'max_hits': args.max_hits
            }
        }
        return hashlib.sha256(json.dumps(key_data, sort_keys=True).encode()).hexdigest()

    @classmethod
    def _check_cache(cls, cache_key, args):
        """Check for existing cached BLAST results."""
        cache_dir = Path(args.output) / cls.CACHE_DIR
        cache_file = cache_dir / f"{cache_key}.xml"
        return cache_file if cache_file.exists() else None

    @classmethod
    def _save_cache(cls, cache_key, result, args):
        """Save BLAST results to cache."""
        cache_dir = Path(args.output) / cls.CACHE_DIR
        cache_dir.mkdir(exist_ok=True)
        cache_file = cache_dir / f"{cache_key}.xml"
        with open(cache_file, "w") as f:
            f.write(result.read())
        result.seek(0)
        return cache_file

    @staticmethod
    def _run_remote(query, args):
        """Execute remote BLAST via NCBI WWW service."""
        logger.info("Running remote BLAST via NCBI...")
        try:
            time.sleep(0.34)  # Rate limit: ~3 requests/sec without API key
            return NCBIWWW.qblast(
                args.program,
                args.database,
                query,
                hitlist_size=args.max_hits,
                alignments=args.max_hits,
                descriptions=args.max_hits
            )
        except Exception as e:
            logger.error(f"Remote BLAST failed: {e}")
            raise

    @staticmethod
    def _run_local(query, args):
        """Execute local BLAST against specified database."""
        logger.info(f"Running local BLAST against {args.local_blast}")
        output_file = Path(args.output) / "local_blast.xml"
        if not Path(args.local_blast).exists():
            raise FileNotFoundError(f"BLAST database {args.local_blast} not found")
        blast_cline = None
        if args.program == "blastn":
            blast_cline = NcbiblastnCommandline(
                query=query,
                db=args.local_blast,
                outfmt=5,
                out=output_file,
                num_threads=args.threads
            )
        elif args.program == "blastx":
            blast_cline = NcbiblastxCommandline(
                query=query,
                db=args.local_blast,
                outfmt=5,
                out=output_file,
                num_threads=args.threads
            )
        elif args.program == "blastp":
            blast_cline = NcbiblastpCommandline(
                query=query,
                db=args.local_blast,
                outfmt=5,
                out=output_file,
                num_threads=args.threads
            )
        else:
            raise ValueError(f"Unsupported BLAST program: {args.program}")
        logger.info(f"Executing: {blast_cline}")
        stdout, stderr = blast_cline()
        if stderr:
            logger.error(f"BLAST error: {stderr}")
            raise RuntimeError(stderr)
        return open(output_file)

class SequenceAnalyzer:
    """Handles sequence fetching, alignment, and analysis."""
    @staticmethod
    def fetch_sequences(hit_ids, args):
        """Fetch sequences from GenBank with parallel requests."""
        if not hit_ids:
            raise ValueError("No valid BLAST hits to process")
        logger.info(f"Fetching {len(hit_ids)} sequences from GenBank")
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            futures = [executor.submit(SequenceAnalyzer._fetch_single, hid) for hid in hit_ids]
            results = []
            for future in futures:
                result = future.result()
                if result:
                    results.append(result)
                    logger.debug(f"Fetched {result.id}")
                else:
                    logger.warning("Failed to fetch sequence")
                time.sleep(0.34)  # Rate limit for NCBI compliance
            return results

    @staticmethod
    def _fetch_single(hit_id):
        """Fetch single sequence from GenBank."""
        try:
            handle = Entrez.efetch(db="nucleotide", id=hit_id, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            return record
        except Exception as e:
            logger.error(f"Failed to fetch {hit_id}: {e}")
            return None

    @staticmethod
    def align(sequences, args):
        """Perform multiple sequence alignment with MUSCLE."""
        if not sequences:
            raise ValueError("No sequences available for alignment")
        logger.info(f"Aligning {len(sequences)} sequences using MUSCLE")
        input_file = Path(args.output) / "alignment_input.fasta"
        output_file = Path(args.output) / "alignment.aln"
        SeqIO.write(sequences, input_file, "fasta")
        muscle_cmd = ["muscle", "-in", str(input_file), "-out", str(output_file)]
        logger.info(f"Executing: {' '.join(muscle_cmd)}")
        try:
            subprocess.run(muscle_cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"MUSCLE failed with error: {e}")
            raise RuntimeError("MUSCLE execution failed.")
        finally:
            if not args.keep_temp:
                input_file.unlink(missing_ok=True)
        return AlignIO.read(output_file, "fasta")

    @staticmethod
    def translate_sequences(sequences):
        """Translate nucleotide sequences to protein sequences."""
        translated_sequences = []
        for seq_record in sequences:
            try:
                translated_seq = seq_record.seq.translate(to_stop=True)
                translated_record = SeqIO.SeqRecord(translated_seq, id=seq_record.id, description=seq_record.description)
                translated_sequences.append(translated_record)
            except Exception as e:
                logger.warning(f"Translation failed for {seq_record.id}: {e}")
        return translated_sequences

class VariantCaller:
    """Handles variant calling using external tools like bcftools."""
    @staticmethod
    def call_variants(reference, bam_file, args):
        """Call variants using bcftools."""
        vcf_file = Path(args.output) / "variants.vcf"
        mpileup_cmd = ["bcftools", "mpileup", "-Ou", "-f", str(reference), str(bam_file)]
        call_cmd = ["bcftools", "call", "-mv", "-Ov", "-o", str(vcf_file)]
        logger.info(f"Executing variant calling: mpileup -> call")
        try:
            mpileup = subprocess.run(mpileup_cmd, capture_output=True, check=True)
            subprocess.run(call_cmd, input=mpileup.stdout, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"Variant calling failed: {e}")
            raise RuntimeError("Variant calling execution failed.")
        return vcf_file

class MotifAnalyzer:
    """Handles motif analysis using external tools like MEME."""
    @staticmethod
    def find_motifs(fasta_file, args):
        """Find motifs using MEME."""
        meme_output_dir = Path(args.output) / "meme_output"
        meme_output_dir.mkdir(exist_ok=True)
        meme_cmd = [
            "meme", str(fasta_file), "-oc", str(meme_output_dir),
            "-nmotifs", str(args.num_motifs), "-minw", str(args.min_motif_width),
            "-maxw", str(args.max_motif_width)
        ]
        logger.info(f"Executing motif analysis: {' '.join(meme_cmd)}")
        try:
            subprocess.run(meme_cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"MEME failed with error: {e}")
            raise RuntimeError("MEME execution failed.")
        motifs_file = meme_output_dir / "meme.txt"
        with open(motifs_file) as f:
            motifs = f.read()
        return motifs

class PhylogeneticAnalyzer:
    """Handles advanced phylogenetic analysis using IQ-TREE."""
    @staticmethod
    def construct_tree(alignment_file, args):
        """Construct phylogenetic tree using IQ-TREE."""
        tree_file = Path(args.output) / "tree.nwk"
        iqtree_cmd = [
            "iqtree", "-s", str(alignment_file), "-st", "PROTEIN",
            "-pre", str(tree_file.stem), "-nt", "AUTO"
        ]
        logger.info(f"Executing phylogenetic tree construction: {' '.join(iqtree_cmd)}")
        try:
            subprocess.run(iqtree_cmd, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"IQ-TREE failed with error: {e}")
            raise RuntimeError("IQ-TREE execution failed.")
        with open(tree_file) as f:
            tree_newick = f.read()
        return tree_newick

class SequenceValidator:
    """Handles sequence validation using FastQC."""
    @staticmethod
    def validate_sequence(fasta_file, args):
        """Validate sequence using FastQC."""
        fastqc_output_dir = Path(args.output) / "fastqc_output"
        fastqc_output_dir.mkdir(exist_ok=True)
        fastqc_cmd = [
            "fastqc", str(fasta_file), "-o", str(fastqc_output_dir), "--extract"
        ]
        logger.info(f"Executing sequence validation: {' '.join(fastqc_cmd)}")
        try:
            subprocess.run(fastqc_cmd, check=True)
            return fastqc_output_dir
        except subprocess.CalledProcessError as e:
            logger.error(f"FastQC failed with error: {e}")
            raise RuntimeError("FastQC execution failed.")

class TreeVisualizer:
    """Generates interactive phylogenetic tree visualizations."""
    @staticmethod
    def generate_interactive_tree(newick_str):
        """Generate Plotly-based phylogenetic tree visualization."""
        if not newick_str:
            return "Interactive tree placeholder"
        try:
            tree = Phylo.read(StringIO(newick_str), "newick")
            fig = go的比例(
                x=[],
                y=[],
                mode='markers+lines',
                text=[],
                hoverinfo='text',
                marker=dict(color='blue', size=10)
            )])

            def draw_clade(clade, parent_x, parent_y, depth, y_offset=0):
                x_values = [parent_x]
                y_values = [parent_y]
                clade_x = parent_x + 1
                clade_y = y_offset - depth
                x_values.append(clade_x)
                y_values.append(clade_y)
                name = clade.name if clade.name else "Internal Node"
                fig.data[0].x += tuple(x_values)
                fig.data[0].y += tuple(y_values)
                fig.data[0].text += (name,)
                if clade.clades:
                    n_clades = len(clade.clades)
                    for i, subclade in enumerate(clade.clades):
                        child_y_offset = y_offset + (n_clades - 1) / 2 - i
                        draw_clade(subclade, clade_x, clade_y, depth + 1, child_y_offset)

            draw_clade(tree.root, 0, 0, 0)
            fig.update_layout(title='Phylogenetic Tree', showlegend=False, autosize=True)
            return fig.to_html(full_html=False)
        except Exception as e:
            logger.error(f"Tree visualization failed: {e}")
            return "Failed to generate interactive tree"

class LogoGenerator:
    """Generates sequence conservation logos."""
    @staticmethod
    def generate_conservation_logo(alignment, args):
        """Generate sequence logo using logomaker."""
        try:
            logo_path = Path(args.output) / "logo.png"
            df = logomaker.alignment_to_matrix([str(rec.seq) for rec in alignment])
            plt.figure(figsize=(20, 5))
            logomaker.Logo(df)
            plt.savefig(logo_path)
            plt.close()
            with open(logo_path, "rb") as f:
                return base64.b64encode(f.read()).decode()
        except Exception as e:
            logger.error(f"Logo generation failed: {e}")
            return None

class VisualizationEngine:
    """Generates interactive HTML reports with visualization components."""
    HTML_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <title>Genomic Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1, h2 { color: #2c3e50; }
        pre { background: #f4f4f4; padding: 10px; border-radius: 5px; }
        img { max-width: 100%; }
    </style>
</head>
<body>
    <h1>Genomic Analysis Report</h1>
    <h2>Phylogenetic Tree</h2>
    {{ tree_html | safe }}
    <h2>Sequence Conservation Logo</h2>
    {% if logo_img %}
        <img src="data:image/png;base64,{{ logo_img }}" alt="Conservation Logo">
    {% else %}
        <p>No conservation logo generated.</p>
    {% endif %}
    <h2>Designed Primers</h2>
    <ul>
        {% for primer in primers %}
        <li>Forward: {{ primer.forward }} | Reverse: {{ primer.reverse }}</li>
        {% endfor %}
    </ul>
    <h2>Variants</h2>
    <pre>{{ variants | tojson(indent=2) }}</pre>
    <h2>Motifs</h2>
    <pre>{{ motifs }}</pre>
    <h2>FastQC Report</h2>
    <pre>{{ fastqc_report | tojson(indent=2) }}</pre>
    <h2>Summary</h2>
    <pre>{{ summary | tojson(indent=2) }}</pre>
</body>
</html>
"""
    @staticmethod
    def generate_report(results, args):
        """Generate comprehensive HTML report."""
        logger.info("Compiling final report")
        report_path = Path(args.output) / "report.html"
        logo_img = LogoGenerator.generate_conservation_logo(results['alignment'], args)
        tree_html = TreeVisualizer.generate_interactive_tree(results.get('tree_newick'))
        template = Template(VisualizationEngine.HTML_TEMPLATE)
        html = template.render(
            tree_html=tree_html,
            logo_img=logo_img,
            primers=results.get('primers', []),
            variants=results.get('variants', {}),
            motifs=results.get('motifs', ''),
            fastqc_report=results.get('fastqc_report', ''),
            summary=results.get('summary', {})
        )
        with open(report_path, "w") as f:
            f.write(html)
        return report_path

class AnalysisPipeline:
    """Main analysis pipeline orchestrator."""
    def __init__(self, args):
        self.args = args
        self.state = PipelineState(args)
        self.results = {
            "summary": {},
            "alignment": None,
            "tree_newick": None,
            "variants": {},
            "motifs": "",
            "fastqc_report": "",
            "primers": []
        }

    def run(self):
        """Execute pipeline with checkpoint recovery."""
        try:
            if self.state.state['stage'] == 'complete' and not self.args.force:
                logger.info("Previous complete run found. Use --force to re-run.")
                return
            self._validate_environment()
            for query in self.args.queries:
                query_seq = self._load_query(query)
                if self.args.validate_sequence:
                    self._validate_sequence(query_seq)
                    self.state.save('validation', {'fastqc_report': self.results['fastqc_report']})
                if self.args.design_primers:
                    self._design_primers(query_seq)
                    self.state.save('primers', self.results)
                blast_hits = self._run_blast(query_seq)
                sequences = self._process_hits(blast_hits)
                self.state.save('fetch', {'sequences': len(sequences)})
                if self.args.program in ["blastx", "blastp"]:
                    sequences = SequenceAnalyzer.translate_sequences(sequences)
                alignment = self._align_sequences(sequences)
                self.state.save('align', {'alignment': str(alignment)})
                self._analyze_alignment(alignment)
                if self.args.call_variants:
                    self._call_variants(alignment)
                    self.state.save('variants', {'variants': self.results['variants']})
                if self.args.find_motifs:
                    self._find_motifs(alignment)
                    self.state.save('motifs', {'motifs': self.results['motifs']})
            self._generate_outputs()
            self._cleanup()
            self.state.save("complete", {})
            logger.info("Pipeline completed successfully")
        except Exception as e:
            logger.critical(f"Pipeline failed: {str(e)}")
            raise

    def _validate_environment(self):
        """Validate dependencies and environment setup."""
        Path(self.args.output).mkdir(exist_ok=True)
        NcbiHandler.set_email(self.args.email)
        required = {
            'muscle': 'MUSCLE required for alignment',
            'samtools': 'samtools required for BAM file generation',
            'minimap2': 'minimap2 required for mapping to reference'
        }
        if self.args.local_blast:
            required['blastn'] = 'BLAST+ required for local searches'
            required['blastx'] = 'BLAST+ required for local BLASTX searches'
            required['blastp'] = 'BLAST+ required for local BLASTP searches'
        if self.args.call_variants:
            required['bcftools'] = 'bcftools required for variant calling'
        if self.args.find_motifs:
            required['meme'] = 'MEME required for motif analysis'
        if self.args.validate_sequence:
            required['fastqc'] = 'FastQC required for sequence validation'
        if self.args.advanced_tree:
            required['iqtree'] = 'IQ-TREE required for advanced phylogenetic analysis'
        missing = [cmd for cmd, msg in required.items() if not shutil.which(cmd)]
        if missing:
            raise EnvironmentError(f"Missing dependencies: {', '.join(missing)}. Please install them and try again.")

    def _load_query(self, query):
        """Load and validate input query sequence."""
        if Path(query).exists():
            try:
                record = SeqIO.read(Path(query), "fasta")
                self.results['summary']['query'] = {
                    "id": record.id,
                    "length": len(record),
                    "description": record.description
                }
                return str(record.seq)
            except Exception as e:
                raise ValueError(f"Invalid FASTA file: {e}")
        elif query.startswith("https://") or query.startswith("ftp://"):
            logger.info(f"Downloading query sequence from URL: {query}")
            response = subprocess.run(["curl", "-sL", query], capture_output=True, text=True)
            if response.returncode != 0:
                raise ValueError(f"Failed to download sequence from URL: {response.stderr}")
            try:
                record = SeqIO.read(StringIO(response.stdout), "fasta")
                self.results['summary']['query'] = {
                    "id": record.id,
                    "length": len(record),
                    "description": record.description
                }
                return str(record.seq)
            except Exception as e:
                raise ValueError(f"Invalid sequence downloaded from URL: {e}")
        else:
            logger.info(f"Fetching query sequence by accession number: {query}")
            try:
                handle = Entrez.efetch(db="nucleotide", id=query, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                self.results['summary']['query'] = {
                    "id": record.id,
                    "length": len(record),
                    "description": record.description
                }
                return str(record.seq)
            except Exception as e:
                raise ValueError(f"Failed to fetch sequence by accession number: {e}")

    def _run_blast(self, query_seq):
        """Execute BLAST with checkpoint recovery."""
        if self.state.state['stage'] in ['blast', 'fetch', 'align', 'tree', 'variants', 'motifs', 'validation']:
            logger.info("Resuming from BLAST results")
            return self.state.state['data'].get('blast_hits', [])
        blast_results = BlastRunner.run(query_seq, self.args)
        parsed_hits = self._parse_blast(blast_results)
        if not parsed_hits:
            raise ValueError("No valid BLAST hits found")
        self.state.save("blast", {"blast_hits": parsed_hits})
        return parsed_hits

    def _parse_blast(self, blast_handle):
        """Parse BLAST results with enhanced filtering."""
        blast_records = NCBIXML.parse(blast_handle)
        hits = []
        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < self.args.e_value and hsp.score > self.args.score_threshold:
                        hits.append({
                            'id': alignment.title.split()[1],
                            'e_value': hsp.expect,
                            'score': hsp.score,
                            'alignment_length': hsp.align_length
                        })
        return hits

    def _process_hits(self, hits):
        """Process BLAST hits with validation."""
        hit_ids = [hit['id'] for hit in hits[:self.args.max_hits]]
        if not hit_ids:
            raise ValueError("No valid BLAST hits to process")
        return SequenceAnalyzer.fetch_sequences(hit_ids, self.args)

    def _align_sequences(self, sequences):
        """Perform sequence alignment with resume support."""
        if self.state.state['stage'] in ['align', 'tree', 'variants', 'motifs', 'validation']:
            return AlignIO.read(Path(self.args.output) / "alignment.aln", "fasta")
        return SequenceAnalyzer.align(sequences, self.args)

    def _analyze_alignment(self, alignment):
        """Analyze alignment to build phylogenetic tree and identify variants."""
        logger.info("Analyzing alignment")
        if self.args.advanced_tree:
            tree_newick = PhylogeneticAnalyzer.construct_tree(Path(self.args.output) / "alignment.aln", self.args)
        else:
            constructor = DistanceTreeConstructor()
            matrix = constructor._calculate_distance(alignment)
            tree = constructor.upgma(matrix)
            tree_newick = tree.format('newick')
        self.results['tree_newick'] = tree_newick

    def _call_variants(self, alignment):
        """Call variants using bcftools."""
        if not self.args.reference:
            raise ValueError("Reference genome required for variant calling")
        reference = self.args.reference
        reference_file = Path(self.args.output) / "reference.fasta"
        if Path(reference).exists():
            reference_file = Path(reference)
        elif reference.startswith("https://") or reference.startswith("ftp://"):
            logger.info(f"Downloading reference genome from URL: {reference}")
            response = subprocess.run(["curl", "-sL", reference], capture_output=True, text=True)
            if response.returncode != 0:
                raise ValueError(f"Failed to download reference genome from URL: {response.stderr}")
            with open(reference_file, "w") as f:
                f.write(response.stdout)
        else:
            logger.info(f"Fetching reference genome by accession number: {reference}")
            try:
                handle = Entrez.efetch(db="nucleotide", id=reference, rettype="fasta", retmode="text")
                with open(reference_file, "w") as f:
                    f.write(handle.read())
            except Exception as e:
                raise ValueError(f"Failed to fetch reference genome by accession number: {e}")
        if not reference_file.exists():
            raise FileNotFoundError(f"Reference file {reference_file} not found")
        fasta_file = Path(self.args.output) / "alignment.fasta"
        AlignIO.write(alignment, fasta_file, "fasta")
        bam_file = Path(self.args.output) / "alignment.bam"
        sam_cmd = ["minimap2", "-a", str(reference_file), str(fasta_file)]
        view_cmd = ["samtools", "view", "-bS", "-o", str(bam_file)]
        logger.info(f"Generating BAM file: minimap2 -> samtools")
        try:
            sam = subprocess.run(sam_cmd, capture_output=True, check=True)
            subprocess.run(view_cmd, input=sam.stdout, check=True)
        except subprocess.CalledProcessError as e:
            logger.error(f"BAM generation failed: {e}")
            raise RuntimeError("BAM file generation failed.")
        self.results['variants'] = VariantCaller.call_variants(reference_file, bam_file, self.args)

    def _find_motifs(self, alignment):
        """Find motifs using MEME."""
        fasta_file = Path(self.args.output) / "alignment.fasta"
        AlignIO.write(alignment, fasta_file, "fasta")
        self.results['motifs'] = MotifAnalyzer.find_motifs(fasta_file, self.args)

    def _validate_sequence(self, query_seq):
        """Validate sequence using FastQC."""
        fasta_file = Path(self.args.output) / "query.fasta"
        with open(fasta_file, "w") as f:
            SeqIO.write(SeqIO.SeqRecord(Seq(query_seq), id="query", description=""), f, "fasta")
        self.results['fastqc_report'] = SequenceValidator.validate_sequence(fasta_file, self.args)

    def _design_primers(self, query_seq):
        """Design primers for the query sequence."""
        logger.info("Designing primers")
        try:
            primer_results = primer3.bindings.designPrimers(
                {
                    'SEQUENCE_ID': 'Query',
                    'SEQUENCE_TEMPLATE': query_seq,
                    'PRIMER_OPT_SIZE': 20,
                    'PRIMER_MIN_SIZE': 18,
                    'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': 60.0,
                    'PRIMER_MIN_TM': 57.0,
                    'PRIMER_MAX_TM': 63.0,
                    'PRIMER_MAX_GC': 50.0,
                    'PRIMER_MIN_GC': 20.0,
                    'PRIMER_MAX_LIBRARY_MISPRIMING': 8.0,
                    'PRIMER_PAIR_MAX_LIBRARY_MISPRIMING': 8.0,
                    'PRIMER_PRODUCT_SIZE_RANGE': [[100, 300]]
                },
                {}
            )
            self.results['primers'] = [
                {
                    'forward': primer_results['PRIMER_LEFT_0_SEQUENCE'],
                    'reverse': primer_results['PRIMER_RIGHT_0_SEQUENCE']
                }
            ]
        except Exception as e:
            logger.error(f"Primer design failed: {e}")
            self.results['primers'] = []

    def _generate_outputs(self):
        """Generate additional output files."""
        logger.info("Generating output files")
        AlignIO.write(self.results['alignment'], Path(self.args.output) / "alignment.aln", "clustal")
        AlignIO.write(self.results['alignment'], Path(self.args.output) / "alignment.fasta", "fasta")
        with open(Path(self.args.output) / "tree.nwk", "w") as f:
            f.write(self.results['tree_newick'])
        with open(Path(self.args.output) / "primers.json", "w") as f:
            json.dump(self.results['primers'], f, indent=2)
        if self.results.get('variants'):
            with open(Path(self.args.output) / "variants.vcf", "w") as f:
                f.write(str(self.results['variants']))
        if self.results.get('motifs'):
            with open(Path(self.args.output) / "motifs.txt", "w") as f:
                f.write(self.results['motifs'])

    def _cleanup(self):
        """Cleanup temporary files."""
        if not self.args.keep_temp:
            cache_dir = Path(self.args.output) / BlastRunner.CACHE_DIR
            if cache_dir.exists():
                shutil.rmtree(cache_dir)

def main():
    """Command-line interface and main execution flow."""
    parser = argparse.ArgumentParser(
        description="Enhanced Genomic Analysis Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show program's version number and exit"
    )
    parser.add_argument(
        "queries",
        nargs="+",
        help="One or more input sequences (FASTA file path, accession number, or URL)"
    )
    parser.add_argument("-o", "--output", default="results", help="Output directory")
    parser.add_argument("--config", help="Path to YAML/JSON config file")
    parser.add_argument("--local_blast", help="Path to local BLAST database")
    parser.add_argument("--program", choices=["blastn", "blastx", "blastp"], default="blastn", help="BLAST program")
    parser.add_argument("--database", default="nt", help="BLAST database")
    parser.add_argument("--max_hits", type=int, default=50, help="Maximum hits to process")
    parser.add_argument("--threads", type=int, default=4, help="Processing threads")
    parser.add_argument("--force", action="store_true", help="Overwrite existing results")
    parser.add_argument("--e_value", type=float, default=1e-5, help="E-value threshold for BLAST")
    parser.add_argument("--score_threshold", type=int, default=50, help="Score threshold for BLAST")
    parser.add_argument("--keep_temp", action="store_true", help="Keep temporary files")
    parser.add_argument("--resume", action="store_true", help="Resume from checkpoint")
    parser.add_argument("--design_primers", action="store_true", help="Design primers for the query sequence")
    parser.add_argument("--log_level", default="INFO", help="Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)")
    parser.add_argument("--email", help="NCBI email address")
    parser.add_argument("--call_variants", action="store_true", help="Call variants using bcftools")
    parser.add_argument("--reference", help="Reference genome (FASTA file path, accession number, or URL) for variant calling")
    parser.add_argument("--find_motifs", action="store_true", help="Find motifs using MEME")
    parser.add_argument("--num_motifs", type=int, default=5, help="Number of motifs to find")
    parser.add_argument("--min_motif_width", type=int, default=6, help="Minimum motif width")
    parser.add_argument("--max_motif_width", type=int, default=12, help="Maximum motif width")
    parser.add_argument("--advanced_tree", action="store_true", help="Use IQ-TREE for advanced phylogenetic analysis")
    parser.add_argument("--validate_sequence", action="store_true", help="Validate sequence using FastQC")
    parser.add_argument("--dry_run", action="store_true", help="Simulate pipeline execution without running commands.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging.")

    args = parser.parse_args()

    if args.config:
        args = ConfigManager.load(args.config, args)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    log_level = logging.DEBUG if args.verbose else getattr(logging, args.log_level.upper())
    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(
        level=log_level,
        format=log_format,
        handlers=[
            logging.FileHandler(output_dir / "pipeline.log"),
            logging.StreamHandler()
        ]
    )

    def _check_dependencies():
        required_tools = ['muscle', 'samtools', 'minimap2']
        if args.local_blast:
            required_tools.extend(['blastn', 'blastx', 'blastp'])
        if args.call_variants:
            required_tools.append('bcftools')
        if args.find_motifs:
            required_tools.append('meme')
        if args.validate_sequence:
            required_tools.append('fastqc')
        if args.advanced_tree:
            required_tools.append('iqtree')
        missing = [tool for tool in required_tools if not shutil.which(tool)]
        if missing:
            logger.critical(f"Missing dependencies: {', '.join(missing)}. Please install them and try again.")
            sys.exit(1)

    _check_dependencies()

    if args.dry_run:
        logger.info("Dry run mode: No commands will be executed.")
        return

    try:
        pipeline = AnalysisPipeline(args)
        pipeline.run()
    except Exception as e:
        logger.critical(f"Pipeline failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
