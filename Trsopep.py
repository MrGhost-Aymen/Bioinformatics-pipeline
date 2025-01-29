"""
Advanced Genomic Analysis Pipeline By Trso
This script performs comprehensive genomic analysis including BLAST queries, sequence alignment,
phylogenetic analysis, and primer design with enhanced configuration management and error handling.
Example:
    $ python script.py input.fasta --config config.yaml --output results
    Expected outputs:
    - results/report.html: Interactive analysis report
    - results/alignment.aln: Multiple sequence alignment
    - results/primers.json: Designed primer sequences
    - results/variants.vcf: SNP and indel calls
    - results/motifs.txt: Detected motifs
"""
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
from io import StringIO
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import yaml  # Requires PyYAML package
from Bio import SeqIO, Entrez, Phylo, AlignIO
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.Align.Applications import ClustalwCommandline
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from jinja2 import Template
import primer3
import logomaker
from Bio.Seq import Seq
from collections import defaultdict
import plotly.graph_objects as go
import subprocess
import fastqc
import meme

# Configure logging
logger = logging.getLogger(__name__)

class PipelineState:
    """Manages pipeline checkpoints and state recovery
    Attributes:
        state (dict): Current pipeline state including stage and data
        checkpoint_file (Path): Path to JSON checkpoint file
    """
    CHECKPOINT_STAGES = ['init', 'blast', 'fetch', 'align', 'tree', 'primers', 'variants', 'motifs', 'validation', 'complete']

    def __init__(self, args):
        self.args = args
        self.checkpoint_file = Path(args.output) / "checkpoint.json"
        self.state = self._load_checkpoint() or {
            "stage": "init",
            "data": {},
            "timestamp": datetime.now().isoformat()
        }

    def _load_checkpoint(self):
        """Load existing checkpoint if available and resume flag is set"""
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
        """Persist current state to checkpoint file"""
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
    """Handles configuration loading from file and command-line
    Args:
        config_path (str): Path to YAML/JSON config file
        args (Namespace): Command-line arguments namespace
    """
    @staticmethod
    def load(config_path, args):
        """Merge config file settings with command-line arguments"""
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
    """Manages NCBI Entrez interactions and credentials"""
    @staticmethod
    def set_email(email):
        """Validate and set NCBI email address"""
        Entrez.email = email or os.getenv("NCBI_EMAIL")
        if not Entrez.email:
            raise ValueError("NCBI email required. Set via --email or NCBI_EMAIL env var")
        logger.info(f"Using NCBI email: {Entrez.email}")

class BlastRunner:
    """Handles BLAST searches with local/remote execution and caching"""
    CACHE_DIR = "blast_cache"

    @classmethod
    def run(cls, query, args):
        """Execute BLAST with caching and resume support"""
        cache_key = cls._cache_key(query, args)
        cached = cls._check_cache(cache_key, args)
        if cached and not args.force:
            logger.info("Using cached BLAST results")
            return cached
        result = cls._run_remote(query, args) if not args.local_blast else cls._run_local(query, args)
        cls._save_cache(cache_key, result, args)
        return result

    @classmethod
    def _cache_key(cls, query, args):
        """Generate SHA256 hash for BLAST query caching"""
        key_data = {
            'query': query,
            'program': args.program,
            'database': args.database,
            'params': {
                'e_value': args.e_value,
                'perc_identity': args.perc_identity,
                'max_hits': args.max_hits
            }
        }
        return hashlib.sha256(json.dumps(key_data, sort_keys=True).encode()).hexdigest()

    @classmethod
    def _check_cache(cls, cache_key, args):
        """Check for existing cached BLAST results"""
        cache_dir = Path(args.output) / cls.CACHE_DIR
        cache_file = cache_dir / f"{cache_key}.xml"
        return cache_file if cache_file.exists() else None

    @classmethod
    def _save_cache(cls, cache_key, result, args):
        """Save BLAST results to cache"""
        cache_dir = Path(args.output) / cls.CACHE_DIR
        cache_dir.mkdir(exist_ok=True)
        cache_file = cache_dir / f"{cache_key}.xml"
        with open(cache_file, "w") as f:
            f.write(result.read())
        result.seek(0)
        return cache_file

    @staticmethod
    def _run_remote(query, args):
        """Execute remote BLAST via NCBI WWW service"""
        logger.info("Running remote BLAST via NCBI...")
        try:
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
        """Execute local BLAST against specified database"""
        logger.info(f"Running local BLAST against {args.local_blast}")
        output_file = Path(args.output) / "local_blast.xml"
        # Validate database path
        if not Path(args.local_blast).exists():
            raise FileNotFoundError(f"BLAST database {args.local_blast} not found")
        blastn_cline = NcbiblastnCommandline(
            query=query,
            db=args.local_blast,
            outfmt=5,
            out=output_file,
            num_threads=args.threads
        )
        logger.info(f"Executing: {blastn_cline}")
        stdout, stderr = blastn_cline()
        if stderr:
            logger.error(f"BLAST error: {stderr}")
            raise RuntimeError(stderr)
        return open(output_file)

class SequenceAnalyzer:
    """Handles sequence fetching, alignment, and analysis"""
    @staticmethod
    def fetch_sequences(hit_ids, args):
        """Fetch sequences from GenBank with parallel requests"""
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
        return results

    @staticmethod
    def _fetch_single(hit_id):
        """Fetch single sequence from GenBank"""
        try:
            handle = Entrez.efetch(db="nucleotide", id=hit_id, rettype="fasta", retmode="text")
            return SeqIO.read(handle, "fasta")
        except Exception as e:
            logger.error(f"Failed to fetch {hit_id}: {e}")
            return None

    @staticmethod
    def align(sequences, args):
        """Perform multiple sequence alignment with ClustalW"""
        if not sequences:
            raise ValueError("No sequences available for alignment")
        logger.info(f"Aligning {len(sequences)} sequences")
        input_file = Path(args.output) / "alignment_input.fasta"
        output_file = Path(args.output) / "alignment.aln"
        SeqIO.write(sequences, input_file, "fasta")
        clustalw_cline = ClustalwCommandline(
            "clustalw2",
            infile=input_file,
            outfile=output_file,
            type="dna",
            quicktree=True
        )
        logger.info(f"Executing: {clustalw_cline}")
        stdout, stderr = clustalw_cline()
        if not args.keep_temp:
            input_file.unlink(missing_ok=True)
        return AlignIO.read(output_file, "clustal")

class VariantCaller:
    """Handles variant calling using external tools like bcftools"""
    @staticmethod
    def call_variants(reference, bam_file, args):
        """Call variants using bcftools"""
        vcf_file = Path(args.output) / "variants.vcf"
        bcftools_cmd = f"bcftools mpileup -Ou -f {reference} {bam_file} | bcftools call -mv -Ov -o {vcf_file}"
        logger.info(f"Executing variant calling: {bcftools_cmd}")
        subprocess.run(bcftools_cmd, shell=True, check=True)
        return vcf_file

class MotifAnalyzer:
    """Handles motif analysis using external tools like MEME"""
    @staticmethod
    def find_motifs(fasta_file, args):
        """Find motifs using MEME"""
        meme_output_dir = Path(args.output) / "meme_output"
        meme_output_dir.mkdir(exist_ok=True)
        meme_cmd = f"meme {fasta_file} -oc {meme_output_dir} -nmotifs {args.num_motifs} -minw {args.min_motif_width} -maxw {args.max_motif_width}"
        logger.info(f"Executing motif analysis: {meme_cmd}")
        subprocess.run(meme_cmd, shell=True, check=True)
        motifs_file = meme_output_dir / "meme.txt"
        with open(motifs_file) as f:
            motifs = f.read()
        return motifs

class PhylogeneticAnalyzer:
    """Handles advanced phylogenetic analysis using IQ-TREE"""
    @staticmethod
    def construct_tree(alignment_file, args):
        """Construct phylogenetic tree using IQ-TREE"""
        tree_file = Path(args.output) / "tree.nwk"
        iqtree_cmd = f"iqtree -s {alignment_file} -st DNA -pre {tree_file.stem} -nt AUTO"
        logger.info(f"Executing phylogenetic tree construction: {iqtree_cmd}")
        subprocess.run(iqtree_cmd, shell=True, check=True)
        with open(tree_file) as f:
            tree_newick = f.read()
        return tree_newick

class SequenceValidator:
    """Handles sequence validation using FastQC"""
    @staticmethod
    def validate_sequence(fasta_file, args):
        """Validate sequence using FastQC"""
        fastqc_output_dir = Path(args.output) / "fastqc_output"
        fastqc_cmd = f"fastqc {fasta_file} -o {fastqc_output_dir} --extract"
        logger.info(f"Executing sequence validation: {fastqc_cmd}")
        subprocess.run(fastqc_cmd, shell=True, check=True)
        return fastqc_output_dir

class VisualizationEngine:
    """Generates interactive HTML reports with visualization components"""
    HTML_TEMPLATE = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Genomic Analysis Report</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    </head>
    <body>
        <h1>Genomic Analysis Report</h1>
        <h2>Phylogenetic Tree</h2>
        <div id="tree"></div>
        <h2>Sequence Conservation Logo</h2>
        <img src="data:image/png;base64,{{ logo_img }}" alt="Conservation Logo">
        <h2>Designed Primers</h2>
        <ul>
            {% for primer in primers %}
            <li>{{ primer }}</li>
            {% endfor %}
        </ul>
        <h2>Variants</h2>
        <pre>{{ variants | tojson(indent=2) }}</pre>
        <h2>Motifs</h2>
        <pre>{{ motifs }}</pre>
        <h2>FastQC Report</h2>
        <iframe src="{{ fastqc_report }}" width="100%" height="600px"></iframe>
        <h2>Summary</h2>
        <pre>{{ summary | tojson(indent=2) }}</pre>
        <script>
            {{ tree_html|safe }}
        </script>
    </body>
    </html>
    """

    @staticmethod
    def generate_report(results, args):
        """Generate comprehensive HTML report"""
        logger.info("Compiling final report")
        report_path = Path(args.output) / "report.html"
        # Generate visualization assets
        logo_img = VisualizationEngine._generate_conservation_logo(results['alignment'])
        tree_html = VisualizationEngine._generate_interactive_tree(results.get('tree_newick'))
        # Render template
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

    @staticmethod
    def _generate_conservation_logo(alignment):
        """Generate sequence logo using logomaker"""
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

    @staticmethod
    def _generate_interactive_tree(newick_str):
        """Generate Plotly-based phylogenetic tree visualization"""
        if not newick_str:
            return "Interactive tree placeholder"
        tree = Phylo.read(StringIO(newick_str), "newick")
        fig = go.Figure(data=[go.Scatter(
            x=[],
            y=[],
            mode='markers+lines',
            text=[],
            hoverinfo='text',
            marker=dict(color='blue', size=10)
        )])

        def draw_clade(clade, parent_x, parent_y, depth):
            x_values = [parent_x]
            y_values = [parent_y]
            clade_x = parent_x + 1
            clade_y = parent_y - depth
            x_values.append(clade_x)
            y_values.append(clade_y)
            fig.data[0].x += tuple(x_values)
            fig.data[0].y += tuple(y_values)
            fig.data[0].text += (clade.name,)
            if clade.clades:
                for subclade in clade.clades:
                    draw_clade(subclade, clade_x, clade_y, depth + 1)

        draw_clade(tree.root, 0, 0, 0)
        fig.update_layout(title='Phylogenetic Tree', showlegend=False)
        return fig.to_html(full_html=False)

class AnalysisPipeline:
    """Main analysis pipeline orchestrator"""
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
        """Execute pipeline with checkpoint recovery"""
        try:
            if self.state.state['stage'] == 'complete' and not self.args.force:
                logger.info("Previous complete run found. Use --force to re-run.")
                return
            self._validate_environment()
            query_seq = self._load_query()
            if self.args.validate_sequence:
                self._validate_sequence(query_seq)
                self.state.save('validation', {'fastqc_report': self.results['fastqc_report']})
            if self.args.design_primers:
                self._design_primers(query_seq)
                self.state.save('primers', self.results)
            blast_hits = self._run_blast(query_seq)
            sequences = self._process_hits(blast_hits)
            self.state.save('fetch', {'sequences': len(sequences)})
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
            logger.error(f"Pipeline failed: {e}")
            raise

    def _validate_environment(self):
        """Validate dependencies and environment setup"""
        Path(self.args.output).mkdir(exist_ok=True)
        NcbiHandler.set_email(self.args.email)
        # Check required executables
        required = {'clustalw2': 'ClustalW required for alignment'}
        if self.args.local_blast:
            required['blastn'] = 'BLAST+ required for local searches'
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
            raise EnvironmentError(f"Missing dependencies: {', '.join(missing)}")

    def _load_query(self):
        """Load and validate input query sequence"""
        query_path = Path(self.args.query)
        if query_path.exists():
            try:
                record = SeqIO.read(query_path, "fasta")
                self.results['summary']['query'] = {
                    "id": record.id,
                    "length": len(record),
                    "description": record.description
                }
                return str(record.seq)
            except Exception as e:
                raise ValueError(f"Invalid FASTA file: {e}")
        return self.args.query

    def _run_blast(self, query_seq):
        """Execute BLAST with checkpoint recovery"""
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
        """Parse BLAST results with enhanced filtering"""
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
        """Process BLAST hits with validation"""
        hit_ids = [hit['id'] for hit in hits[:self.args.max_hits]]
        if not hit_ids:
            raise ValueError("No valid BLAST hits to process")
        return SequenceAnalyzer.fetch_sequences(hit_ids, self.args)

    def _align_sequences(self, sequences):
        """Perform sequence alignment with resume support"""
        if self.state.state['stage'] in ['align', 'tree', 'variants', 'motifs', 'validation']:
            return AlignIO.read(Path(self.args.output) / "alignment.aln", "clustal")
        return SequenceAnalyzer.align(sequences, self.args)

    def _analyze_alignment(self, alignment):
        """Analyze alignment to build phylogenetic tree and identify variants"""
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
        """Call variants using bcftools"""
        if not self.args.reference:
            raise ValueError("Reference genome required for variant calling")
        reference_file = Path(self.args.reference)
        if not reference_file.exists():
            raise FileNotFoundError(f"Reference file {reference_file} not found")
        # Convert alignment to BAM format (this is a placeholder, actual conversion needed)
        bam_file = Path(self.args.output) / "alignment.bam"
        # Placeholder for BAM conversion logic
        self.results['variants'] = VariantCaller.call_variants(reference_file, bam_file, self.args)

    def _find_motifs(self, alignment):
        """Find motifs using MEME"""
        fasta_file = Path(self.args.output) / "alignment.fasta"
        AlignIO.write(alignment, fasta_file, "fasta")
        self.results['motifs'] = MotifAnalyzer.find_motifs(fasta_file, self.args)

    def _validate_sequence(self, query_seq):
        """Validate sequence using FastQC"""
        fasta_file = Path(self.args.output) / "query.fasta"
        with open(fasta_file, "w") as f:
            SeqIO.write(SeqIO.SeqRecord(Seq(query_seq), id="query", description=""), f, "fasta")
        self.results['fastqc_report'] = SequenceValidator.validate_sequence(fasta_file, self.args)

    def _design_primers(self, query_seq):
        """Design primers for the query sequence"""
        logger.info("Designing primers")
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

    def _generate_outputs(self):
        """Generate additional output files"""
        logger.info("Generating output files")
        # Save alignment in ALN format
        AlignIO.write(self.results['alignment'], Path(self.args.output) / "alignment.aln", "clustal")
        # Save alignment in FASTA format
        AlignIO.write(self.results['alignment'], Path(self.args.output) / "alignment.fasta", "fasta")
        # Save tree in NEWICK format
        with open(Path(self.args.output) / "tree.nwk", "w") as f:
            f.write(self.results['tree_newick'])
        # Save primers in JSON format
        with open(Path(self.args.output) / "primers.json", "w") as f:
            json.dump(self.results['primers'], f, indent=2)
        # Save variants in VCF format
        if self.results.get('variants'):
            with open(Path(self.args.output) / "variants.vcf", "w") as f:
                f.write(self.results['variants'])
        # Save motifs in TXT format
        if self.results.get('motifs'):
            with open(Path(self.args.output) / "motifs.txt", "w") as f:
                f.write(self.results['motifs'])

    def _cleanup(self):
        """Cleanup temporary files"""
        if not self.args.keep_temp:
            cache_dir = Path(self.args.output) / BlastRunner.CACHE_DIR
            if cache_dir.exists():
                shutil.rmtree(cache_dir)

def main():
    """Command-line interface and main execution flow"""
    parser = argparse.ArgumentParser(
        description="Enhanced Genomic Analysis Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # Core arguments
    parser.add_argument("query", help="Input sequence or FASTA file path")
    parser.add_argument("-o", "--output", default="results", help="Output directory")
    parser.add_argument("--config", help="Path to YAML/JSON config file")
    # BLAST parameters
    parser.add_argument("--local_blast", help="Path to local BLAST database")
    parser.add_argument("--program", default="blastn", help="BLAST program")
    parser.add_argument("--database", default="nt", help="BLAST database")
    # Advanced parameters
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
    # Variant calling
    parser.add_argument("--call_variants", action="store_true", help="Call variants using bcftools")
    parser.add_argument("--reference", help="Reference genome for variant calling")
    # Motif analysis
    parser.add_argument("--find_motifs", action="store_true", help="Find motifs using MEME")
    parser.add_argument("--num_motifs", type=int, default=5, help="Number of motifs to find")
    parser.add_argument("--min_motif_width", type=int, default=6, help="Minimum motif width")
    parser.add_argument("--max_motif_width", type=int, default=12, help="Maximum motif width")
    # Advanced phylogenetic analysis
    parser.add_argument("--advanced_tree", action="store_true", help="Use IQ-TREE for advanced phylogenetic analysis")
    # Sequence validation
    parser.add_argument("--validate_sequence", action="store_true", help="Validate sequence using FastQC")
    args = parser.parse_args()
    # Load configuration file if specified
    if args.config:
        args = ConfigManager.load(args.config, args)
    # Configure logging
    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format=log_format,
        handlers=[
            logging.FileHandler(Path(args.output) / "pipeline.log"),
            logging.StreamHandler()
        ]
    )
    try:
        pipeline = AnalysisPipeline(args)
        pipeline.run()
    except Exception as e:
        logger.critical(f"Pipeline failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
