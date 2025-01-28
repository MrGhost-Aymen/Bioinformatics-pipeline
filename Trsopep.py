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

# Configure logging
logger = logging.getLogger(__name__)

class PipelineState:
    """Manages pipeline checkpoints and state recovery
    
    Attributes:
        state (dict): Current pipeline state including stage and data
        checkpoint_file (Path): Path to JSON checkpoint file
    """
    
    CHECKPOINT_STAGES = ['init', 'blast', 'fetch', 'align', 'tree', 'primers', 'complete']

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

class VisualizationEngine:
    """Generates interactive HTML reports with visualization components"""
    
    HTML_TEMPLATE = """
    <!-- Interactive report template with enhanced visualization -->
    ... (extended template with Plotly integration)
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
        # Implementation using Plotly would go here
        return "<div>Interactive tree placeholder</div>"

class AnalysisPipeline:
    """Main analysis pipeline orchestrator"""
    
    def __init__(self, args):
        self.args = args
        self.state = PipelineState(args)
        self.results = {
            "summary": {},
            "alignment": None,
            "tree_newick": None,
            "variants": [],
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
            
            if self.args.design_primers:
                self._design_primers(query_seq)
                self.state.save('primers', self.results)
                
            blast_hits = self._run_blast(query_seq)
            sequences = self._process_hits(blast_hits)
            self.state.save('fetch', {'sequences': len(sequences)})
            
            alignment = self._align_sequences(sequences)
            self.state.save('align', {'alignment': str(alignment)})
            
            self._analyze_alignment(alignment)
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
        if self.state.state['stage'] in ['blast', 'fetch', 'align', 'tree']:
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
        ... (existing parsing logic with additional filters)

    def _process_hits(self, hits):
        """Process BLAST hits with validation"""
        hit_ids = [hit['id'] for hit in hits[:self.args.max_hits]]
        if not hit_ids:
            raise ValueError("No valid BLAST hits to process")
        return SequenceAnalyzer.fetch_sequences(hit_ids, self.args)

    def _align_sequences(self, sequences):
        """Perform sequence alignment with resume support"""
        if self.state.state['stage'] in ['align', 'tree']:
            return AlignIO.read(Path(self.args.output) / "alignment.aln", "clustal")
        return SequenceAnalyzer.align(sequences, self.args)

    ... (rest of existing methods with enhanced docstrings and error handling)

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
    
    # Add other parameters from original script
    ... (remaining arguments from original implementation)

    args = parser.parse_args()
    
    # Load configuration file if specified
    if args.config:
        args = ConfigManager.load(args.config, args)

    # Configure logging
    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(
        level=args.log_level,
        format=log_format,
        handlers=[
            logging.FileHandler(Path(args.output)/"pipeline.log"),
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
