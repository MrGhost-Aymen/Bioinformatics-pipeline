import argparse
import base64
import json
import logging
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO, Entrez, Phylo, AlignIO
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.Align.Applications import ClustalwCommandline
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from jinja2 import Template
import primer3
import RNA
import logomaker
from Bio.Seq import Seq
from collections import defaultdict

# Configure logging
logger = logging.getLogger(__name__)

class PipelineState:
    def __init__(self, args):
        self.args = args
        self.checkpoint_file = Path(args.output) / "checkpoint.json"
        self.state = self._load_checkpoint() or {
            "stage": "init",
            "data": {},
            "timestamp": datetime.now().isoformat()
        }

    def _load_checkpoint(self):
        if self.args.resume and self.checkpoint_file.exists():
            try:
                with open(self.checkpoint_file) as f:
                    return json.load(f)
            except Exception as e:
                logger.warning(f"Failed to load checkpoint: {e}")
        return None

    def save(self, stage, data):
        self.state.update({
            "stage": stage,
            "data": data,
            "timestamp": datetime.now().isoformat()
        })
        with open(self.checkpoint_file, "w") as f:
            json.dump(self.state, f)

class NcbiHandler:
    @staticmethod
    def set_email(email):
        Entrez.email = email or os.getenv("NCBI_EMAIL") or input("Enter NCBI email: ")
        logger.info(f"Using NCBI email: {Entrez.email}")

class BlastRunner:
    @staticmethod
    def run(query, args):
        if args.local_blast:
            return BlastRunner._run_local(query, args)
        return BlastRunner._run_remote(query, args)

    @staticmethod
    def _run_remote(query, args):
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
            logger.error(f"BLAST failed: {e}")
            raise

    @staticmethod
    def _run_local(query, args):
        logger.info(f"Running local BLAST against {args.local_blast}")
        output_file = Path(args.output) / "local_blast.xml"
        blastn_cline = NcbiblastnCommandline(
            query=query,
            db=args.local_blast,
            outfmt=5,
            out=output_file,
            num_threads=args.threads
        )
        logger.info(f"Executing: {blastn_cline}")
        blastn_cline()
        return open(output_file)

class SequenceAnalyzer:
    @staticmethod
    def fetch_sequences(hit_ids, args):
        logger.info(f"Fetching {len(hit_ids)} sequences from GenBank")
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            futures = [executor.submit(SequenceAnalyzer._fetch_single, hid) for hid in hit_ids]
            return [f.result() for f in futures if f.result()]

    @staticmethod
    def _fetch_single(hit_id):
        try:
            handle = Entrez.efetch(db="nucleotide", id=hit_id, rettype="fasta", retmode="text")
            return SeqIO.read(handle, "fasta")
        except Exception as e:
            logger.warning(f"Failed to fetch {hit_id}: {e}")
            return None

    @staticmethod
    def align(sequences, args):
        logger.info("Performing multiple sequence alignment")
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
        clustalw_cline()
        
        if not args.keep_temp:
            input_file.unlink(missing_ok=True)
        
        return AlignIO.read(output_file, "clustal")

class VisualizationEngine:
    HTML_TEMPLATE = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Analysis Report</title>
        <script src="https://d3js.org/d3.v7.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/d3-phylogram/2.0.1/d3-phylogram.min.js"></script>
        <style>
            .node circle { fill: #999; }
            .node text { font: 10px sans-serif; }
            .link { fill: none; stroke: #333; stroke-width: 1.5px; }
        </style>
    </head>
    <body>
        <h1>Sequence Analysis Report</h1>
        {% if tree_newick %}
        <div id="phylogram"></div>
        {% endif %}
        {% if logo_img %}<img src="data:image/png;base64,{{ logo_img }}">{% endif %}
        {% if primers %}
        <h2>Designed Primers</h2>
        {% for pair in primers %}
        <div>
            <h3>Primer Pair {{ loop.index }}</h3>
            <p>Left: {{ pair.left.sequence }} ({{ pair.left.position[0] }}-{{ pair.left.position[0] + pair.left.position[1] }})</p>
            <p>Right: {{ pair.right.sequence }} ({{ pair.right.position[0] }}-{{ pair.right.position[0] + pair.right.position[1] }})</p>
            <p>Product Size: {{ pair.product_size }}bp</p>
        </div>
        {% endfor %}
        {% endif %}
        <script>
            {% if tree_newick %}
            (function() {
                const newick = '{{ tree_newick }}';
                const width = 1000;
                const height = 800;
                
                d3.phylogeny()
                    .size([width, height])
                    .newick(newick)
                    .labels(true)
                    .render('#phylogram');
            })();
            {% endif %}
        </script>
    </body>
    </html>
    """

    @staticmethod
    def generate_report(results, args):
        logger.info("Generating HTML report")
        report_path = Path(args.output) / "report.html"
        
        logo_img = VisualizationEngine._generate_logo(results['alignment'], args)
        tree_newick = results.get('tree_newick')
        primers = results.get('primers', [])
        
        template = Template(VisualizationEngine.HTML_TEMPLATE)
        html = template.render(
            tree_newick=tree_newick,
            logo_img=logo_img,
            primers=primers,
            summary=json.dumps(results['summary'], indent=2)
        )
        
        with open(report_path, "w") as f:
            f.write(html)
        return report_path

    @staticmethod
    def _generate_logo(alignment, args):
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

class AnalysisPipeline:
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
        try:
            if self.state.state['stage'] == 'complete':
                logger.info("Previous complete run found. Use --force to re-run.")
                return

            self._setup_environment()
            query_seq = self._load_query()
            if self.args.design_primers:
                self._design_primers(query_seq)
            blast_results = self._run_blast(query_seq)
            sequences = self._process_hits(blast_results)
            alignment = self._align_sequences(sequences)
            self._analyze_alignment(alignment)
            self._generate_outputs()
            self._cleanup()
            self.state.save("complete", {})

        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise

    def _setup_environment(self):
        Path(self.args.output).mkdir(exist_ok=True)
        NcbiHandler.set_email(self.args.email)

    def _load_query(self):
        if Path(self.args.query).exists():
            record = SeqIO.read(self.args.query, "fasta")
            self.results['summary']['query'] = {
                "id": record.id,
                "length": len(record),
                "description": record.description,
                "sequence": str(record.seq)
            }
            return str(record.seq)
        return self.args.query

    def _run_blast(self, query_seq):
        if self.state.state['stage'] == 'blast':
            logger.info("Resuming from BLAST stage")
            return self.state.state['data'].get('blast_results')

        blast_results = BlastRunner.run(query_seq, self.args)
        parsed_hits = self._parse_blast(blast_results)
        self.state.save("blast", {"blast_results": parsed_hits})
        return parsed_hits

    def _parse_blast(self, blast_handle):
        blast_records = NCBIXML.parse(blast_handle)
        hits = []
        for record in blast_records:
            for alignment in record.alignments:
                best_hsp = None
                for hsp in alignment.hsps:
                    if best_hsp is None or hsp.expect < best_hsp.expect:
                        best_hsp = hsp
                if not best_hsp:
                    continue
                
                perc_ident = (best_hsp.identities / best_hsp.align_length) * 100
                if (best_hsp.expect > self.args.e_value or 
                    perc_ident < self.args.perc_identity):
                    continue

                hits.append({
                    "id": alignment.hit_id,
                    "accession": alignment.accession,
                    "e_value": best_hsp.expect,
                    "perc_identity": perc_ident,
                    "length": alignment.length
                })
                if len(hits) >= self.args.max_hits:
                    break
        return hits

    def _process_hits(self, hits):
        hit_ids = [hit['id'] for hit in hits[:self.args.max_hits]]
        sequences = SequenceAnalyzer.fetch_sequences(hit_ids, self.args)
        return [s for s in sequences if s is not None]

    def _align_sequences(self, sequences):
        if self.state.state['stage'] == 'alignment':
            return AlignIO.read(Path(self.args.output) / "alignment.aln", "clustal")
        return SequenceAnalyzer.align(sequences, self.args)

    def _analyze_alignment(self, alignment):
        self.results['alignment'] = alignment
        self._find_variants(alignment)
        if self.args.tree:
            self._build_tree(alignment)

    def _find_variants(self, alignment):
        summary = AlignInfo.SummaryInfo(alignment)
        consensus = summary.dumb_consensus()
        self.results['variants'] = []
        for i in range(len(consensus)):
            variants = defaultdict(int)
            for rec in alignment:
                variants[str(rec.seq[i])] += 1
            if len(variants) > 1:
                self.results['variants'].append({
                    "position": i+1,
                    "consensus": str(consensus[i]),
                    "variants": dict(variants)
                })

    def _build_tree(self, alignment):
        logger.info("Building phylogenetic tree")
        constructor = DistanceTreeConstructor()
        dm = constructor.get_distance(alignment)
        
        if self.args.tree_method == 'nj':
            tree = constructor.nj(dm)
        else:
            tree = constructor.upgma(dm)
            
        self.results['tree_newick'] = Phylo.to_newick(tree)
        Phylo.write(tree, Path(self.args.output)/"tree.nwk", "newick")

    def _design_primers(self, query_seq):
        logger.info("Designing primers")
        try:
            primer_params = {
                'SEQUENCE_ID': 'query',
                'SEQUENCE_TEMPLATE': query_seq,
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_MIN_SIZE': 18,
                'PRIMER_MAX_SIZE': 25,
                'PRIMER_NUM_RETURN': 3,
                'PRIMER_PRODUCT_SIZE_RANGE': [[100, 300]],
            }
            results = primer3.bindings.design_primers(primer_params)
            
            primers = []
            num_pairs = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
            for i in range(num_pairs):
                primers.append({
                    'left': {
                        'sequence': results[f'PRIMER_LEFT_{i}_SEQUENCE'],
                        'position': results[f'PRIMER_LEFT_{i}']
                    },
                    'right': {
                        'sequence': results[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                        'position': results[f'PRIMER_RIGHT_{i}']
                    },
                    'product_size': results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
                })
            self.results['primers'] = primers
        except Exception as e:
            logger.error(f"Primer design failed: {e}")

    def _generate_outputs(self):
        VisualizationEngine.generate_report(self.results, self.args)
        with open(Path(self.args.output)/"results.json", "w") as f:
            json.dump(self.results, f, indent=2)

    def _cleanup(self):
        if not self.args.keep_temp:
            temp_files = ["alignment_input.fasta", "temp_tree.png"]
            for f in temp_files:
                Path(f).unlink(missing_ok=True)

def main():
    parser = argparse.ArgumentParser(description="Advanced Genomic Analysis Pipeline")
    parser.add_argument("query", help="Input sequence or FASTA file")
    parser.add_argument("-o", "--output", default="results", help="Output directory")
    parser.add_argument("--local_blast", help="Path to local BLAST database")
    parser.add_argument("--max_hits", type=int, default=50, help="Maximum hits to process")
    parser.add_argument("--threads", type=int, default=4, help="Processing threads")
    parser.add_argument("--email", help="NCBI email address")
    parser.add_argument("--tree", action="store_true", help="Build phylogenetic tree")
    parser.add_argument("--tree_method", choices=['nj', 'upgma'], default='nj', 
                      help="Tree building method")
    parser.add_argument("--e_value", type=float, default=0.001, 
                      help="E-value threshold for BLAST hits")
    parser.add_argument("--perc_identity", type=float, default=30.0,
                      help="Minimum percent identity for BLAST hits")
    parser.add_argument("--design_primers", action="store_true",
                      help="Design primers for query sequence")
    parser.add_argument("--resume", action="store_true", help="Resume failed run")
    parser.add_argument("--keep_temp", action="store_true", help="Keep temporary files")
    parser.add_argument("--log_level", default="INFO", 
                      choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                      help="Logging verbosity level")
    
    args = parser.parse_args()

    # Configure logging
    Path(args.output).mkdir(exist_ok=True)
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(Path(args.output)/"pipeline.log"),
            logging.StreamHandler()
        ]
    )

    try:
        pipeline = AnalysisPipeline(args)
        pipeline.run()
        logger.info("Analysis completed successfully")
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
