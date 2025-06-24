"""
Author: Matthew Schmitz, Allen Institute, 2024

Append mitochondrial annotations (via MitoFinder) to a genome assembly.

usage:
    python add_mito.py \
        -a NC_012920.1 \
        -r /path/to/mito_ref.gbk \
        -g /path/to/genome.fa \
        -t /path/to/annotations.gtf \
        -s /path/to/sif_directory \
        -n rhemac10 \
        [-e your.email@example.com] \
        [--api_key YOUR_NCBI_API_KEY]
"""

import sys
import os
import argparse

def download_fasta_http(accession: str,
                        out_path: str = None,
                        email: str = "your.email@example.com",
                        api_key: str = None) -> str:
    """
    Fetches the FASTA via HTTP GET to NCBI E-utilities.
    """
    
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nucleotide",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        "email": email,
    }
    if api_key:
        params["api_key"] = api_key
    r = requests.get(base, params=params)
    r.raise_for_status()
    if out_path is None:
        out_path = accession + ".fasta"
    with open(out_path, "w") as f:
        f.write(r.text)
    return out_path

def gff_to_gtf(gff_path: str, gtf_path: str):
    """
    Convert a GFF from MitoFinder into a GTF, adding gene_id, gene, gene_name
    all set to the value of the Name= tag.
    """
    with open(gff_path) as gf, open(gtf_path, "w") as out:
        for line in gf:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = cols[8]
            # split on whitespace or semicolon to find Name=
            tokens = re.split(r"[;\s]+", attrs)
            name = None
            for tok in tokens:
                if tok.startswith("Name="):
                    name = tok.split("=", 1)[1]
                    break
            if not name:
                name = "."
            # build GTF‐style attr field
            cols[8] = f'gene_id "{name}"; gene "{name}"; gene_name "{name}";'
            out.write("\t".join(cols) + "\n")

def run_mitofinder_and_append(
    mito_accession: str,
    ref_gbk: str,
    genomic_fasta: str,
    transcript_gtf: str,
    sif_dir: str,
    spec_name: str,
    email: str = None,
    api_key: str = None
) -> (str, str):
    """
    1) Download mitochondrial FASTA by accession.
    2) Run MitoFinder (singularity).
    3) Convert its GFF → GTF, append to transcript_gtf → *_mito.gtf.
    4) Append mito.fa to genomic_fasta → *_mito.fasta.
    Returns (new_genome_fasta, new_transcript_gtf).
    """
    # --- 1) download mito fasta
    mito_fasta = download_fasta_http(
        mito_accession,
        out_path=f"{mito_accession}.fasta",
        email=email,
        api_key=api_key
    )

    # --- 2) run MitoFinder
    sif = os.path.join(sif_dir, "remiallio_default_mitofinder.sif")
    cmd = [
        "singularity", "run", sif,
        "-j", spec_name,
        "-a", mito_fasta,
        "-r", ref_gbk,
        "-o", "2",
        "--blast-identity-prot", "30",
        "--blast-size", "20",
        "--blast-identity-nucl", "30"
    ]
    subprocess.run(cmd, check=True)

    # --- 3) locate and convert its GFF
    results_dir = os.path.join(
        spec_name, f"{spec_name}_MitoFinder_mitfi_Final_Results"
    )
    gff_file = os.path.join(results_dir, f"{spec_name}_mtDNA_contig.gff")
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"Expected {gff_file}")

    mito_gtf = os.path.join(results_dir, f"{spec_name}_mtDNA_contig.gtf")
    gff_to_gtf(gff_file, mito_gtf)

    # --- 4) write new transcriptome GTF
    base_gtf = os.path.splitext(transcript_gtf)[0]
    out_gtf = f"{base_gtf}_mito.gtf"
    with open(out_gtf, "w") as w:
        # original annotations
        with open(transcript_gtf) as rf:
            w.write(rf.read().rstrip("\n") + "\n")
        # mito annotations
        with open(mito_gtf) as mf:
            w.write(mf.read())

    # --- 5) write new genomic FASTA
    base_fa = os.path.splitext(genomic_fasta)[0]
    out_fa = f"{base_fa}_mito.fasta"
    with open(out_fa, "w") as w:
        # original genome
        with open(genomic_fasta) as gf:
            w.write(gf.read().rstrip("\n") + "\n")
        # append mito-chromosome
        with open(mito_fasta) as mf:
            w.write(mf.read())

    return out_fa, out_gtf

def main():
    parser = argparse.ArgumentParser(
        description="Download mito‐chromosome, run MitoFinder, and append "
                    "mitochondrial annotations to a genome FASTA and GTF."
    )
    parser.add_argument(
        '-a', '--accession',
        required=True,
        help="NCBI GenBank accession for the mitochondrial chromosome (e.g. NC_012920.1)"
    )
    parser.add_argument(
        '-r', '--ref_gbk',
        type=os.path.abspath,
        required=True,
        help="Path to the reference GenBank .gbk file for the mitochondrion"
    )
    parser.add_argument(
        '-g', '--genome_fasta',
        type=os.path.abspath,
        required=True,
        help="Path to the genomic FASTA to which the mito sequence will be appended"
    )
    parser.add_argument(
        '-t', '--transcript_gtf',
        type=os.path.abspath,
        required=True,
        help="Path to the transcriptome GTF to which the mito GTF will be appended"
    )
    parser.add_argument(
        '-s', '--sif_dir',
        type=os.path.abspath,
        required=True,
        help="Directory containing remiallio_default_mitofinder.sif"
    )
    parser.add_argument(
        '-n', '--spec_name',
        required=True,
        help="Sample/species name (used by MitoFinder to name its output folder)"
    )
    parser.add_argument(
        '-e', '--email',
        help="NCBI Entrez email address (for EFetch)"
    )
    parser.add_argument(
        '-k', '--api_key',
        help="NCBI Entrez API key (optional, raises rate limits)"
    )

    args = parser.parse_args()

    print("Running MitoFinder and appending mito sequence/annotations…")
    try:
        new_fasta, new_gtf = run_mitofinder_and_append(
            mito_accession=args.accession,
            ref_gbk=args.ref_gbk,
            genomic_fasta=args.genome_fasta,
            transcript_gtf=args.transcript_gtf,
            sif_dir=args.sif_dir,
            spec_name=args.spec_name,
            email=args.email,
            api_key=args.api_key
        )
    except Exception as e:
        print(f"✖ Error: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Combined FASTA written to: {new_fasta}")
    print(f"Combined GTF written to:   {new_gtf}")

if __name__ == "__main__":
    main()








 