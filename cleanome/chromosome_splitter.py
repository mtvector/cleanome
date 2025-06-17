"""
Author: Matthew Schmitz, Allen Institute, 2024
"""
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv
import pyBigWig

def zcumsum(iterable):
    cumulative_sum = [0]
    total = 0
    for item in iterable:
        total += item
        cumulative_sum.append(total)
    cumulative_sum.append(float("inf"))
    return cumulative_sum


class ChromosomeSplitter:
    """
    Splits contigs greater than a threshold into multiple contigs, and corrects associated GTF file.

    Breakpoints are either provided or detected by splitting evenly and
    choosing intergenic midpoints.  All metadata about splits is stored
    in each SeqRecord.description field, so no special FASTA header block
    is required.
    """

    def __init__(
        self,
        fasta_file: str,
        gtf_file: str,
        new_fasta_file: str,
        new_gtf_file: str,
        length_threshold: int = int(5e8),
        split_locations: dict[str, list[int]] = None,
    ):
        # 1) determine split locations if needed
        if split_locations is None:
            split_locations = self.find_split_points(fasta_file, gtf_file, length_threshold)
        self.split_locations = split_locations

        # 2) save original lengths
        self.orig_lengths = self.read_fasta(fasta_file)

        # 3) build new SeqRecord list with descriptions
        records = self.read_fasta_and_split(fasta_file)

        # 4) write the new FASTA
        self.write_new_fasta(records, new_fasta_file)

        # 5) rewrite the GTF
        self.gtf_chrom_name_change(gtf_file, split_locations, new_gtf_file)

    @staticmethod
    def read_fasta(fasta_file: str) -> dict[str, int]:
        """Return dict of {chrom: length} for each record in FASTA."""
        seqs = {}
        for rec in SeqIO.parse(fasta_file, "fasta"):
            seqs[rec.id] = len(rec.seq)
        return seqs

    @staticmethod
    def read_gtf_as_bed(gtf_file: str, gene_only: bool = True) -> pd.DataFrame:
        """Load a GTF, extract genes as BED-style dataframe for intergenic detection."""
        df = pd.read_csv(
            gtf_file,
            sep="\t",
            comment="#",
            header=None,
            names=[
                "seqname",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute",
            ],
        )
        df["start"] = df["start"] - 1  # 0-based for BED
        df["name"] = df["attribute"].str.extract(r'gene_id "([^"]+)"')
        missing = df["name"].isna()
        df.loc[missing, "name"] = df.loc[missing, "attribute"].str.extract(r'gene_name "([^"]+)"')
        df["extra"] = df["attribute"]
        if gene_only:
            df = df[df["feature"] == "gene"]
        return df[["seqname", "start", "end", "name", "score", "strand", "extra"]]

    @staticmethod
    def find_nearest_intergenic_regions(
        genes: pd.DataFrame, split_point: int, num_regions: int = 50
    ) -> list[tuple[int, int]]:
        """Among the nearest `num_regions` genes, find all intergenic windows."""
        genes = genes.copy()
        genes["distance"] = (genes["start"] - split_point).abs()
        nearest = genes.nsmallest(num_regions, "distance").sort_values("start")
        sites = np.concatenate([nearest["start"].values, nearest["end"].values])
        intergenic = []
        prev_end = nearest["end"].values[0]
        for _, row in nearest.iterrows():
            if row["start"] > prev_end and row["start"] == sites[sites > prev_end].min():
                intergenic.append((prev_end, row["start"]))
            prev_end = row["end"]
        return intergenic

    @classmethod
    def find_split_points(
        cls, fasta_file: str, gtf_file: str, length_threshold: int
    ) -> dict[str, list[int]]:
        """Auto-detect breakpoints so that no piece exceeds `length_threshold`."""
        seq_lens = cls.read_fasta(fasta_file)
        genes = cls.read_gtf_as_bed(gtf_file)
        splits: dict[str, list[int]] = {}

        for chrom, length in seq_lens.items():
            if length <= length_threshold:
                continue
            n = int(np.ceil(length / length_threshold))
            theor = [i * length // n for i in range(1, n)]
            chrom_genes = genes[genes["seqname"] == chrom]
            opts: list[int] = []
            for t in theor:
                ig = cls.find_nearest_intergenic_regions(chrom_genes, t)
                best = max(ig, key=lambda x: x[1] - x[0])
                opts.append((best[0] + best[1]) // 2)
            splits[chrom] = opts
        return splits

    def read_fasta_and_split(self, fasta_file: str) -> list[SeqRecord]:
        """
        Builds SeqRecord objects with .description carrying:
          - orig=<chrom>:1-<full_length>
          - chunk=<chrom>__part<i>:<start>-<end> (1-based coords)
        """
        records: list[SeqRecord] = []
        for rec in SeqIO.parse(fasta_file, "fasta"):
            chrom = rec.id
            full = self.orig_lengths[chrom]
            if chrom in self.split_locations:
                edges = [0] + self.split_locations[chrom] + [full]
                for i in range(len(edges) - 1):
                    s0, s1 = edges[i], edges[i + 1]
                    pid = f"{chrom}__part{ i + 1 }"
                    # convert to 1-based inclusive for human readability
                    desc = (
                        f"orig={chrom}:1-{full} "
                        f"chunk={pid}:{s0+1}-{s1}"
                    )
                    records.append(SeqRecord(rec.seq[s0:s1], id=pid, description=desc))
            else:
                desc = f"orig={chrom}:1-{full}"
                records.append(SeqRecord(rec.seq, id=chrom, description=desc))
        return records

    def write_new_fasta(self, records: list[SeqRecord], out_fa: str):
        """Dump exactly our SeqRecord list to FASTA (descriptions become > headers)."""
        with open(out_fa, "w") as h:
            SeqIO.write(records, h, "fasta")

    @staticmethod
    def gtf_chrom_name_change(
        gtf_file: str, split_locations: dict[str, list[int]], output_file: str
    ):
        """Rewrite GTF rows to new __partN names and adjust coords by chunk offsets."""
        starts = {c: zcumsum(ls) for c, ls in split_locations.items()}
        with open(gtf_file) as gh, open(output_file, "w") as oh:
            rdr = csv.reader(gh, delimiter="\t")
            wtr = csv.writer(oh, delimiter="\t")
            for row in rdr:
                if not row or len(row) < 5:
                    continue
                if row[0].startswith("#"):
                    wtr.writerow(row)
                    continue
                chrom, s, e = row[0], int(row[3]), int(row[4])
                if chrom in starts:
                    offs = starts[chrom]
                    # find which chunk contains this start
                    for idx in range(len(offs) - 1):
                        if offs[idx] <= s < offs[idx + 1]:
                            row[0] = f"{chrom}__part{idx+1}"
                            row[3] = str(s - offs[idx])
                            row[4] = str(e - offs[idx])
                            break
                wtr.writerow(row)

    def read_split_locations_from_fasta_header(self, fasta_file: str) -> dict[str, list[int]]:
        """
        Scans each SeqRecord.description in FASTA, pulls out all
        chunk-end coordinates, and reconstructs split_locations.
        """
        locs: dict[str, list[int]] = {}
        for rec in SeqIO.parse(fasta_file, "fasta"):
            parts = {f.split("=", 1)[0]: f.split("=", 1)[1] for f in rec.description.split()}
            if "chunk" not in parts:
                continue
            orig_chr, _ = parts["orig"].split(":", 1)
            # chunk looks like "chr__partN:A-B"
            _, coord = parts["chunk"].split(":", 1)
            _, end = coord.split("-", 1)
            locs.setdefault(orig_chr, []).append(int(end))
        # drop the last end==full-length
        for c, ends in locs.items():
            ends.sort()
            locs[c] = [e for e in ends[:-1]]
        return locs

    def uncorrect_bed_positions(
        self,
        bed_file: str,
        output_file: str,
        split_locations: dict[str, list[int]] = None,
        split_fasta_file: str = None,
    ):
        """Map BED from split contigs back to original coords."""
        if split_locations is None:
            if split_fasta_file is None:
                raise ValueError("Must provide split_fasta_file if split_locations is None")
            split_locations = self.read_split_locations_from_fasta_header(split_fasta_file)

        # build chunk-starts for each chrom
        starts = {c: [0] + sorted(ls) for c, ls in split_locations.items()}

        with open(bed_file) as bh, open(output_file, "w") as oh:
            rdr = csv.reader(bh, delimiter="\t")
            wtr = csv.writer(oh, delimiter="\t")
            for row in rdr:
                if row[0].startswith("#"):
                    wtr.writerow(row)
                    continue
                chrom, s, e = row[0], int(row[1]), int(row[2])
                if "__part" in chrom:
                    orig, part = chrom.split("__part")
                    idx = int(part) - 1
                    off = starts[orig][idx]
                    s += off
                    e += off
                    row[0], row[1], row[2] = orig, str(s), str(e)
                wtr.writerow(row)

    def uncorrect_bigwig_positions(
        self,
        bigwig_file: str,
        output_file: str,
        split_locations: dict[str, list[int]] = None,
        split_fasta_file: str = None,
    ):
        """Map BigWig from split contigs back to original coordinates."""
        if split_locations is None:
            if split_fasta_file is None:
                raise ValueError("Must provide split_fasta_file if split_locations is None")
            split_locations = self.read_split_locations_from_fasta_header(split_fasta_file)

        starts = {c: [0] + sorted(ls) for c, ls in split_locations.items()}

        bw_in = pyBigWig.open(bigwig_file)
        # rebuild per-chrom size map
        chrom_sizes: dict[str, int] = {}
        for ch, size in bw_in.chroms().items():
            if "__part" in ch:
                orig, part = ch.split("__part")
                idx = int(part) - 1
                s0 = starts[orig][idx]
                s1 = (starts[orig] + [float("inf")])[idx + 1]
                chrom_sizes[orig] = chrom_sizes.get(orig, 0) + int(s1 - s0)
            else:
                chrom_sizes[ch] = size

        bw_out = pyBigWig.open(output_file, "w")
        bw_out.addHeader(list(chrom_sizes.items()))

        for ch in bw_in.chroms():
            intervals = bw_in.intervals(ch)
            if not intervals:
                continue
            if "__part" in ch:
                orig, part = ch.split("__part")
                idx = int(part) - 1
                off = starts[orig][idx]
                adj = [(s + off, e + off, v) for s, e, v in intervals]
                bw_out.addEntries(
                    [orig] * len(adj),
                    [i[0] for i in adj],
                    ends=[i[1] for i in adj],
                    values=[i[2] for i in adj],
                )
            else:
                bw_out.addEntries(
                    [ch] * len(intervals),
                    [i[0] for i in intervals],
                    ends=[i[1] for i in intervals],
                    values=[i[2] for i in intervals],
                )

        bw_in.close()
        bw_out.close()
