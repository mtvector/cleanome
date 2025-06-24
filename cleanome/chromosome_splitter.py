import os
import csv
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pyBigWig

def zcumsum(iterable):
    cumulative = [0]
    total = 0
    for x in iterable:
        total += x
        cumulative.append(total)
    cumulative.append(float("inf"))
    return cumulative


class ChromosomeSplitter:
    """
    Splits contigs larger than length_threshold into parts,
    writes a sidecar “.splitchr” file recording each contig’s breakpoints,
    and corrects GTF / BED / BigWig downstream using that sidecar.
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
        # 1) detect or accept split locations
        if split_locations is None:
            split_locations = self.find_split_points(fasta_file, gtf_file, length_threshold)
        self.split_locations = split_locations

        # 2) read original lengths
        self.orig_lengths = self._read_lengths(fasta_file)

        # 3) build split SeqRecords
        records = self._make_split_records(fasta_file)

        # 4) write new FASTA and sidecar
        self._write_fasta_and_sidecar(records, new_fasta_file)

        # 5) correct GTF
        self._gtf_chrom_name_change(gtf_file, split_locations, new_gtf_file)

    @staticmethod
    def _read_lengths(fasta_file: str) -> dict[str, int]:
        """Return {chrom: length} for each record in FASTA."""
        lengths = {}
        for rec in SeqIO.parse(fasta_file, "fasta"):
            lengths[rec.id] = len(rec.seq)
        return lengths

    @staticmethod
    def _read_gtf_as_bed(gtf_file: str) -> pd.DataFrame:
        """Load GTF and return gene BED for intergenic calculations."""
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
        df["start"] -= 1
        df["name"] = df["attribute"].str.extract(r'gene_id "([^"]+)"')
        missing = df["name"].isna()
        df.loc[missing, "name"] = df.loc[missing, "attribute"].str.extract(r'gene_name "([^"]+)"')
        return df[df["feature"] == "gene"][["seqname", "start", "end", "name", "score", "strand", "attribute"]]

    @staticmethod
    def _find_intergenic(genes: pd.DataFrame, point: int, n: int = 50) -> list[tuple[int,int]]:
        g = genes.copy()
        g["dist"] = (g["start"] - point).abs()
        near = g.nsmallest(n, "dist").sort_values("start")
        sites = np.concatenate([near["start"].values, near["end"].values])
        inter = []
        prev = near["end"].values[0]
        for _, row in near.iterrows():
            if row["start"] > prev and row["start"] == sites[sites > prev].min():
                inter.append((prev, row["start"]))
            prev = row["end"]
        return inter

    @classmethod
    def find_split_points(
        cls, fasta_file: str, gtf_file: str, length_threshold: int
    ) -> dict[str, list[int]]:
        """Auto-detect breakpoints so each piece ≤ length_threshold."""
        lens = cls._read_lengths(fasta_file)
        genes = cls._read_gtf_as_bed(gtf_file)
        splits: dict[str, list[int]] = {}
        for chrom, L in lens.items():
            if L <= length_threshold:
                continue
            n = int(np.ceil(L / length_threshold))
            theory = [i * L // n for i in range(1, n)]
            cg = genes[genes["seqname"] == chrom]
            points = []
            for t in theory:
                ig = cls._find_intergenic(cg, t)
                best = max(ig, key=lambda x: x[1] - x[0])
                points.append((best[0] + best[1]) // 2)
            splits[chrom] = points
        return splits

    def _make_split_records(self, fasta_file: str) -> list[SeqRecord]:
        """
        Return SeqRecord list where each record.description carries no metadata,
        splitting each chromosome at the breakpoints in self.split_locations.
        """
        recs: list[SeqRecord] = []
        for rec in SeqIO.parse(fasta_file, "fasta"):
            chrom = rec.id
            full = self.orig_lengths[chrom]
            if chrom in self.split_locations:
                edges = [0] + self.split_locations[chrom] + [full]
                for i in range(len(edges) - 1):
                    s, e = edges[i], edges[i + 1]
                    pid = f"{chrom}__part{i+1}"
                    # no description needed
                    recs.append(SeqRecord(rec.seq[s:e], id=pid, description=""))
            else:
                recs.append(SeqRecord(rec.seq, id=chrom, description=""))
        return recs

    def _write_fasta_and_sidecar(self, records: list[SeqRecord], out_fasta: str):
        """
        Dump FASTA, then write a .splitchr sidecar alongside it
        recording self.split_locations.
        """
        # 1) fasta
        with open(out_fasta, "w") as h:
            SeqIO.write(records, h, "fasta")

        # 2) sidecar
        base, _ = os.path.splitext(out_fasta)
        sidecar = base + ".splitchr"
        with open(sidecar, "w") as h:
            for chrom, pts in self.split_locations.items():
                h.write(f"{chrom}\t{','.join(map(str, pts))}\n")

    def _gtf_chrom_name_change(
        self,
        gtf_file: str,
        split_locations: dict[str, list[int]],
        output_file: str,
    ):
        """Rewrite GTF: remap contig names and adjust coords by chunk offsets."""
        starts = {c: zcumsum(locs) for c, locs in split_locations.items()}
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
                    for idx in range(len(offs) - 1):
                        if offs[idx] <= s < offs[idx + 1]:
                            row[0] = f"{chrom}__part{idx+1}"
                            row[3] = str(s - offs[idx])
                            row[4] = str(e - offs[idx])
                            break
                wtr.writerow(row)

    def _read_sidecar(self, fasta_or_splitfile: str) -> dict[str, list[int]]:
        """
        Read .splitchr file next to fasta_or_splitfile.
        """
        base, _ = os.path.splitext(fasta_or_splitfile)
        sidecar = base + ".splitchr"
        locs: dict[str, list[int]] = {}
        with open(sidecar) as h:
            for line in h:
                chrom, pts = line.strip().split("\t")
                if pts:
                    locs[chrom] = [int(x) for x in pts.split(",")]
        return locs

    def uncorrect_bed_positions(
        self,
        bed_file: str,
        output_file: str,
        split_locations: dict[str, list[int]] = None,
        split_fasta_file: str = None,
    ):
        """Map BED from split contigs back to original coords using sidecar."""
        if split_locations is None:
            if split_fasta_file is None:
                raise ValueError("Must provide split_fasta_file or split_locations")
            split_locations = self._read_sidecar(split_fasta_file)

        starts = {c: [0] + sorted(locs) for c, locs in split_locations.items()}

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
        """Map BigWig from split contigs back to original coordinates using sidecar."""
        if split_locations is None:
            if split_fasta_file is None:
                raise ValueError("Must provide split_fasta_file or split_locations")
            split_locations = self._read_sidecar(split_fasta_file)

        starts = {c: [0] + sorted(locs) for c, locs in split_locations.items()}

        bw_in = pyBigWig.open(bigwig_file)
        # rebuild chrom sizes
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
            ivals = bw_in.intervals(ch) or []
            if "__part" in ch:
                orig, part = ch.split("__part")
                idx = int(part) - 1
                off = starts[orig][idx]
                adj = [(s + off, e + off, v) for s, e, v in ivals]
                bw_out.addEntries(
                    [orig] * len(adj),
                    [i[0] for i in adj],
                    ends=[i[1] for i in adj],
                    values=[i[2] for i in adj],
                )
            else:
                bw_out.addEntries(
                    [ch] * len(ivals),
                    [i[0] for i in ivals],
                    ends=[i[1] for i in ivals],
                    values=[i[2] for i in ivals],
                )

        bw_in.close()
        bw_out.close()
