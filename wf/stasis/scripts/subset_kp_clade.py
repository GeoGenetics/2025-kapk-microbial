#!/usr/bin/env python3
"""
Subset the concatenated codon alignment and ML tree to the Kap K clade
for ClonalFrameML recombination analysis.

Kap K clade = taxa with ≥98% ANI to ancient_consensus, plus the focal genomes.
These are the SOIL, NAY, GCA054* and KapK genomes.
"""

import sys
import re
from pathlib import Path

from Bio import SeqIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Kap K clade taxa — all genomes with >98% ANI and the focal KapK genomes.
# Based on FastANI results and phylogenomics tree topology.
KAP_K_CLADE = {
    "ancient_consensus",
    "modern_consensus",
    "s17_ancient_relative",
    "SOIL100000032",
    "SOIL100000033",
    "SOIL100000034",
    "SOIL100000035",
    "SOIL100000037",
    "SOIL100000038",
    "SOIL100000040",
    "SOIL100000041",
    "SOIL100000042",
    "SOIL100000043",
    "NAY3300025461b7",
    "NAY3300025446b16",
    "NAY3300025496b14",
    "GCA054338395",
    "GCA054338475",
    "GCA054313015",
    "GCA054313355",
    "GCA054313435",
    "GCA054328175",
    "GCA054328195",
}


def strip_marker_suffix(name: str) -> str:
    for sep in ["_COG", "_DNGNGWU", "_FAM"]:
        idx = name.find(sep)
        if idx >= 0:
            return name[:idx]
    return name


def subset_fasta(fasta_in: Path, fasta_out: Path, keep_taxa: set) -> list:
    """Write subset FASTA; returns list of taxa actually written."""
    kept = []
    records = []
    for rec in SeqIO.parse(fasta_in, "fasta"):
        taxon = strip_marker_suffix(rec.id)
        if taxon in keep_taxa:
            records.append(SeqRecord(rec.seq, id=taxon, description=""))
            kept.append(taxon)
    SeqIO.write(records, fasta_out, "fasta")
    return kept


def prune_tree(tree_in: Path, tree_out: Path, keep_taxa: set):
    """Prune tree to keep_taxa; writes Newick."""
    tree = Phylo.read(tree_in, "newick")

    # Collect all terminal names in the tree
    all_leaves = {c.name for c in tree.get_terminals()}
    to_prune = all_leaves - keep_taxa

    # Prune by collapsing terminals not in keep set
    for leaf_name in to_prune:
        target = tree.find_any(leaf_name)
        if target:
            try:
                tree.prune(target)
            except Exception:
                pass

    Phylo.write(tree, tree_out, "newick")


def main():
    if len(sys.argv) < 5:
        print(f"Usage: {sys.argv[0]} <concat_codon.fna> <ml_tree.treefile> "
              f"<out_fasta> <out_tree>")
        sys.exit(1)

    fasta_in  = Path(sys.argv[1])
    tree_in   = Path(sys.argv[2])
    fasta_out = Path(sys.argv[3])
    tree_out  = Path(sys.argv[4])

    fasta_out.parent.mkdir(parents=True, exist_ok=True)

    print(f"Subsetting {fasta_in.name} to {len(KAP_K_CLADE)} Kap K clade taxa",
          file=sys.stderr)

    kept = subset_fasta(fasta_in, fasta_out, KAP_K_CLADE)
    print(f"  Written {len(kept)} sequences to {fasta_out}", file=sys.stderr)

    missing = KAP_K_CLADE - set(kept)
    if missing:
        print(f"  Warning: {len(missing)} taxa not found in alignment: "
              f"{sorted(missing)}", file=sys.stderr)

    prune_tree(tree_in, tree_out, set(kept))
    print(f"  Written pruned tree to {tree_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
