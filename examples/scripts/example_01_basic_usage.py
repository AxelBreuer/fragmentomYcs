# -*- coding: utf-8 -*-

from pathlib import Path

from fragmentomYcs.plot_size_distribution import plot_size_distribution
from fragmentomYcs.run_fragmentomYcs import run_fragmentomYcs

script_name = "fragmentomYcs.example_01_basic_usage"

# input files paths (relative to this script's location)
_extdata = Path(__file__).parent.parent.parent / "inst" / "extdata"
mut_file   = _extdata / "mutation" / "cfdna-egfr-del_chr7_55241864_55243064_10k.mutations.tsv"
bam_file   = _extdata / "bam"      / "cfdna-egfr-del_chr7_55241864_55243064_10k.bam"
fasta_file = _extdata / "fasta"    / "hg19_chr7_55231864_55253064.fa"

print(script_name + f": mut_file = {mut_file}")
print(script_name + f": bam_file = {bam_file}")
print(script_name + f": fasta_file = {fasta_file}")

df_fragments = run_fragmentomYcs(
    mut=mut_file,
    bam=bam_file,
    fasta=fasta_file,
    sample_id="cfdna-egfr-del",
    apply_bcftools_norm=True,
    n_cores=1,
)

print(script_name + ": df_fragments =")
print(df_fragments)

plot_size_distribution(
    df_fragments,
    vals_z=["MUT", "WT"],
    show_histogram=True,
    show_density=True,
    x_limits=(100, 420),
    histo_args={"alpha": 0.25},
    density_args={"linewidth": 2},
    histogram_binwidth=10,
    colors_z=["#F6BD60", "#84A59D"],
)