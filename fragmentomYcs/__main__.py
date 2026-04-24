"""Command-line interface for the fragmentomYcs Python pipeline.

One-time setup (run from the project root)::

    pip install -e .

After that, run from any directory::

    fragmentomYcs --mut <mut> --bam <bam> --fasta <fasta> [options]

Or equivalently::

    python -m fragmentomYcs --mut <mut> --bam <bam> --fasta <fasta> [options]
"""

import argparse

from fragmentomYcs.run_fragmentomYcs import run_fragmentomYcs


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m python",
        description="Fragment-level allelic analysis of cfDNA sequencing data.",
    )

    # --- Required ---
    req = p.add_argument_group("required")
    req.add_argument("--mut", required=True,
                     help="Mutation file (.vcf/.tsv) or 'chr:pos:ref:alt' string.")
    req.add_argument("--bam", required=True,
                     help="Indexed BAM file.")
    req.add_argument("--fasta", required=True,
                     help="Reference FASTA file (.fa / .fasta).")

    # --- Optional ---
    p.add_argument("--sample-id", default=None,
                   help="Sample identifier added to every output row.")
    p.add_argument("--output", default=None, dest="output_path",
                   help="Output TSV file path.  Prints summary to stdout if omitted.")
    p.add_argument("--n-cores", type=int, default=1,
                   help="Number of parallel worker processes (default: 1).")
    p.add_argument("--neg-offset", type=int, default=-600,
                   dest="neg_offset_mate_search",
                   help="Upstream BAM search window offset (default: -600).")
    p.add_argument("--pos-offset", type=int, default=600,
                   dest="pos_offset_mate_search",
                   help="Downstream BAM search window offset (default: 600).")
    p.add_argument("--zero-based", action="store_false", dest="one_based",
                   help="Treat positions in --mut as 0-based (default: 1-based).")
    p.add_argument("--report-5p-3p-bases", type=int, default=5,
                   dest="report_5p_3p_bases_fragment",
                   help="Number of fragment-end bases to report (default: 5).")
    p.add_argument("--report-bam-info", action="store_true",
                   help="Include raw BAM fields in the output.")
    p.add_argument("--report-softclip", action="store_true",
                   help="Include soft-clip counts at fragment extremities.")
    p.add_argument("--remove-softclip", action="store_true",
                   help="Trim soft-clipped bases before analysis.")
    p.add_argument("--retain-fail-qc", action="store_true",
                   help="Keep QC-failed fragments in the output.")
    p.add_argument("--apply-norm", action="store_true",
                   dest="apply_bcftools_norm",
                   help="Left-align variants (backend set in fragmentomYcs.cfg).")
    p.add_argument("--tmp-folder", default=None,
                   help="Temporary directory (defaults to system temp).")
    p.add_argument("--verbose", action="store_true",
                   help="Print progress messages.")

    return p


def main(argv=None) -> None:
    args = _build_parser().parse_args(argv)

    result = run_fragmentomYcs(
        mut=args.mut,
        bam=args.bam,
        fasta=args.fasta,
        sample_id=args.sample_id,
        neg_offset_mate_search=args.neg_offset_mate_search,
        pos_offset_mate_search=args.pos_offset_mate_search,
        one_based=args.one_based,
        report_bam_info=args.report_bam_info,
        report_softclip=args.report_softclip,
        report_5p_3p_bases_fragment=args.report_5p_3p_bases_fragment,
        remove_softclip=args.remove_softclip,
        retain_fail_qc=args.retain_fail_qc,
        apply_bcftools_norm=args.apply_bcftools_norm,
        tmp_folder=args.tmp_folder,
        output_path=args.output_path,
        verbose=args.verbose,
        n_cores=args.n_cores,
    )

    # result is None when output_path was provided (written to file)
    if result is not None:
        print(result.to_csv(sep="\t", index=False), end="")


if __name__ == "__main__":
    main()
