"""Compare the integration run's EN-TEx genotypes with ENCODE's GRCh38 calls.

ENCODE called these donors from 10x linked reads with freeBayes, so the call
sets agree closely but not exactly. Both sides are split into biallelic
records, left-aligned against the reference and restricted to the calling
regions. Run from the repository root after the integration run:

    pixi run --manifest-path workflow/envs/pixi.toml python .tests/integration/compare_encode.py
"""

import argparse
import subprocess
import tempfile
from pathlib import Path

RESOURCES = Path("/data/zusers/ramirezc/testing/MOHD-DACC-WGS-Pipeline/resources")
ENCODE_DATA = Path("/data/projects/encode/data")
# Sample ID of the integration run -> ENCODE's released GRCh38 variant calls.
ENCODE_CALLS = {
    "ENCSR410ALS_rep1": ENCODE_DATA / "ENCSR410ALS/ENCFF501WAL.vcf.gz",
    "ENCSR456SNK_rep1": ENCODE_DATA / "ENCSR456SNK/ENCFF568LEP.vcf.gz",
}


def build_parser():
    """Build the compare_encode.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the compare_encode.py flags.
    """
    p = argparse.ArgumentParser(description="Compare genotypes with ENCODE's.")
    p.add_argument("--run", default="results/encode-integration", help="Run directory.")
    p.add_argument(
        "--fasta",
        default=str(RESOURCES / "GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"),
        help="Reference FASTA of both call sets.",
    )
    p.add_argument(
        "--intervals",
        default=str(RESOURCES / "wgs_calling_regions.hg38.interval_list"),
        help="Calling regions, as a Picard interval list.",
    )
    return p


def interval_list_to_bed(interval_list, bed):
    """Write a Picard interval list as a BED file.

    Parameters
    ----------
    interval_list : str
        Picard interval list: 1-based, closed intervals after `@` header lines.
    bed : str
        BED file to write.
    """
    with open(interval_list) as src, open(bed, "w") as dst:
        for line in src:
            if not line.startswith("@"):
                contig, start, end = line.split("\t")[:3]
                dst.write(f"{contig}\t{int(start) - 1}\t{end}\n")


def genotypes(vcf, fasta, bed, contigs, include):
    """Return the alternate genotypes of a VCF's first sample.

    Parameters
    ----------
    vcf : str
        Gzipped VCF or gVCF.
    fasta : str
        Reference FASTA, for left-alignment.
    bed : str
        Regions to keep.
    contigs : str
        `##contig` header lines of the reference. They replace the VCF's own.
    include : str
        Filter expression, in bcftools syntax, for the records to keep.

    Returns
    -------
    dict of tuple to str
        `(chrom, pos, ref, alt)` to the unphased genotype, `het` or `hom`.
    """
    pipeline = (
        # ENCODE's VCFs are plain gzip, lack the ##contig lines and use INFO tags
        # their header does not define. Only the genotypes are compared.
        f"zcat {vcf}"
        " | grep -v '^##contig='"
        " | awk -F '\\t' -v OFS='\\t'"
        f' \'/^#CHROM/ {{ while ((getline c < "{contigs}") > 0) print c }}'
        ' /^#/ { print; next } { $8 = "."; print }\''
        f" | bcftools view -T {bed} -i '{include}'"
        f" | bcftools norm -c x -m -any -f {fasta}"
        ' | bcftools view -e \'ALT="<NON_REF>" || ALT="*"\''
        " | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n'"
    )
    out = subprocess.run(
        ["bash", "-o", "pipefail", "-c", pipeline],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout
    calls = {}
    for line in out.splitlines():
        chrom, pos, ref, alt, gt = line.split("\t")
        alleles = gt.replace("|", "/").split("/")
        if "1" in alleles:
            zygosity = "hom" if alleles.count("1") == len(alleles) else "het"
            calls[(chrom, int(pos), ref, alt)] = zygosity
    return calls


def summarize(ours, encode):
    """Return recall, agreement and genotype concordance for one variant class.

    Parameters
    ----------
    ours, encode : dict of tuple to str
        Genotypes from `genotypes`.

    Returns
    -------
    tuple of int and float
        ENCODE count, our count, ENCODE calls we make, our calls ENCODE makes,
        and the share of shared sites with the same zygosity.
    """
    shared = ours.keys() & encode.keys()
    same = sum(ours[k] == encode[k] for k in shared)
    return (
        len(encode),
        len(ours),
        len(shared) / len(encode),
        len(shared) / len(ours),
        same / len(shared),
    )


def main(argv=None):
    """Print one comparison row per sample and variant class.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with tempfile.TemporaryDirectory() as tmp:
        bed = f"{tmp}/calling_regions.bed"
        interval_list_to_bed(args.intervals, bed)
        contigs = f"{tmp}/contigs.txt"
        with open(f"{args.fasta}.fai") as fai, open(contigs, "w") as out:
            for line in fai:
                name, length = line.split("\t")[:2]
                out.write(f"##contig=<ID={name},length={length}>\n")
        print("sample\tclass\tencode\tours\trecall\tagreement\tgenotype_match")
        for sample, encode_vcf in ENCODE_CALLS.items():
            gvcf = f"{args.run}/{sample}/reblock/{sample}.rb.g.vcf.gz"
            ours = genotypes(gvcf, args.fasta, bed, contigs, 'GT="alt"')
            encode = genotypes(encode_vcf, args.fasta, bed, contigs, 'FILTER="PASS"')
            for name, is_class in [
                ("snv", lambda k: len(k[2]) == len(k[3]) == 1),
                ("indel", lambda k: len(k[2]) != len(k[3])),
            ]:
                row = summarize(
                    {k: v for k, v in ours.items() if is_class(k)},
                    {k: v for k, v in encode.items() if is_class(k)},
                )
                print(
                    f"{sample}\t{name}\t{row[0]}\t{row[1]}\t"
                    + "\t".join(f"{x:.4f}" for x in row[2:])
                )


if __name__ == "__main__":
    main()
