from itertools import chain
from functools import reduce
from pathlib import Path


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 64


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


# def get_pmd_files(wildcards, weight, dmg_mode):
#     files = expand(
#         config["rdir"]
#         + "/taxonomic-profiling-dmg/" + "{indata}" + "/k{tp_bowtie2_k}" + "/{smp}.tp-mdmg."
#         + ".pmd"
#         + ".weight-"
#         + weight
#         + "--"
#         + dmg_mode
#         + ".txt.gz",
#         smp=sample_label_read,
#     )
#     return files


def get_dmg_files_dr(wildcards, weight, indata, k, ani, mag_db, dr_fraction):
    files = expand(
        config["rdir"]
        + "/reads-dmg-distribution-mags-metaDMG/"
        + indata
        + "/k"
        + k
        + "/"
        + ani
        + "/local"
        + "/"
        + mag_db
        + "/"
        + dr_fraction
        + "/{smp}.tp-mdmg.weight-"
        + weight
        + ".csv.gz",
        smp=sample_label_read,
    )
    return files


def read_dmg_tables(file):
    if is_non_zero_file(file):
        df = pd.read_csv(file, sep=",")
        return df


def read_pmd_tables(file):
    if is_non_zero_file(file):
        smp = Path(file).stem.split(".")[0]
        df = pd.read_csv(file, sep="\t", header=None)
        df.columns = ["read", "pmd_score"]
        df["label"] = sample_label_dict_read[smp]
        return df


rule reads_dmg_distribution_metaDMG_summary:
    wildcard_constraints:
        indata="prefilter|standard|raw",
        weight="0|1",
        dr_fraction="damaged|non-damaged",
    input:
        dmg_tables=lambda wc: get_dmg_files_dr(
            wc,
            weight=wc.weight,
            indata=wc.indata,
            k=wc.tp_bowtie2_k,
            ani=wc.read_ani,
            mag_db=wc.mag_db,
            dr_fraction=wc.dr_fraction,
        ),
        # pmd_tables=lambda wc: get_pmd_files(
        #     wc, weight=wc.weight, dmg_mode=wc.dmg_mode
        # ),
    output:
        dmg_summary=(
            config["rdir"]
            + "/reads-dmg-distribution-mags-metaDMG/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/local"
            + "/{mag_db}"
            + "/{dr_fraction}"
            + "/tp-mdmg.weight-{weight}.csv.gz"
        ),
        # dmg_pmd_summary=(
        #     config["rdir"]
        #     + "/taxonomic-profiling-dmg/tp-mdmg.{dmg_mode}.pmd.weight-{weight}.csv.gz"
        # ),
    threads: 1
    log:
        config["rdir"]
        + "/logs/reads-dmg-distribution-mags-metaDMG/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/local"
        + "/{mag_db}"
        + "/{dr_fraction}"
        + "/tp-mdmg.weight-{weight}.log",
    benchmark:
        (
            config["rdir"]
            + "/reads-dmg-distribution-mags-metaDMG/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/local"
            + "/{mag_db}"
            + "/{dr_fraction}"
            + "/tp-mdmg.weight-{weight}.bmk"
        )
    message:
        """--- Summarize metaDMG results."""
    run:
        df = pd.concat(
            map(
                lambda file: read_dmg_tables(file),
                fast_flatten(list({input.dmg_tables})),
            )
        )
        df["label"] = df["sample"].astype(str).map(sample_label_dict_read)
        df.to_csv(
            output.dmg_summary,
            sep=",",
            compression="gzip",
            header=True,
            index=False,
        )

        # df = pd.concat(
        #     map(
        #         lambda file: read_pmd_tables(file),
        #         fast_flatten(list({input.pmd_tables})),
        #     )
        # )
        # df.to_csv(
        #     output.dmg_pmd_summary,
        #     sep=",",
        #     compression="gzip",
        #     header=True,
        #     index=False,
        # )
