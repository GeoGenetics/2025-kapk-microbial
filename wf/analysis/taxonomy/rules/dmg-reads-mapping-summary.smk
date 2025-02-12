from itertools import chain
from functools import reduce
from pathlib import Path
import re


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def get_tp_mapping_inputs_bm(
    wildcards, indata, k, ani, weight, dr_fraction, bm_read_ani
):
    files = expand(
        config["rdir"]
        + "/reads-dmg-mapping/"
        + indata
        + "/k"
        + k
        + "/"
        + ani
        + "/weight-"
        + weight
        + "/"
        + dr_fraction
        + "/{smp}.dmg-reads.dedup.metrics",
        smp=sample_label_read,
    )
    return files


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


def remove_suffix(input_string, suffix):
    if suffix and input_string.endswith(suffix):
        return input_string[: -len(suffix)]
    return input_string


def read_picard_tp_mapping_tables_dmg_dist(file):
    if is_non_zero_file(file):
        try:
            smp = remove_suffix(pathlib.Path(file).stem, ".dmg-reads.dedup")
            df = pd.read_csv(file, sep="\t", comment="#")
            df["TOTAL_READS"] = (
                df["UNPAIRED_READS_EXAMINED"]
                + df["READ_PAIRS_EXAMINED"]
                + df["SECONDARY_OR_SUPPLEMENTARY_RDS"]
                + df["UNMAPPED_READS"]
            )
            df["TOTAL_READS_DEDUP"] = (
                df["UNPAIRED_READ_DUPLICATES"] + df["READ_PAIR_OPTICAL_DUPLICATES"]
            )
            df["nREADS"] = df["TOTAL_READS"] - df["TOTAL_READS_DEDUP"]
            df["label"] = sample_label_dict_read[smp]
            return df
        except pd.errors.EmptyDataError:
            pass


def get_tp_stats_files_bm(wildcards, indata, k, ani, weight, dr_fraction, bm_read_ani):
    files = expand(
        config["rdir"]
        + "/reads-dmg-filtering/"
        + indata
        + "/k"
        + k
        + "/"
        + ani
        + "/weight-"
        + weight
        + "/"
        + dr_fraction
        + "/"
        + bm_read_ani
        + "/{smp}.dmg-reads.dedup_stats.tsv.gz",
        smp=sample_label_read,
    )
    return files


def get_tp_stats_filtered_files_bm(
    wildcards, indata, k, ani, weight, dr_fraction, bm_read_ani
):
    files = expand(
        config["rdir"]
        + "/reads-dmg-filtering/"
        + indata
        + "/k"
        + k
        + "/"
        + ani
        + "/weight-"
        + weight
        + "/"
        + dr_fraction
        + "/"
        + bm_read_ani
        + "/{smp}.dmg-reads.dedup_stats-filtered.tsv.gz",
        smp=sample_label_read,
    )
    return files


def read_tp_stats_tables_dmg_dist(file, suffix):
    if is_non_zero_file(file):
        smp = remove_suffix(pathlib.Path(file).stem, suffix)
        df = pd.read_csv(file, sep="\t")
        df["label"] = sample_label_dict_read[smp]
        return df


rule reads_dmg_summary:
    wildcard_constraints:
        indata="prefilter|standard|raw",
        weight="0|1",
        dr_fraction="damaged|non-damaged",
    input:
        tp_map_stats=lambda wc: get_tp_mapping_inputs_bm(
            wc,
            indata=wc.indata,
            k=wc.tp_bowtie2_k,
            ani=wc.read_ani,
            weight=wc.weight,
            dr_fraction=wc.dr_fraction,
            bm_read_ani=wc.bm_read_ani,
        ),
        tp_stats=lambda wc: get_tp_stats_files_bm(
            wc,
            indata=wc.indata,
            k=wc.tp_bowtie2_k,
            ani=wc.read_ani,
            weight=wc.weight,
            dr_fraction=wc.dr_fraction,
            bm_read_ani=wc.bm_read_ani,
        ),
        tp_stats_filtered=lambda wc: get_tp_stats_filtered_files_bm(
            wc,
            indata=wc.indata,
            k=wc.tp_bowtie2_k,
            ani=wc.read_ani,
            weight=wc.weight,
            dr_fraction=wc.dr_fraction,
            bm_read_ani=wc.bm_read_ani,
        ),
    output:
        tp_summary_map_stats=(
            config["rdir"]
            + "/reads-dmg-summary/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{dr_fraction}"
            + "/{bm_read_ani}"
            + "/tp-summary-mapping.tsv.gz"
        ),
        tp_stats_summary=(
            config["rdir"]
            + "/reads-dmg-summary/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{dr_fraction}"
            + "/{bm_read_ani}"
            + "/tp-mapping.summary.tsv.gz"
        ),
        tp_stats_filtered_summary=(
            config["rdir"]
            + "/reads-dmg-summary/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{dr_fraction}"
            + "/{bm_read_ani}"
            + "/tp-mapping-filtered.summary.tsv.gz"
        ),
    threads: 1
    log:
        config["rdir"]
        + "/logs/reads-dmg-summary/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/weight-{weight}"
        + "/{dr_fraction}"
        + "/{bm_read_ani}"
        + "/tp-summary.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/reads-dmg-summary/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/weight-{weight}"
            + "/{dr_fraction}"
            + "/{bm_read_ani}"
            + "/tp-summary.bmk"
        )
    message:
        """--- Summarize taxonomic profiling."""
    run:
        df = pd.concat(
            map(
                lambda file: read_picard_tp_mapping_tables_dmg_dist(file),
                fast_flatten(list({input.tp_map_stats})),
            )
        )
        df.to_csv(
            output.tp_summary_map_stats,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )

        df = pd.concat(
            map(
                lambda file: read_tp_stats_tables_dmg_dist(
                    file, suffix=".dmg-reads.dedup_stats.tsv"
                ),
                fast_flatten(list({input.tp_stats})),
            )
        )
        df.to_csv(
            output.tp_stats_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )
        df = pd.concat(
            map(
                lambda file: read_tp_stats_tables_dmg_dist(
                    file, suffix=".dmg-reads.dedup_stats-filtered.tsv"
                ),
                fast_flatten(list({input.tp_stats_filtered})),
            )
        )
        df.to_csv(
            output.tp_stats_filtered_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )
