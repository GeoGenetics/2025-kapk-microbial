from itertools import chain
from functools import reduce
from pathlib import Path
import re


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def get_tp_mapping_inputs(wildcards, indata, k, ani):
    files = expand(
        config["rdir"]
        + "/taxonomic-profiling-mapping/"
        + indata
        + "/k"
        + k
        + "/{smp}.dedup.metrics",
        smp=sample_label_read,
    )
    return files


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


def remove_suffix(input_string, suffix):
    if suffix and input_string.endswith(suffix):
        return input_string[: -len(suffix)]
    return input_string


def read_picard_tp_mapping_tables(file):
    if is_non_zero_file(file):
        try:
            smp = remove_suffix(pathlib.Path(file).stem, ".dedup")
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


def get_tp_stats_files(wildcards, indata, k, ani):
    files = expand(
        config["rdir"]
        + "/taxonomic-profiling-filtering/"
        + indata
        + "/k"
        + k
        + "/"
        + ani
        + "/{smp}.dedup_stats.tsv.gz",
        smp=sample_label_read,
    )
    return files


def get_tp_stats_filtered_files(wildcards, indata, k, ani):
    files = expand(
        config["rdir"]
        + "/taxonomic-profiling-filtering/"
        + indata
        + "/k"
        + k
        + "/"
        + ani
        + "/{smp}.dedup_stats-filtered.tsv.gz",
        smp=sample_label_read,
    )
    return files


def read_tp_stats_tables(file, suffix):
    if is_non_zero_file(file):
        smp = remove_suffix(pathlib.Path(file).stem, suffix)
        try:
            df = pd.read_csv(file, sep="\t")
            df["label"] = sample_label_dict_read[smp]
            return df
        except pd.errors.EmptyDataError:
            pass


rule taxonomic_profiling_mapping_summary:
    wildcard_constraints:
        indata="prefilter|standard|raw",
    input:
        tp_map_stats=lambda wc: get_tp_mapping_inputs(
            wc, indata=wc.indata, k=wc.tp_bowtie2_k, ani=wc.read_ani
        ),
        tp_stats=lambda wc: get_tp_stats_files(
            wc, indata=wc.indata, k=wc.tp_bowtie2_k, ani=wc.read_ani
        ),
        tp_stats_filtered=lambda wc: get_tp_stats_filtered_files(
            wc, indata=wc.indata, k=wc.tp_bowtie2_k, ani=wc.read_ani
        ),
    output:
        tp_summary_map_stats=(
            config["rdir"]
            + "/taxonomic-profiling-summary/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/tp-summary-mapping.tsv.gz"
        ),
        tp_stats_summary=(
            config["rdir"]
            + "/taxonomic-profiling-summary/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/tp-mapping.summary.tsv.gz"
        ),
        tp_stats_filtered_summary=(
            config["rdir"]
            + "/taxonomic-profiling-summary/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/tp-mapping-filtered.summary.tsv.gz"
        ),
    threads: 1
    log:
        config["rdir"]
        + "/logs/taxonomic-profiling-summary/"
        + "{indata}"
        + "/k{tp_bowtie2_k}"
        + "/{read_ani}"
        + "/tp-summary.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/taxonomic-profiling-summary/"
            + "{indata}"
            + "/k{tp_bowtie2_k}"
            + "/{read_ani}"
            + "/tp-summary.bmk"
        )
    message:
        """--- Summarize taxonomic profiling."""
    run:
        df = pd.concat(
            map(
                lambda file: read_picard_tp_mapping_tables(file),
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
                lambda file: read_tp_stats_tables(file, suffix=".dedup_stats.tsv"),
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
                lambda file: read_tp_stats_tables(
                    file, suffix=".dedup_stats-filtered.tsv"
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
