from itertools import chain
import pathlib


def get_stats_read(wildcards, file_paths):
    files = expand(
        file_paths,
        smp=sample_label_read,
    )
    return files


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


def remove_suffix(input_string, suffix):
    if suffix and input_string.endswith(suffix):
        return input_string[: -len(suffix)]
    return input_string


rule stats_summary_reads:
    input:
        stats_derep=lambda wc: get_stats_read(
            wc, file_paths=config["rdir"] + "/stats/{smp}.stats-derep.txt"
        ),
        stats_initial=lambda wc: get_stats_read(
            wc, file_paths=config["rdir"] + "/stats/{smp}.stats-initial.txt"
        ),
        stats_extension=lambda wc: get_stats_read(
            wc, file_paths=config["rdir"] + "/stats/{smp}.stats-extension.txt"
        ),
    output:
        stats_derep_summary=config["rdir"] + "/stats/all.stats-derep-summary.tsv.gz",
        stats_initial_summary=(
            config["rdir"] + "/stats/all.stats-initial-summary.tsv.gz"
        ),
        stats_extension_summary=(
            config["rdir"] + "/stats/all.stats-extension-summary.tsv.gz"
        ),
    threads: 1
    log:
        config["rdir"] + "/logs/stats/stats-summary-reads.log",
    benchmark:
        config["rdir"] + "/benchmarks/stats/stats-summary-reads.bmk"
    message:
        """--- Summarize stats for reads"""
    run:
        df = pd.concat(
            map(
                lambda file: pd.read_csv(file, sep="\t"),
                fast_flatten(list({input.stats_derep})),
            )
        )
        df["label"] = df["file"].apply(
            lambda x: remove_suffix(pathlib.Path(x).stem, ".fa")
        )
        df["label"] = df["label"].astype(str).map(sample_label_dict_read)
        df.to_csv(
            output.stats_derep_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )

        df = pd.concat(
            map(
                lambda file: pd.read_csv(file, sep="\t"),
                fast_flatten(list({input.stats_initial})),
            )
        )
        df["label"] = df["file"].apply(
            lambda x: remove_suffix(pathlib.Path(x).stem, ".fq")
        )
        df["label"] = df["label"].astype(str).map(sample_label_dict_read)
        df.to_csv(
            output.stats_initial_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )

        df = pd.concat(
            map(
                lambda file: pd.read_csv(file, sep="\t"),
                fast_flatten(list({input.stats_extension})),
            )
        )
        df["label"] = df["file"].apply(
            lambda x: remove_suffix(pathlib.Path(x).stem, ".extended.fastq")
        )
        df["label"] = df["label"].astype(str).map(sample_label_dict_read)
        df.to_csv(
            output.stats_extension_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )
