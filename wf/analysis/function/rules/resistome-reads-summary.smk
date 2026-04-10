from itertools import chain
from functools import reduce
from pathlib import Path


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def get_resistome_files(wildcards, path, suffix):
    files = expand(
        config["rdir"] + "/" + path + "/{smp}" + suffix,
        smp=sample_label_read,
    )
    return files


def read_resistome_tables(file, suffix):
    if is_non_zero_file(file):
        smp = remove_suffix(pathlib.Path(file).stem, suffix)
        df = pd.read_csv(file, sep="\t")
        df["label"] = sample_label_dict_read[smp]
        return df


# def read_tables(file):
#     if is_non_zero_file(file):
#         df = pd.read_csv(file, sep=",")
#         return df


rule resistome_reads_summary:
    input:
        stats_tables=lambda wc: get_resistome_files(
            wc, path="resistome-reads", suffix=".megares_stats.tsv.gz"
        ),
        stats_filtered_tables=lambda wc: get_resistome_files(
            wc, path="resistome-reads", suffix=".megares_stats-filtered.tsv.gz"
        ),
    output:
        stats_summary=(
            config["rdir"] + "/resistome-reads/resistome_reads-stats-summary.tsv.gz"
        ),
        stats_filtered_summary=(
            config["rdir"]
            + "/resistome-reads/resistome_reads-stats-filtered-summary.tsv.gz"
        ),
    threads: 1
    log:
        config["rdir"] + "/logs/resistome-reads/resistome.summary.log",
    benchmark:
        config["rdir"] + "/benchmarks/resistome-reads/resistome.summary.bmk"
    message:
        """--- Summarize resistome results."""
    run:
        # Stats tables
        df = pd.concat(
            map(
                lambda file: read_resistome_tables(file, suffix=".megares_stats.tsv"),
                fast_flatten(list({input.stats_tables})),
            )
        )
        df.to_csv(
            output.stats_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )
        # Stats filtered table
        df = pd.concat(
            map(
                lambda file: read_resistome_tables(
                    file, suffix=".megares_stats-filtered.tsv"
                ),
                fast_flatten(
                    list({input.stats_filtered_tables}),
                ),
            )
        )
        df.to_csv(
            output.stats_filtered_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )
