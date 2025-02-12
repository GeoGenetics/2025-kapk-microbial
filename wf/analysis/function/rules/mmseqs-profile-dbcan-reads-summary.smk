from itertools import chain
from functools import reduce
from pathlib import Path


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def get_dbcan_files(wildcards, suffix):
    files = expand(
        config["rdir"] + "/dbcan-reads/{smp}" + suffix,
        smp=sample_label_read,
    )
    return files


def read_dbcan_tables(file, suffix):
    if is_non_zero_file(file):
        smp = remove_suffix(pathlib.Path(file).stem, suffix)
        df = pd.read_csv(file, sep="\t")
        df["label"] = sample_label_dict_read[smp]
        return df


# def read_tables(file):
#     if is_non_zero_file(file):
#         df = pd.read_csv(file, sep=",")
#         return df


rule mmseqs_profile_dbcan_reads_summary:
    input:
        hits_tables=lambda wc: get_dbcan_files(
            wc, suffix=".profile.dbcan_group-abundances-agg.tsv.gz"
        ),
    output:
        group_abundance_summary=(
            config["rdir"] + "/dbcan-reads/dbcan.group-abundances-agg.tsv.gz"
        ),
    threads: 1
    log:
        config["rdir"] + "/logs/dbcan-reads/dbcan.summary.log",
    benchmark:
        config["rdir"] + "/benchmarks/dbcan-reads/dbcan.summary.bmk"
    message:
        """--- Summarize dbcan results."""
    run:
        df = pd.concat(
            map(
                lambda file: read_dbcan_tables(
                    file, suffix=".profile.dbcan_group-abundances-agg.tsv"
                ),
                fast_flatten(list({input.hits_tables})),
            )
        )
        df.to_csv(
            output.group_abundance_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )
