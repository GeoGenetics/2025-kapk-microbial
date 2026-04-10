from itertools import chain
from functools import reduce
from pathlib import Path


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def get_kegg_damaged_files(wildcards, suffix):
    files = expand(
        config["rdir"] + "/kegg-reads-anvio-damaged/{smp}" + suffix,
        smp=sample_label_read,
    )
    return files


def read_tables(file, suffix):
    if is_non_zero_file(file):
        smp = remove_suffix(pathlib.Path(file).stem, suffix)
        df = pd.read_csv(file, sep="\t")
        df["label"] = sample_label_dict_read[smp]
        if suffix == "_modules":
            df.columns = [
                "gene_coverage" if "gene_coverage" in x else x for x in df.columns
            ]
            df.columns = [
                "gene_detection" if "gene_detection" in x else x for x in df.columns
            ]
            df.columns = [
                "avg_coverage" if "avg_coverage" in x else x for x in df.columns
            ]
            df.columns = [
                "avg_detection" if "avg_detection" in x else x for x in df.columns
            ]
        else:
            df.columns = ["coverage" if "coverage" in x else x for x in df.columns]
            df.columns = ["detection" if "detection" in x else x for x in df.columns]
        return df


# def read_tables(file):
#     if is_non_zero_file(file):
#         df = pd.read_csv(file, sep=",")
#         return df


rule mmseqs_profile_kegg_reads_damaged_anvio_summary:
    input:
        hits_tables=lambda wc: get_kegg_damaged_files(wc, suffix="_hits.txt"),
        hits_modules_tables=lambda wc: get_kegg_damaged_files(
            wc, suffix="_module_paths.txt"
        ),
        modules_tables=lambda wc: get_kegg_damaged_files(wc, suffix="_modules.txt"),
    output:
        hits_summary=(
            config["rdir"] + "/kegg-reads-anvio-damaged/kegg-hits-summary.tsv.gz"
        ),
        hits_modules_summary=(
            config["rdir"]
            + "/kegg-reads-anvio-damaged/kegg-hits-modules-summary.tsv.gz"
        ),
        modules_summary=(
            config["rdir"] + "/kegg-reads-anvio-damaged/kegg-modules-summary.tsv.gz"
        ),
    threads: 1
    log:
        config["rdir"] + "/logs/kegg-reads-anvio-damaged/kegg.summary.log",
    benchmark:
        config["rdir"] + "/benchmarks/kegg-reads-anvio-damaged/kegg.summary.bmk"
    message:
        """--- Summarize KEGG results."""
    run:
        df = pd.concat(
            map(
                lambda file: read_tables(file, suffix="_hits"),
                fast_flatten(list({input.hits_tables})),
            )
        )
        df.to_csv(
            output.hits_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )

        df = pd.concat(
            map(
                lambda file: read_tables(file, suffix="_module_paths"),
                fast_flatten(list({input.hits_modules_tables})),
            )
        )
        df.to_csv(
            output.hits_modules_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )

        df = pd.concat(
            map(
                lambda file: read_tables(file, suffix="_modules"),
                fast_flatten(list({input.modules_tables})),
            )
        )
        df.to_csv(
            output.modules_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )
