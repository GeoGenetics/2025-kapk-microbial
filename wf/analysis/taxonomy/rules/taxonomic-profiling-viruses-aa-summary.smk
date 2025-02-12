from itertools import chain
from functools import reduce
from pathlib import Path


def is_non_zero_file(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


def get_viral_files(wildcards, suffix):
    files = expand(
        config["rdir"]
        + "/taxonomic-profiling-filtering-aa/"
        + "{indata}"
        + "/{viruses_db_aa}"
        + "/{smp}"
        + suffix,
        smp=sample_label_read,
        indata=wildcards.indata,
        viruses_db_aa=wildcards.viruses_db_aa,
    )
    return files


def read_viral_tables(file, suffix):
    if is_non_zero_file(file):
        smp = remove_suffix(pathlib.Path(file).stem, suffix)
        df = pd.read_csv(file, sep="\t")
        df["label"] = sample_label_dict_read[smp]
        return df


# def read_tables(file):
#     if is_non_zero_file(file):
#         df = pd.read_csv(file, sep=",")
#         return df


rule mmseqs_profile_viruses_aa_summary:
    input:
        group_tables=lambda wc: get_viral_files(
            wc, suffix=".profile.viral_group-abundances-agg.tsv.gz"
        ),
        hits_tables=lambda wc: get_viral_files(
            wc, suffix=".profile.viral_group-abundances.tsv.gz"
        ),
    output:
        group_abundance_summary=(
            config["rdir"]
            + "/taxonomic-profiling-filtering-aa/"
            + "{indata}"
            + "/{viruses_db_aa}/{viruses_db_aa}-profiling.group-abundances-agg.tsv.gz"
        ),
        abundance_summary=(
            config["rdir"]
            + "/taxonomic-profiling-filtering-aa/"
            + "{indata}"
            + "/{viruses_db_aa}/{viruses_db_aa}-profiling.group-abundances.tsv.gz"
        ),
    threads: 1
    log:
        config["rdir"]
        + "/taxonomic-profiling-filtering-aa/"
        + "{indata}"
        + "/{viruses_db_aa}"
        + "/logs/taxonomic-profiling-filtering-aa/viral.summary.log",
    benchmark:
        config["rdir"]
        +"/taxonomic-profiling-filtering-aa/"
        +"{indata}"
        +"/{viruses_db_aa}"
        +"/logs/taxonomic-profiling-filtering-aa/viral.summary.bmk"
    message:
        """--- Summarize dbcan results."""
    run:
        df = pd.concat(
            map(
                lambda file: read_viral_tables(
                    file, suffix=".profile.viral_group-abundances-agg.tsv"
                ),
                fast_flatten(list({input.group_tables})),
            )
        )
        df.to_csv(
            output.group_abundance_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )

        df = pd.concat(
            map(
                lambda file: read_viral_tables(
                    file, suffix=".profile.viral_group-abundances.tsv"
                ),
                fast_flatten(list({input.hits_tables})),
            )
        )
        df.to_csv(
            output.abundance_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )