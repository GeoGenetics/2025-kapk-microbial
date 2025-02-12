from itertools import chain
import os, gzip


def get_prefilter_coarse_inputs_summary(wildcards):
    files = expand(
        config["rdir"]
        + "/prefilter-profiling-coarse/{smp}.prefilter-coarse.tax-summary.tsv.gz",
        smp=sample_label_read,
    )
    return files


def get_prefilter_coarse_inputs(wildcards):
    files = expand(
        config["rdir"]
        + "/prefilter-profiling-coarse/{smp}.prefilter-coarse.tax.tsv.gz",
        smp=sample_label_read,
    )
    return files


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


rule prefilter_profiling_coarse_summary:
    input:
        k2_out_summ=get_prefilter_coarse_inputs_summary,
        k2_out=get_prefilter_coarse_inputs,
    output:
        k2_coarse_summary=(
            config["rdir"]
            + "/prefilter-profiling-coarse/prefilter-coarse-summary.tsv.gz"
        ),
        k2_coarse_tax_summary=(
            config["rdir"]
            + "/prefilter-profiling-coarse/prefilter-coarse-tax-summary.tsv.gz"
        ),
    threads: 1
    log:
        config["rdir"] + "/logs/prefilter-profiling-coarse/prefilter-coarse-summary.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/prefilter-profiling-coarse/prefilter-coarse-summary.bmk"
        )
    message:
        """--- Summarize coarse taxonomic profiling with prefilter."""
    run:
        df = pd.concat(
            map(
                lambda file: pd.read_csv(file, sep="\t"),
                fast_flatten(list({input.k2_out_summ})),
            )
        )
        df["label"] = df["label"].astype(str).map(sample_label_dict_read)
        df.to_csv(
            output.k2_coarse_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )
        if os.path.exists(output.k2_coarse_tax_summary):
            os.remove(output.k2_coarse_tax_summary)
        files = fast_flatten(list({input.k2_out}))
        # get first element of list
        file = files.pop(0)
        df = pd.read_csv(
            file, sep="\t", names=["label", "read", "tax_id", "tax_string"]
        )
        df["label"] = df["label"].astype(str).map(sample_label_dict_read)
        df.to_csv(
            output.k2_coarse_tax_summary,
            sep="\t",
            compression="gzip",
            header=True,
            index=False,
        )

        for file in files:
            df = pd.read_csv(
                file, sep="\t", names=["label", "read", "tax_id", "tax_string"]
            )
            df["label"] = df["label"].astype(str).map(sample_label_dict_read)
            print(df.head())
            with gzip.open(output.k2_coarse_tax_summary, "a") as f:
                f.write(
                    df.to_csv(
                        sep="\t",
                        header=False,
                        index=False,
                    ).encode("utf-8")
                )
