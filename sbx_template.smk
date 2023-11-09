try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")


localrules:
    all_template,


rule all_template:
    input:
        QC_FP / "mush" / "big_file.txt",


rule example_rule:
    """Takes in cleaned .fastq.gz and mushes them all together into a file"""
    input:
        expand(QC_FP / "cleaned" / "{sample}_{rp}.fastq.gz", sample=Samples, rp=Pairs),
    output:
        QC_FP / "mush" / "big_file1.txt",
    log:
        LOG_FP / "example_rule.log",
    benchmark:
        BENCHMARK_FP / "example_rule.tsv"
    params:
        opts=Cfg["sbx_template"]["example_rule_options"],
    conda:
        "envs/sbx_template_env.yml"
    shell:
        "cat {params.opts} {input} >> {output} 2> {log}"


rule example_with_script:
    """Take in big_file1 and copy it to big_file using a python script"""
    input:
        QC_FP / "mush" / "big_file1.txt",
    output:
        QC_FP / "mush" / "big_file.txt",
    log:
        LOG_FP / "example_with_script.log",
    benchmark:
        BENCHMARK_FP / "example_with_script.tsv"
    conda:
        "envs/sbx_template_env.yml"
    script:
        "scripts/example_with_script.py"
