try:
    SBX_TEMPLATE_VERSION = get_ext_version("sbx_template")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_TEMPLATE_VERSION = "0.0.0"

try:
    logger = get_extension_logger("sbx_template")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    import logging

    logger = logging.getLogger("sunbeam.pipeline.extensions.sbx_template")


logger.info("Doing some extension specific setup...")
logger.info(f"Using sbx_template version {SBX_TEMPLATE_VERSION}.")
logger.error("Don't worry, this isn't a real error.")


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
    container:
        f"docker://sunbeamlabs/sbx_template:{SBX_TEMPLATE_VERSION}"
    shell:
        "(cat {params.opts} {input} > {output}) > {log} 2>&1"


rule example_with_script:
    """Take in big_file1 and then ignore it and write the results of `samtools --help` to the output using a python script"""
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
    container:
        f"docker://sunbeamlabs/sbx_template:{SBX_TEMPLATE_VERSION}"
    script:
        "scripts/example_with_script.py"
