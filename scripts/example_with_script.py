import sys

#sys.stderr(f"Copying {snakemake.input} to {snakemake.output}...")
snakemake.output.write(snakemake.input.readlines())
with open(snakemake.log, "w") as f:
    f.write("File copied!")
