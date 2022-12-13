import shutil
import sys

sys.stderr(f"Copying {snakemake.input} to {snakemake.output}...")
shutil.copyfile(snakemake.input, snakemake.output)
with open(snakemake.log, "w") as f:
    f.write("File copied!")
