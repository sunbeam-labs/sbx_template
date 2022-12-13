import sys

sys.stderr(f"Copying {type(snakemake.input)} to {type(snakemake.output)}...")
sys.stderr(f"{snakemake.input[0]}")
snakemake.output.write(snakemake.input.readlines())
with open(snakemake.log, "w") as f:
    f.write("File copied!")
