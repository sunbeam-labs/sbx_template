from pathlib import Path

Path(snakemake.output[0]).symlink_to(snakemake.input[0])
with open(snakemake.log, "w") as f:
    f.write("File copied!")
