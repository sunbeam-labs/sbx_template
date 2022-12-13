import shutil

shutil.copyfile(snakemake.input, snakemake.output)
with open(snakemake.log, "w") as f:
    f.write("File copied!")
