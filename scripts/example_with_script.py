import subprocess as sp


args = ["samtools", "--help"]
res = sp.check_output(args)
with open(snakemake.output[0], "wb") as f:
    f.write(res)

with open(snakemake.log[0], "w") as f:
    f.write("File created!")
