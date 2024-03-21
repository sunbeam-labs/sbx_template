FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_template_env

COPY envs/sbx_template_env.yml ./

# Install environment
RUN conda env create --file sbx_template_env.yml --name sbx_template

ENV PATH="/opt/conda/envs/sbx_template/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "sbx_template", "/bin/bash", "-c"]

# Run
CMD "bash"