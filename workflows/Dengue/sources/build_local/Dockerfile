FROM condaforge/mambaforge:23.3.1-1 as conda

# Install snakemake
RUN mamba install -y snakemake-minimal==7.32 -c bioconda

# Pass UID to container (ensures appropriate permissions for created files/directories)
ARG HOST_UID
RUN useradd --uid ${HOST_UID} --create-home user
USER user

# Copy workflow files to container
COPY --chown=user . /home/user/
WORKDIR /home/user/

# Set-up the workflow conda environments
ENV TERM=xterm-color
RUN snakemake --conda-create-envs-only --use-conda --cores 1 extractphylotree__test plotexportsandimports__test

# Run the workflow script
CMD ["bash", "-i", "-c", "./run.sh"]
