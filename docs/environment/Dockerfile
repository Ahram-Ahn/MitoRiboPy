# MitoRiboPy v0.3.0 Docker image.
#
# Build:
#   docker build -t mitoribopy:0.3.0 -f docs/environment/Dockerfile .
#
# Run:
#   docker run --rm -v $PWD:/work -w /work mitoribopy:0.3.0 \
#     mitoribopy align --help
#
# The image is based on the miniconda3 image so we can pull the
# bioinformatics stack (cutadapt, bowtie2, samtools, umi_tools) from
# bioconda in one shot, exactly matching docs/environment/environment.yml.
FROM continuumio/miniconda3:latest

# Avoid interactive prompts from apt; the only system package we need
# is a small set of build-time dependencies for wheels that might fall
# back to source builds.
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        git \
 && rm -rf /var/lib/apt/lists/*

# Create the environment from the repo's environment.yml. We copy only
# the env yml first so that rebuilds triggered by source changes do not
# invalidate this slow layer.
WORKDIR /opt/mitoribopy
COPY docs/environment/environment.yml /opt/mitoribopy/environment.yml
RUN conda env create -f /opt/mitoribopy/environment.yml \
 && conda clean -afy

# Install MitoRiboPy into the env. The 'pip install -e ..' line in
# environment.yml points at '..' relative to the env yml; we replicate
# that here after copying the full source tree.
COPY . /opt/mitoribopy
RUN conda run -n mitoribopy pip install --no-deps -e /opt/mitoribopy

# Make 'mitoribopy' available without activating the env manually.
ENV PATH=/opt/conda/envs/mitoribopy/bin:$PATH

# Default work directory for user mounts.
WORKDIR /work

# Smoke-test the CLI during the build so a broken image fails loudly.
RUN mitoribopy --version

ENTRYPOINT ["mitoribopy"]
CMD ["--help"]
