FROM rocker/r-ver:4.3.0

# Install system dependencies needed by TCGAbiolinks, maftools, and other
# Bioconductor/CRAN packages used in this pipeline
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    default-jdk \
    && rm -rf /var/lib/apt/lists/*

# Install renv and restore the project's exact R package environment
RUN R -e "install.packages('renv')"

WORKDIR /project
COPY renv.lock renv.lock
RUN R -e "renv::restore()"

COPY . .

CMD ["R"]