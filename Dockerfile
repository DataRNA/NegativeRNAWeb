FROM rocker/shiny:4.4.1

# System dependencies for RNAfold and R packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        vienna-rna \
        libgsl-dev \
        libcurl4-openssl-dev \
        libxml2-dev \
        libssl-dev && \
    rm -rf /var/lib/apt/lists/*

# Verify RNAfold is installed and accessible
RUN which RNAfold && RNAfold --version

# Install required R packages (CRAN + Bioconductor)
RUN R -e "install.packages(c('shiny','shinydashboard','DT','seqinr','fmsb','future','future.apply','progressr','RColorBrewer','openxlsx','plotly','gridExtra','ggplot2','BiocManager'), repos='https://cloud.r-project.org')" && \
    R -e "BiocManager::install('Biostrings', ask = FALSE, update = FALSE)"

# Copy application
WORKDIR /srv/shiny-server/NeRNAWeb
COPY . /srv/shiny-server/NeRNAWeb

# Expose Shiny port
EXPOSE 3838

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('.', host='0.0.0.0', port=3838)"]

