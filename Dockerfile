FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=UTC

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget curl git build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

ENV PATH="/opt/conda/bin:$PATH"
RUN conda init bash

# Accept conda TOS
RUN conda config --set channel_priority flexible && \
    echo "yes" | conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    echo "yes" | conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Create conda environment and install bioinformatics tools
RUN conda create -n amr_project python=3.9.19 -y && \
    conda install -n amr_project -c bioconda -c conda-forge \
    abricate=1.0.1 blast=2.16.0 snippy=4.6.0 -y

SHELL ["conda", "run", "-n", "amr_project", "/bin/bash", "-c"]

WORKDIR /app

# Environment variables
ENV WORK_DIR=/app/work \
    MODELS_DIR=/app/models \
    SCRIPTS_DIR=/app/scripts \
    FEATURE_TEMPLATES_DIR=/app/feature_templates \
    REFERENCE_GENOME=/app/reference/salmonella_LT2.gbff \
    CARD_PROTEIN_FILE=/app/card_db/card_all_proteins.fasta

# Install core Python dependencies
RUN conda run -n amr_project pip install --no-cache-dir --default-timeout=1000 \
    numpy==1.26.4 \
    fastapi==0.115.0 \
    pandas==2.2.2 \
    scikit-learn==1.4.2 \
    uvicorn==0.32.0 \
    uvloop==0.22.1 \
    httptools==0.7.1 \
    python-multipart==0.0.12

# Install remaining dependencies from requirements.txt
COPY requirements.txt /app/
RUN conda run -n amr_project pip install --no-cache-dir -r requirements.txt

# Cleanup
RUN conda clean --all -f -y && \
    rm -rf /opt/conda/pkgs/* && \
    rm -rf /root/.cache/pip && \
    find /opt/conda -name "*.pyc" -delete && \
    find /opt/conda -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true

# Copy project files
COPY scripts/ /app/scripts/
COPY models/ /app/models/
COPY feature_templates/ /app/feature_templates/
COPY reference/ /app/reference/
COPY card_db/ /app/card_db/
COPY api/ /app/api/

# Setup directories and permissions
RUN mkdir -p /app/logs /app/work /app/uploads /app/results && \
    chmod +x /app/scripts/*.sh

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import requests; requests.get('http://localhost:8000/health')" || exit 1

CMD ["conda", "run", "--no-capture-output", "-n", "amr_project", \
     "python", "-m", "uvicorn", "api.app:app", "--host", "0.0.0.0", "--port", "8000"]
