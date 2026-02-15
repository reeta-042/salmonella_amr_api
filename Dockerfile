FROM salmonella_amr:latest

WORKDIR /app

# Set all environment variables
ENV WORK_DIR=/app/work \
    MODELS_DIR=/app/models \
    SCRIPTS_DIR=/app/scripts \
    FEATURE_TEMPLATES_DIR=/app/feature_templates \
    REFERENCE_GENOME=/app/reference/salmonella_LT2.gbff \
    CARD_PROTEIN_FILE=/app/card_db/card_all_proteins.fasta

RUN conda run -n amr_project pip install --no-cache-dir --default-timeout=1000 \
    numpy==1.26.4 \
    fastapi==0.115.0 \
    pandas==2.2.2 \
    scikit-learn==1.4.2 \
    uvicorn==0.32.0 \
    uvloop==0.22.1 \
    httptools==0.7.1 \
    python-multipart==0.0.12

# 2. THE REQUIREMENTS LAYER:
COPY requirements.txt /app/
RUN conda run -n amr_project pip install --no-cache-dir -r requirements.txt

# 3. THE CODE LAYER:
COPY scripts/ /app/scripts/
COPY models/ /app/models/
COPY feature_templates/ /app/feature_templates/
COPY reference/ /app/reference/
COPY card_db/ /app/card_db/
COPY api/ /app/api/

# 4. Setup permissions and directories
RUN mkdir -p /app/logs /app/work /app/uploads /app/results && \
    chmod +x /app/scripts/*.sh

EXPOSE 8000

# 5. Launch the API
CMD ["conda", "run", "--no-capture-output", "-n", "amr_project", \
     "python", "-m", "uvicorn", "api.app:app", "--host", "0.0.0.0", "--port", "8000"]