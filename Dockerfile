# Use Python 3.8 as base image
FROM python:3.8-slim

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    libxml2-dev \
    libxslt-dev \
    libz-dev \
    ncbi-blast+ \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Create directories for data mounting
RUN mkdir -p /data /config

# Set the config file location
ENV CONFIG_PATH=/config/config.ini

# Default command shows help
CMD ["python", "AddedAnnotations.py", "--help"]
