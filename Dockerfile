# Minimal, relies on prebuilt manylinux wheels for numpy/scipy
FROM python:3.11-slim

# Make Python faster and quieter in containers
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

COPY requirements.txt .
RUN pip install --upgrade pip && pip install -r requirements.txt

COPY . .
EXPOSE 8080
CMD ["python", "main.py"]

