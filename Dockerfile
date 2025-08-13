# Minimal image to run the 2D FDTD TE simulation in headless mode (frames only by default)

FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    HEADLESS=1

WORKDIR /app

# System deps kept minimal; Agg backend needs no GUI
RUN pip install --no-cache-dir --upgrade pip

COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy project
COPY . /app

RUN mkdir -p /app/outputs
VOLUME ["/app/outputs"]

CMD ["python", "examples/run_simulation.py"]

# Minimal image to run the 2D FDTD TE simulation in headless mode

FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    HEADLESS=1

WORKDIR /app

# System deps
RUN apt-get update \
    && apt-get install -y --no-install-recommends ffmpeg \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache-dir --upgrade pip

COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy project
COPY . /app

# Ensure outputs directory exists
RUN mkdir -p /app/outputs
VOLUME ["/app/outputs"]

# Default command: run the headless example; override HEADLESS=0 to try GUI (not recommended in container)
CMD ["python", "examples/run_simulation.py"]


