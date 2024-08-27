# For more information, please refer to https://aka.ms/vscode-docker-python
FROM python:3-slim

# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE=1

# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED=1

# Define default paths to data files
ENV XO_SNPS=/data/BSP_TIGER.marker_dataframe.pickle.gzip
ENV XO_INTERVALS=/data/BSP_TIGER.intervals_dataframe.pickle.gzip
ENV XO_PEAKS=/data/peaks.csv

WORKDIR /app
COPY . /app

# Install using pyproject.toml (copied in with the rest of the app)
RUN python -m pip install .

# Creates a non-root user with an explicit UID and adds permission to access the /app folder
RUN adduser -u 5678 --disabled-password --gecos "" appuser && chown -R appuser /app
USER appuser

EXPOSE 5006

CMD ["xo", "view"]
