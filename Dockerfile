# For more information, please refer to https://aka.ms/vscode-docker-python
FROM python:3-slim

LABEL name="libudalab/crossovers"
LABEL maintainer="conery@uoregon.edu"

# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE=1

# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED=1

# Define default paths to data files
ENV XO_SNPS=/data/BSP_TIGER.marker_dataframe.pickle.gzip
ENV XO_INTERVALS=/data/BSP_TIGER.intervals_dataframe.pickle.gzip
ENV XO_CROSSOVERS=/data/BSP_COs_final_set.pickle.gzip
ENV XO_PEAKS=/data/peaks.csv
ENV XO_SAVE=/data/summary.csv

COPY requirements.txt .
RUN python -m pip install -r requirements.txt

WORKDIR /app
COPY . /app

RUN python -m pip install .

# Creates a non-root user with an explicit UID and adds permission to access the /app folder
RUN adduser -u 5678 --disabled-password --gecos "" appuser && chown -R appuser /app
USER appuser

EXPOSE 5006

CMD ["xo", "gui"]
