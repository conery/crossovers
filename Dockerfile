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

# Define default params for postprocessor
ENV XO_POST_BLOCK_SIZE='0 100'
ENV XO_POST_BLOCK_LENGTH='0 1000'
ENV XO_POST_COVERAGE=2
ENV XO_POST_MATCH=True
ENV XO_POST_HIGH_Z=0.9
ENV XO_DELTA_HIGH_Z=0.1

WORKDIR /app
COPY . /app

RUN python -m pip install .

# Creates a non-root user with an explicit UID and adds permission to access the /app folder
RUN adduser -u 5678 --disabled-password --gecos "" appuser && chown -R appuser /app
USER appuser

EXPOSE 5006

CMD ["xo", "gui"]
