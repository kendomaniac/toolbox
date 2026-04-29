FROM ubuntu
RUN apt -y update
RUN apt -y upgrade
RUN apt -y install sudo
RUN apt -y install vim
RUN apt -y install nano
ENV DEBIAN_FRONTEND=noninteractive
RUN sudo apt-get update && apt-get install -y --no-install-recommends r-base r-base-dev build-essential
RUN Rscript -e "install.packages('data.table', repos='https://cloud.r-project.org')"
RUN Rscript -e "install.packages('Matrix', repos='https://cloud.r-project.org')"
RUN Rscript -e "install.packages('ggplot2', repos='https://cloud.r-project.org')"
RUN Rscript -e "install.packages('R.utils', repos='https://cloud.r-project.org')"
COPY Ribomito_v3.6.3.R /home/Ribomito_v3.6.3.R
RUN chmod 644 /home/Ribomito_v3.6.3.R
RUN apt-get update
RUN apt-get install -y --no-install-recommends python3-pip python3-venv 
RUN python3 -m venv /opt/venv
RUN /opt/venv/bin/pip install --no-cache-dir parquet-tools 
RUN ln -s /opt/venv/bin/parquet-tools /usr/local/bin/parquet-tools 
RUN rm -rf /var/lib/apt/lists/*






