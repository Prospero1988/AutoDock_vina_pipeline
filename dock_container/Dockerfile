FROM debian:latest

# Step 1: Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Step 2: Set the working directory
WORKDIR /root

# Step 3: Copy data in to /root/dock
COPY dock /root/dock
COPY entrypoint.sh /entrypoint.sh

# Verify the contents of copied directories and files
RUN ls -la /root/dock && ls -la /root/dock/installation
RUN test -f /root/dock/installation/modules.py || (echo "modules.py not found!" && exit 1)
RUN test -f /entrypoint.sh || (echo "Plik nie istnieje" && exit 1)

# Step 3: Install necessary tools
RUN apt-get update && apt-get install -y \
    \
    vim \
    nano \
    wget \
    curl \
    sudo \
    tar \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/

# Step 5: Download and install Miniconda
RUN echo "Downloading Miniconda..." && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    echo "Installing Miniconda..." && \
    bash miniconda.sh -b -p /root/miniconda && \
    rm miniconda.sh

# Step 6: Configure Miniconda
RUN echo "source /root/miniconda/etc/profile.d/conda.sh" >> /root/.bashrc && \
    /root/miniconda/bin/conda init bash

# Step 7: Create a conda environment
RUN /root/miniconda/bin/conda create -y -n auto_dock python=3.11

# Step 8: Run the Python script using the conda environment
RUN /root/miniconda/bin/conda run -n auto_dock python /root/dock/installation/modules.py \
    /root/miniconda/bin/conda clean -afy \
    rm -rf /root/.cache

# Step 9: Install and configure SLURM
RUN apt update && apt upgrade -y 
RUN apt install -y build-essential slurm-wlm slurm-client
RUN usermod -d /var/lib/slurm slurm
RUN usermod -s /bin/false slurm
RUN mkdir -p /var/spool/slurmd /var/log/slurm
RUN chown -R slurm: /var/spool/slurmd /var/log/slurm
RUN chmod -R 755 /var/spool/slurmd /var/log/slurm
COPY slurm.conf /etc/slurm/
RUN test -f /etc/slurm/slurm.conf || (echo "Plik nie istnieje" && exit 1)
RUN chown slurm: /etc/slurm/slurm.conf
RUN chmod 644 /etc/slurm/slurm.conf

# Step 10: Remove non used packages

RUN apt-get remove --purge -y \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# Step XX: Expose the Streamlit port and HTML server port
EXPOSE 8001
EXPOSE 8501

# Step XX: Set the default command to bash
RUN chmod +x /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]

