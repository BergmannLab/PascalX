FROM ubuntu:20.04
ENV DEBIAN_FRONTEND="noninteractive"

# Install system dependencies
# Added: gfortran to provide quadmath.h
RUN apt-get update && apt-get install -y \
    python3 \
    python3-dev \
    python3-setuptools \
    python3-pip \
    g++ \
    gfortran \
    make \
    libboost-all-dev \
    wget \
    unzip
# Add PascalX build output to ldconfig
RUN mkdir -p /PascalX/build/lib
RUN echo "/PascalX/build/lib" > /etc/ld.so.conf.d/pascalx.conf
# Copy PascalX source code
COPY . /PascalX
# Build
RUN cd /PascalX && make all && ldconfig && make test
# Python dependencies
# Fixed: Quoted pip arguments to avoid shell errors
# Create and activate a virtualenv, install into that clean environment
RUN python3 -m pip install --upgrade pip setuptools wheel && \
    python3 -m pip install virtualenv && \
    virtualenv /venv && \
    /venv/bin/pip install --no-cache-dir --prefer-binary \
        --index-url https://pypi.org/simple \
        "numpy==1.24.4" \
        "scipy==1.10.1" \
        "seaborn==0.13.2" \
        "matplotlib==3.7.5" \
        "tqdm==4.67.1" \
        "progressbar2==4.5.0" \
        "python-utils==3.8.2" \
        "sortedcontainers==2.4.0" \
        "pandas==2.0.3"
# Python install
# Changed: use pip install instead of deprecated setup.py install
RUN cd /PascalX/python/ && pip3 install .
# Optional: Install jupyter (uncomment if needed)
RUN pip3 install jupyter
