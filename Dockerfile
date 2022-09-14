# Start from Ubuntu
FROM ubuntu:20.04

# Copy over PascalX
COPY . /PascalX

# Install dependencies
RUN mkdir -p /PascalX/build/lib
ENV DEBIAN_FRONTEND="noninteractive"
RUN apt-get update && apt-get install -y python3 python3-dev python3-setuptools python3-pip g++ make libboost-all-dev wget
RUN echo "/PascalX/build/lib" > /etc/ld.so.conf.d/pascalx.conf

# Build
RUN cd /PascalX && make all && ldconfig && make test
RUN cd /PascalX/python/ && python3 setup.py install

# Install jupyter
RUN pip3 install jupyter

