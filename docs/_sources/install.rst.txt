Installation
============

.. _DockerSec:

Multi-platform (Docker)
-----------------------

Build a `Docker <https://docker.com>`_ image via

::

    git clone https://github.com/BergmannLab/PascalX.git
    cd PascalX
    docker build . -t pascalx:latest

The image can be run in interactive mode with the host directory ``/your/workdir`` mounted as ``/data`` using the command 

::

    docker run --mount src=/your/workdir,target=/data,type=bind -p 8888:8888 -it pascalx bash
    
For data persistance after stopping the image, you should work exlusively in the ``/data`` directory. 

A jupyter notebook running in the image can be started with 

::

    jupyter notebook --ip 0.0.0.0 --allow-root
    
Open on your host machine a browser and visit `http://localhost:8888/ <http://localhost:8888/>`_ to use PascalX.

    
.. note::

    You should allocate sufficient system resources to the docker runtime. A minimum of 8GB RAM should be made available. For making full use of parallel computation capabilities, at least 64GB of RAM is recommended. However, the precise numbers depend on the size of the GWAS and reference panel used.



Debian/Ubuntu
-------------

Use the Docker_ image or compile and install by yourself as follows:

Requirements:
    * Python3 with development headers
    * GNU g++ with libquadmath, make
    * BOOST_ libraries
    
.. _BOOST: https://www.boost.org 

Install of requirements (for sudoers):

::

    sudo apt install python3 python3-dev python3-setuptools python3-pip g++ make libboost-all-dev


Set library path:

::

    export LD_LIBRARY_PATH="/yourpath/PascalX/build/lib:$LD_LIBRARY_PATH"

Or, as sudoer, add to ``/etc/ld.so.conf.d/``.



Install of PascalX:

::

    git clone https://github.com/BergmannLab/PascalX.git
    cd PascalX
    make all
    cd python
    python3 setup.py install



Mac
---

As the default LLVM compiler used by Xcode does not support quadmath at the time being, we recommend to use Docker_ on MacOS.


Windows
-------

Please use the Docker_ image. 


.. toctree:
	:maxdepth: 2

