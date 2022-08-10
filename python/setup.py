#!/usr/bin/python3

# Setup script for PascalX

from setuptools import setup, find_packages, Extension
import os, io, re

os.environ["CC"] = "gcc"

def read(*names, **kwargs):
        with io.open(
            os.path.join(os.path.dirname(__file__), *names),
            encoding=kwargs.get("encoding", "utf8")
        ) as fp:
            return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

setup(
    name='PascalX',
    description='work in progress',
    version=find_version("PascalX", "__init__.py"),
    author='D. Krefl',
    author_email='daniel.krefl@unil.ch',
    packages=find_packages(),
    url='https://github.com/BergmannLab/PascalX',
    license='LICENSE',
    long_description=open('README.md').read(),
    install_requires=[
	"cffi>=1.0.0",
#	"numba>=0.51.2",
	"matplotlib>=3.1.0",
	"sortedcontainers>=2.1.0",
	"tqdm>=4.43.0",
	"scipy>=1.4.0",
	"numpy>=1.18.0",
    "fastnumbers>=3.1.0",
	"seaborn>=0.11.0",
	"progressbar>=2.5",
	"sphinx>=3.2.1",
	"sphinx-rtd-theme>=0.5.0",
	"fastnumbers>=3.1.0"
    ],
    setup_requires=["cffi","path.py"],
    cffi_modules=["hpstats.py:ffibuilder","wchissum.py:ffibuilder"],
    zip_safe=False,
)

