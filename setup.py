from setuptools import setup

setup(
    name="mavedbconvert",
    version="0.6.0-alpha",
    packages=["mavedbconvert", "mavedbconvert.tests"],
    url="https://github.com/FowlerLab/mavedb-convert",
    license="AGPLv3",
    author="Daniel Esposito",
    author_email="esposito.d@wehi.edu.au",
    description=(
        "A command line tool for converting alternate "
        "file formats into a MaveDB compliant format."
    ),
    #install_requires=open("requirements/install.txt", "rt").read().split("\n"),
    entry_points={"console_scripts": ["mavedb-convert=mavedbconvert.main:main"]},
)
