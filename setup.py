import setuptools
import sys

with open("README.md", "r") as fh:
    long_description = fh.read()

requirements = [
    "tables",
    "pandas",
    "openpyxl",
    "tqdm",
    "docopt",
    "hgvsp @ git+https://github.com/VariantEffect/hgvs-patterns",
    "requests",
    "numpy",
    "scipy",
    "joblib",
    "xlsxwriter",
    "fqfa>=1.1.0",
]
# fqfa requires backported dataclasses in Python 3.6
if sys.version_info.major == 3 and sys.version_info.minor == 6:
    requirements.append("dataclasses")

setuptools.setup(
    name="mavedbconvert",
    version="0.1.1-beta",
    author="Alan F Rubin, Daniel Esposito",
    author_email="alan.rubin@wehi.edu.au",
    description=(
        "A command line tool for converting Multiplex Assay of Variant Effect datasets into a MaveDB-ready format."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/VariantEffect/mavedbconvert",
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=requirements,
    entry_points={"console_scripts": ["mavedbconvert=mavedbconvert.main:main"]},
    test_suite="tests",
)
