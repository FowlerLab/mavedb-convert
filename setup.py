from setuptools import setup

setup(
    name='mavedbconvert',
    version='v0.1-alpha',
    packages=[
        'mavedbconvert',
        'mavedbconvert.tests',
    ],
    url='https://github.com/FowlerLab/mavedb-convert',
    license='AGPLv3',
    author='Daniel Esposito',
    author_email='esposito.d@wehi.edu.au',
    description='A command line tool for converting alternate '
                'file formats into a MaveDB compliant format.',
    entry_points={
        'console_scripts': [
            'mavedb-convert=mavedbconvert.main:main',
        ],
    },
)
