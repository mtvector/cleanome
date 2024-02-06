
from setuptools import setup, find_packages

setup(
    name="cleanome",
    author='Matthew Schmitz',
    version="0.1",
    packages=find_packages(),
    python_requires='>=3.6',
    install_requires=[
        'gtfparse',
        'biopython',
        'pandas',
        'numpy',
        'tqdm',
        'requests',
        'ete3'
]
    ],
)
