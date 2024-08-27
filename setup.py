from setuptools import setup, find_packages

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
    name="cleanome",
    author='Matthew Schmitz',
    description="Downloading genome annotations and turning them into genomics reference annotations",
    long_description=long_description,
    long_description_content_type='text/markdown',  # Ensure your README is displayed correctly on PyPI
    version="0.1",
    packages=find_packages(),
    python_requires='>=3.6',
    # install_requires=[#If uncommented, pip might try to do something silly
    #     'gtfparse',
    #     'polars-lts-cpu',
    #     'biopython',
    #     'pandas',
    #     'numpy',
    #     'tqdm',
    #     'requests',
    #     'ete3'
    # ],
    entry_points={
        'console_scripts': [
            'download_genomes=cleanome.download_genomes:main',
            'get_genomes_and_stats=cleanome.get_genomes_and_stats:main',
            'make_cellranger_arc_sh=cleanome.make_cellranger_arc_sh:main',
            'debug_gtf=cleanome.debug_gtf:main',
        ],
    },
)