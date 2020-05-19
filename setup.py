import os
from setuptools import find_packages, setup

# find tutorial for this here: https://packaging.python.org/tutorials/packaging-projects/

version_py = os.path.join(os.path.dirname(__file__), 'ctseq', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','').strip()

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as f:
    requires = f.read().splitlines()

setup(
    name="ctseq",
    version=version,
    description="pipeline to analyze ctDNA data",
    long_description=long_description,
    long_description_content_type='text/markdown',
    author="Ryan Miller",
    author_email="miller.ryan.h@gmail.com",
    url="https://github.com/ryanhmiller/ctseq",
    packages=['ctseq'],
    package_data={"": ["LICENSE", "README.md","ctseq/*.R","ctseq/*.r"]},
    #data_files=[("samplot", ["samplot/templates/samplot_vcf.html"])],
    include_package_data=True,
    install_requires=requires,
    license="MIT",
    zip_safe=False,
    entry_points={"console_scripts": ["ctseq = ctseq.__main__:main"]},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.7'
)
