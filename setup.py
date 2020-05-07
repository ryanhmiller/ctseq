from setuptools import find_packages, setup
from .version import __version__

# find tutorial for this here: https://packaging.python.org/tutorials/packaging-projects/

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as f:
    requires = f.read().splitlines()

setup(
    name="ctseq",
    version=__version__,
    description="pipeline to analyze ctDNA data",
    long_description=long_description,
    long_description_content_type='text/markdown',
    author="Ryan Miller",
    author_email="miller.ryan.h@gmail.com",
    url="https://github.com/ryanhmiller/ctseq",
    packages=['ctseq'],
    package_data={"": ["LICENSE", "README.md"]},
    #data_files=[("samplot", ["samplot/templates/samplot_vcf.html"])],
    #include_package_data=True,
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
