from setuptools import setup, find_packages

with open("biscrosslinker/README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="biscrosslinker",
    version="0.1.0",
    author="Your Name",
    description="A package for simulating crosslinked polymer networks",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "MDAnalysis>=2.0.0",
        "numpy>=1.20.0",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)