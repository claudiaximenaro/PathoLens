from setuptools import setup, find_packages

setup(
    name="patholens",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "biopython",
        "openpyxl",
        "requests",
    ],
    python_requires=">=3.9",
)