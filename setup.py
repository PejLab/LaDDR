from setuptools import setup, find_packages

setup(
    name="latent-rna",
    version="0.1.0",
    description="Extract latent transcriptomic phenotypes",
    author="Daniel Munro",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    entry_points={
        "console_scripts": [
            "latent-rna=latent_rna.cli:cli",
        ],
    },
    install_requires=[
        "bx-python",
        "gtfparse==2.5.0",
        "numpy<2",
        "pandas",
        "pyBigWig",
        "scikit-fda",
        "scikit-learn",
        "statsmodels",
        "tqdm",
        "PyYAML",
    ],
    package_data={
        "latent_rna": ["config.default.yaml"],
    },
)
