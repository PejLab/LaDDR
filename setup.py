from setuptools import setup, find_packages

setup(
    name="laddr",
    version="0.3.0",
    description="Latent Data-Driven RNA phenotyping",
    author="Daniel Munro",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    entry_points={
        "console_scripts": [
            "laddr=laddr.cli:cli",
        ],
    },
    install_requires=[
        "bx-python",
        "numpy<2",
        "pandas",
        "pyBigWig",
        "PyYAML",
        "scikit-fda",
        "scikit-learn",
        "statsmodels",
        "tqdm",
    ],
    package_data={
        "laddr.resources": [
            "config.default.yaml",
            "config.example.yaml",
            "config.extended.yaml",
            "coverage_manifest.tsv",
            "Snakefile",
            "run.sh",
        ],
    },
)
