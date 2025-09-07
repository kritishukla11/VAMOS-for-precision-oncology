from setuptools import setup, find_packages

setup(
    name="vamos-precision",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy","pandas","scikit-learn","matplotlib","seaborn","scipy","anndata","scanpy","pyyaml"
    ],
    python_requires=">=3.9",
    entry_points={"console_scripts": ["vamos-precision=vamos_precision.cli:main"]},
)
