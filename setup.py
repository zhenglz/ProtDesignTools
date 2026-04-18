
from setuptools import setup, find_packages

setup(
    name="protdesigntools",
    version="0.1.0",
    description="A Modular Protein Design Toolkit",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
    ],
    entry_points={
        "console_scripts": [
            "protdesign=protdesigntools.main:main",
        ]
    }
)
