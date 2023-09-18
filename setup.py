from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="fenics_arclength",
    version="0.1.0",
    author="Peerasait Prachaseree",
    author_email="pprachas@bu.edu",
    description="A numerical arclength solver written on top of FEniCS",
    long_description=long_description,
    long_description_content_type="text/markdown",
    project_urls={
        "Source": "https://github.com/pprachas/fenics_arclength",
        "Documentation": "https://fenics-arclength.readthedocs.io/en/latest/",
    },
    maintainer=[
        ("Peerasait Prachaseree", "pprachas@bu.edu"),
        ("Saeed Mohammadzadeh", "saeedmhz@bu.edu"),
        ("Emma Lejeune", "elejeune@bu.edu"),
    ],
    packages=["arc_length"],
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
    ]
)