from setuptools import setup, find_packages

with open("README.md") as f:
    LONG_DESC = f.read()

setup(
    name="pylsmc",
    version=0.1,
    description="Lattice Switch Monte Carlo simulations in Python",
    author="Joe Marsh Rossney",
    url="https://github.com/marshrossney/pylsmc",
    long_description=LONG_DESC,
    package=find_packages(),
    #entry_points={
    #    "console_scripts": [
    #    ]
    },
)
