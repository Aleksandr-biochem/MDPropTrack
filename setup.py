from setuptools import setup

setup(
    name="MDPropTrack",
    version="0.1.0",
    description="MDPropTrack is a Python3 mini-library for quick property extraction from Molecular Dynamics trajectories, time series plotting and assessment of convergence.",
    url="https://github.com/Aleksandr-biochem/MDPropTrack/",
    author="Aleksandr Kovalenko",
    author_email="alekskov1102@gmail.com",
    license="MIT Licenses",
    packages=["MDPropTrack"],
    install_requires=[
        "numpy<2.0.0",
        "pandas>=2.0.0",
        "MDAnalysis>=2.7.0",
        "panedr>=0.8.0",
        "seaborn>=0.13.2",
        "matplotlib>=3.10.0",
        "tqdm>=4.67.1"
    ],
    entry_points={
        "console_scripts": [
            "mdpt = MDPropTrack:MDPropTrack.main",
        ]
    },
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
    ],
)