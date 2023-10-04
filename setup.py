import os

from setuptools import find_packages, setup

if os.path.exists("requirements.txt"):
    with open('requirements.txt') as f:
        requires = f.read().splitlines()
else:
    requires = []

setup(
    name="gogstools",
    version='0.1',
    description="GenOuest tools for manipulating Official Gene Sets",
    long_description="GenOuest tools for manipulating Official Gene Sets",
    author="Anthony Bretaudeau",
    author_email="anthony.bretaudeau@inrae.fr",
    url="https://github.com/genouest/ogs-tools",
    install_requires=requires,
    packages=find_packages(),
    license='GPLv3',
    platforms="Posix; MacOS X; Windows",
    entry_points="",
    scripts=[
        'gff2embl/gff2embl',
        'ogs_check/ogs_check',
        'ogs_merge/ogs_merge',
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python :: 3.9",
    ]
)
