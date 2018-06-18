#!/usr/bin/env python3

import setuptools import setup

## Description will come from README. md file

## Version will come from version.py


with open("README.md", "r") as f:
    long_descr = f.read()


version = {}
with open('iliner/version.py') as fv:
    exec(fv.read(), version)


setup(name='iliner',
    packages=['iliner'],
    version=version['__version__'],
    description='Run iLASH algorithm'
    long_description=long_descr,
    entry_points={
    "console_scripts": ['iliner=iliner.iliner:main']
    },
    url='http://github.com/hpoisner/iliner',
    author='Hannah Poisner',
    author_email='hannah.poisner@icahn.mssm.edu',
    install_requires=[
    'pandas', 'numpy',
    'matplotlib', 'seaborn'
    ],
    python_requires='>=3',
    package_dir={'iliner': 'iliner'},
    keywords=['iLASH', 'IBD', 'identity-by-descent', 'genotype'],
    classifiers=[
    ## Classifiers defined in http://pypi.python.org/pypi?%3Aaction=list_classifiers
    # 'Development Status :: 1 - Planning',
    # 'Development Status :: 2 - Pre-Alpha',
    'Development Status :: 3 - Alpha',
    # 'Development Status :: 4 - Beta',
    # 'Development Status :: 5 - Production/Stable',
    # 'Development Status :: 6 - Mature',
    # 'Development Status :: 7 - Inactive'
    'Intended Audience :: Science/Research',
    'License :: (TBD),
    'Programming Language :: Python :: 3 :: Only',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Medical Science Apps.'],
    license=tbd)


"""

Sources:

- https://github.com/frichter/ore/blob/master/setup.py
- https://packaging.python.org/tutorials/distributing-packages/
"""
