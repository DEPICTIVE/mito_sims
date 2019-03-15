import os
from setuptools import setup, find_packages

NAME='mito_sims_py3'
ROOT = os.path.abspath(os.path.dirname(__file__))

DESCRIPTION = "A tool for simulating apoptosis dynamics."

try:
    with open(os.path.join(ROOT, 'README.md'), encoding='utf-8') as fid:
        LONG_DESCRIPT = fid.read()
except FileNotFoundError:
    LONG_DESCRIPT = ''
try:
    with open(os.path.join(ROOT, NAME,'__version__'), 'r') as fid:
        VERSION=fid.read().strip()
except FileNotFoundError:
    VERSION='0.0.0error'

setup(name=NAME,
    version=VERSION,
    python_requires='>=3.6.0',
    long_description=LONG_DESCRIPT,
    long_description_content_type='text/markdown',
    description = DESCRIPTION,
    packages=find_packages(exclude=('pars',)),
    package_data={
        '':['__version__']
    },
    setup_requires=[
        'numpy',
        'scipy',
        'matplotlib'],
    install_requires=[
        'matplotlib',
        'scipy',
        'numpy'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
        'Intended Audience :: Developer',
        'Intended Audience :: End Users',
        'Topic :: biological data analysis'
    ]
)
# 'Development Status :: Alpha',
