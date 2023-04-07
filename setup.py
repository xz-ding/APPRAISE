"""
Setup file for APPRAISE package.
Author: Xiaozhe Ding
Email: xding@caltech.edu, dingxiaozhe@gmail.com
"""

import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name='appraise',
    version='1.2.4',
    author='Xiaozhe Ding',
    author_email='xding@caltech.edu',
    description='Scripts and functions needed for for APPARISE.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX'
    ],
    install_requires=[
          'scipy>=1.4.1',
          'numpy>=1.18.2',
          'pandas>=1.1.5',
          'matplotlib>=3.2.1',
          'seaborn>=0.11.2'
    ])
