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
    version='1.2',
    author='Xiaozhe Ding',
    author_email='xding@caltech.edu',
    description='Scripts needed for for APPARISE.',
    long_description=long_description,
    long_description_content_type='ext/markdown',
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        'Operating System :: MacOS :: MacOS X',
    ),
)
