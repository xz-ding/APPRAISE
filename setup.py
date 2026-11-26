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
    version='1.3.2',
    author='Xiaozhe Ding',
    author_email='xding@caltech.edu',
    description='Protein competition analysis with APPRAISE.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=setuptools.find_packages(),
    package_data={'appraise': ['examples/b2ar/*.gz', 'examples/b2ar/panel.json', 'examples/b2ar/SOURCE.md']},
    python_requires='>=3.10',
    extras_require={'structures': ['pymol-open-source>=3.1.0a0']},
    classifiers=[
        "Programming Language :: Python :: 3",
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX'
    ],
    install_requires=[
          'scipy>=1.10',
          'numpy>=1.23.5',
          'pandas>=1.5',
          'matplotlib>=3.6',
          'seaborn>=0.11.2'
    ])
