# APPRAISE

***A***utomated ***P***air-wise ***P***eptide-***R***eceptor binding model ***A***nalys***I***s for ***S***creening ***E***ngineered proteins (***APPRAISE***) is a method that predicts the receptor binding propensity of engineered proteins based on high-precision protein structure prediction tools, such as AlphaFold2-multimer. This Python package includes tools for preparing input files and analyzing the modeled structures.

Author: Xiaozhe Ding
Email: dingxiaozhe@gmail.com, xding@caltech.edu
Twitter: @DingXiaozhe

## Environment

APPRAISE 1.2 was tested with the following environment:

 - MacOS 10.14.6

 - Python 3.6.10

 - Alphafold-colabfold 2.1.14 (Available [here](https://github.com/sokrypton/ColabFold))

 - PyMOL 2.3.3 (Schrodinger LLC.)

 - Python packages:

    - scipy 1.4.1

    - numpy 1.18.2

    - pandas 1.1.5

    - matplotlib 3.2.1

    - seaborn 0.11.2



## Installation

### Option 1 (recommended)
Install the distribution from PyPI. In the terminal, run:

```
pip install appraise
```

### Option 2 (back-up)
Download the repository to your local computer and unzip. In the terminal, [change the working folder](https://ss64.com/osx/cd.html) to the directory containing the appraise package folder, and run the following line:

```
pip install -e appraise
```

## Get started

The demo jupyter notebook (./demo/appraise_demo.ipynb) will serve as a guide to help you go through the APPRAISE workflow.


## References

Github repository for APPRAISE: https://github.com/xz-ding/APPRAISE
(Contains the latest version of APPRAISE package and demo notebooks.)

Github repository for ColabFold: https://github.com/sokrypton/ColabFold
(ColabFold provides a panel of user-friendly tools for structure modeling.)

The APPRAISE manuscript: Pending
