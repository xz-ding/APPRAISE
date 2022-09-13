#!/usr/bin/env python3.9
"""
HullRad
Version 8.1
Calculates hydrodynamic coefficients for a biomolecular structure.
AUTHOR: Patrick Fleming

Release Notes:
Version 5 (The program as described in Fleming and Fleming, BJ, 2018)
    Protein
    Nucleic Acids
Version 6 (2019)
    Includes many saccharides for glycosylation
    Rg is calculated from all atom model, not reduced model as in version 5
    Includes option to use numpy and scipy
    Asphericity is calculated from gyration tensor (requires numpy)
    Option to use either qconvex from qhull or ConvexHull from scipy
        (Many thanks to Chad Brautigam for implementation)
Version 7 (2020)
    Ported to Python 3 (See HullRadV7_3.py, this edition is Python 2)
    Addition of chain ID
    Includes several detergents for analysis of protein/detergent micelles.
        (See DT_data below for list.)
Version 8 (2021)
    Accepts both PDB format and mmCIF format files
    Error checking for presence of all backbone atoms
      (Needed to calculate hydration layer)
Version 8.1 (2022)
    Correct bugs introduced with mmCif reader

If you publish work that uses this code please cite:
Fleming, PJ and Fleming, KG "HullRad: Fast Calculations of Folded and Disordered
Protein and Nucleic Acid Hydrogynamic Properties", Biophysical Journal, 2018,
114:856-869
DOI: 10.1016/j.bpj.2018.01.002


The author would like to acknowledge the work of Jose Garcia de la Torre who led
the way in development of algorithms for calculating macromolecular
hydrodynamic properties.

--------------------------------------------------------------------------------
Output looks like:

 pdb8rat.ent
 HYDROLASE (NUCLEIC ACID,RNA)            13-AUG-91   8RAT
  #Amino Acids    :          124
  #Nucleotides    :            0
  #Saccharides    :            0
  #Detergents     :            0
  M               :        13692     g/mol
  v_bar           :        0.710     mL/g
  R(Anhydrous)    :        15.68     Angstroms
  Axial Ratio     :         1.31
  f/fo            :         1.21
  Dt              :       1.13e-06   cm^2/s
  R(Translation)  :        18.99     Angstroms
  s               :       1.84e-13   sec
  [eta]           :         3.16     cm^3/g
  Dr              :       1.91e+07   s^-1
  R(Rotation)     :        20.36     Angstroms
  tauC            :         8.84     ns (from R_rotation)
  Asphericity     :         0.14     (from Gyration Tensor)

To write out the unified side chain psuedo-atom model uncomment the line
near the bottom labeled, # Write out model (~line 1514)

--------------------------------------------------------------------------------

Uses PDB or mmCIF file format as input.

If your python installation includes numpy and scipy YOU DON'T NEED ANYTHING ELSE.
This is the preferred way to go.

If you don't have numpy and scipy installed, you will need a separate program to
    calculate the convex hull - called qconvex.
This is a program in the qhull suite of programs.

Qhull may be downloaded from http://www.qhull.org/download.
	For UNIX download "Qhull_2015.2 for Unix"
	tar xzvf qhull-2015-src-7.2.0.tgz
	cd qhull-2015.2/
	make
	sudo make install
	and the path for qconvex will be /usr/local/bin/qconvex

OS X qhull binaries are available from MacPorts package manager.
	type "sudo port install qhull".
	and the path for qconvex will be /opt/local/bin/qconvex

OS X binaries are also available from Fink package manager.

The default in this script is the UNIX path (/usr/local/bin/qconvex) but you may
have to change it if you are using OS X. (Note: This is NOT used if you have
numpy and scipy installed.)

After installing qhull go to: ### Edit the path to qconvex ### below (~line 476) and
change /usr/local/bin/qconvex to /opt/local/bin/qconvex if that is where qconvex is.

If you are using OS X but are familiar with UNIX and have gcc installed, the UNIX
instructions above work fine on a Mac.

FOR WINDOWS USERS:
There are precompiled versions of the qhull programs available for Windows 7 and up
at http://www.qhull.org/download
Unpack qhull in a directory of choice and ### Edit the path to qconvex ### below (~line 330)
to include the path where your executable was unpacked.

--------------------------------------------------------------------------------
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org>
"""

USAGE = """
python HullRadV8.1.py [PDBcode].pdb
    or
python HullRadV8.1.py [PDBcode].cif
"""

import os
import sys
import math
import string
import shlex, subprocess
import gzip
import mimetypes
import re

class MMCIFWrapperSyntaxError(Exception):

    def __init__(self, category):
        self.category = category

    def __str__(self):
        return "More items than values for category " + repr(self.category) + "!"


class MMCIF2Dict:
    """
    MMCIF2Dict is a purely algorithmic parser that takes as input public
    mmCIF files and creates a python dictionary from them.

    Because this parser is highly optimised for public mmCIF format, it is
    highly unlikely that it will work successfully on any other formatted mmCIF
    file.

    MMCIF2Dict will not work on mmCIF dictionaries!

    AUTHOR: Glen van Ginkel (Protein Data Bank in Europe; http://pdbe.org)

    """

    loopRE = re.compile(r"^\s*[L|l][O|o][O|o][P|p]_.*")
    # commentsRE = re.compile(r'(.*?)\s#.*$')
    dataRE = re.compile(r"^\s*[D|d][A|a][T|t][A|a]_(?P<data_heading>.*)\s*")
    saveRE = re.compile(r"^\s*[S|s][A|a][V|v][E|e]_(?P<save_heading>.*)\s*")
    dataNameRE = re.compile(r"^\s*(?P<data_category>_[\S]+)(?P<remainder>.*)")
    dataCategoryItem = re.compile(
        r"^\s*(?P<data_category>_[\S]+)(?:\.)(?P<category_item>\S+)"
    )
    dataValueRE = re.compile(r'\s*(\'[\S\s]+?\'(?=\s)|"[\S\s]+?"(?=\s)|[\S]+)', re.M)
    header = ""
    data_map = None
    file_path = None
    reserve_token_order = False

    def parse(
        self,
        file_path,
        ignoreCategories=[],
        preserve_token_order=False,
        onlyCategories=[],
    ):
        """Public method which only functions to check the existence of
        the mmCIF file in preparation for reading in the private parseFile
        method.
        """
        if os.path.exists(file_path) and os.path.isfile(file_path):
            return self._parseFile(
                file_path, ignoreCategories, preserve_token_order, onlyCategories
            )
        else:
            print("The file provided does not exist or is not a file.")
            return None

    def _tokenizeData(self, line):
        """Private method that will do the work of parsing the mmCIF data file
        return Dictionary"""
        if "'" in line or '"' in line:
            return [
                x[1:-1] if x[0] == x[-1] and x[0] in ["'", '"'] else x
                for x in self.dataValueRE.findall(line)
            ]
        else:
            return line.strip().split()

    def _parseFile(
        self, file_path, ignoreCategories, preserve_token_order, onlyCategories
    ):
        """Private method that will do the work of parsing the mmCIF data file
        return Dictionary"""

        if preserve_token_order:
            try:
                from collections import OrderedDict as _dict
            except ImportError:
                # fallback: try to use the ordereddict backport when using python 2.6
                try:
                    from ordereddict import OrderedDict as _dict
                except ImportError:
                    # backport not installed: use local OrderedDict
                    from mmCif.ordereddict import OrderedDict as _dict
        else:
            _dict = dict

        mmcif_like_file = _dict()
        data_block = _dict()
        save_block = _dict()

        data_heading = ""
        line_num = 0
        try:
            with openGzip(file_path, "rt") as f1:
                table_names = []
                table_values = []
                table_values_array = []
                isLoop = False
                multiLineValue = False
                skipCategory = False
                for line in f1:
                    line_num += 1
                    if skipCategory:
                        flag = False
                        while line:
                            check = (
                                line.strip().startswith("_")
                                or self.loopRE.match(line.strip()[:5])
                                or self.saveRE.match(line.strip()[:5])
                                or self.dataRE.match(line.strip()[:5])
                            )
                            if flag:
                                if check:
                                    isLoop = False
                                    break
                            else:
                                if not check:
                                    flag = True
                            if not (
                                self.saveRE.match(line.strip()[:5])
                                or self.dataRE.match(line.strip()[:5])
                            ):
                                try:
                                    line = next(f1)
                                    line_num += 1
                                except StopIteration:
                                    break
                            else:
                                break
                        skipCategory = False

                    if (
                        isLoop is True
                        and table_values_array != []
                        and (
                            self.loopRE.match(line) is not None
                            or (line.strip().startswith("_"))
                        )
                    ):
                        isLoop = False
                        num_item = len(table_names)
                        if len(table_values_array) % num_item != 0:
                            raise MMCIFWrapperSyntaxError(category)
                        for val_index, item in enumerate(table_names):
                            data_block[category][item] = table_values_array[
                                val_index::num_item
                            ]
                        table_values_array = []

                    if line.strip() == "":
                        continue
                    if line.startswith("#"):
                        continue
                    if "\t#" in line or " #" in line and not line.startswith(";"):
                        new_line = ""
                        for tok in self.dataValueRE.findall(line):
                            if not tok.startswith("#"):
                                new_line += tok + " "
                            else:
                                break
                        # make sure to preserve the fact that ';' was not the first character
                        line = (
                            new_line if not new_line.startswith(";") else " " + new_line
                        )
                        # Fails for entries "3snv", "1kmm", "1ser", "2prg", "3oqd"
                        # line = re.sub(r'\s#.*$', '', line)
                    if line.startswith(";"):
                        while "\n;" not in line:
                            try:
                                line += next(f1)
                                line_num += 1
                            except StopIteration:
                                break
                        multiLineValue = True
                    if self.dataRE.match(line):
                        if data_block != {}:
                            if table_values_array != []:
                                isLoop = False
                                num_item = len(table_names)
                                if len(table_values_array) % num_item != 0:
                                    raise mmCifSyntaxError(category)
                                for val_index, item in enumerate(table_names):
                                    data_block[category][item] = table_values_array[
                                        val_index::num_item
                                    ]
                                table_names = []
                                table_values_array = []
                            mmcif_like_file[data_heading] = data_block
                            data_block = _dict()
                        data_heading = self.dataRE.match(line).group("data_heading")
                    elif self.saveRE.match(line):
                        while line.strip() != "save_":
                            try:
                                line = next(f1)
                                line_num += 1
                            except StopIteration:
                                break
                        continue
                    elif self.loopRE.match(line):
                        # Save and clear the table_values_array buffer from the
                        # previous loop that was read
                        if table_values_array != []:
                            for itemIndex, name in enumerate(table_names):
                                data_block[category].update(
                                    {
                                        name: [
                                            row[itemIndex] for row in table_values_array
                                        ]
                                    }
                                )
                            table_values_array = []
                        isLoop = True
                        category, item, value = None, None, None
                        # Stores items of a category listed in loop blocks
                        table_names = []
                        # Stores values of items in a loop as a single row
                        table_values = []
                    elif self.dataNameRE.match(line):
                        # Two step process STAR does not know contept of categories
                        m = self.dataNameRE.match(line)
                        flag = m.group("data_category")

                        tmp_category = self.dataCategoryItem.match(flag)
                        if tmp_category:
                            category = tmp_category.group("data_category")
                            item = tmp_category.group("category_item")
                        else:
                            category = ""
                            item = flag

                        remainder = m.group("remainder")
                        value = None
                        if isLoop and remainder != "":
                            """Append any data values following the last loop
                            category.item tag should any exist"""
                            table_values += self._tokenizeData(remainder)
                            line = ""
                        else:
                            line = remainder + "\n"
                        if not isLoop:
                            if line.strip() != "":
                                value = self._tokenizeData(line)
                            else:
                                # For cases where values are on the following
                                # line
                                try:
                                    line = next(f1)
                                    line_num += 1
                                except StopIteration:
                                    break
                            while value is None:
                                char_start = 1 if line.startswith(";") else 0
                                while line.startswith(
                                    ";"
                                ) and not line.rstrip().endswith("\n;"):
                                    try:
                                        line += next(f1)
                                        line_num += 1
                                    except StopIteration:
                                        break
                                value = (line[char_start : line.rfind("\n;")]).strip()
                                if char_start > 0:
                                    value = (
                                        line[char_start : line.rfind("\n;")]
                                    ).strip()
                                else:
                                    value = self._tokenizeData(" " + line)
                            if (ignoreCategories and category in ignoreCategories) or (
                                onlyCategories and category not in onlyCategories
                            ):
                                pass
                            else:
                                if category in data_block:
                                    data_block[category].update(
                                        {item: value if len(value) > 1 else value[0]}
                                    )
                                else:
                                    data_block.setdefault(
                                        category,
                                        _dict(
                                            {
                                                item: value
                                                if len(value) > 1
                                                else value[0]
                                            }
                                        ),
                                    )  # OrderedDict here preserves item order
                        else:
                            if (ignoreCategories and category in ignoreCategories) or (
                                onlyCategories and category not in onlyCategories
                            ):
                                skipCategory = True
                            else:
                                data_block.setdefault(
                                    category, _dict()
                                )  # OrderedDict here preserves item order
                                table_names.append(item)
                    else:
                        if multiLineValue is True:
                            table_values.append((line[1 : line.rfind("\n;")]).strip())
                            multiLineValue = False
                            line = line[line.rfind("\n;") + 2 :]
                            if line.strip() != "":
                                table_values += self._tokenizeData(line)
                        else:
                            table_values += self._tokenizeData(line)

                        if table_values != []:
                            table_values_array += table_values
                            table_values = []
                if isLoop is True and table_values_array != []:
                    isLoop = False
                    num_item = len(table_names)
                    for val_index, item in enumerate(table_names):
                        data_block[category][item] = table_values_array[
                            val_index::num_item
                        ]
                    table_values_array = []
                if data_block != {}:
                    mmcif_like_file[data_heading] = data_block
            return mmcif_like_file
        except KeyError as key_err:
            print("KeyError [line %i]: %s" % (line_num, str(key_err)))
        except IOError as io_err:
            print("IOException [line %i]: %s" % (line_num, str(io_err)))

#If useNumpy stays False, then it will look for the qconvex executable
useNumpy = False
useScipy = False
try:
    import numpy as np
    useNumpy = True
except:
    pass

if useNumpy:
    try:
        from scipy.spatial import ConvexHull
        useScipy = True
    except:
        pass

print('useNumpy: ', useNumpy)
print('useScipy: ', useScipy)

# Amino Acid mass, volumes,partial specific volumes and unified side chain sphere radius
# Mass is minus water (from SEDNTERP database)
# Volumes from Cohn and Edsall (1943) as corrected by Perkins_EJB_1986
# vbar is 0.60224(V/M)
#          Mass		Volume	   vbar		Radius of equiv. sphere of side chain volume
AA_data = {
 'ALA': (  71.0939,	87.2,      0.74,	1.852),
 'ARG': ( 156.203,	188.2,     0.73,	3.123),
 'ASN': ( 114.119,	120.1,     0.63,	2.422),
 'ASP': ( 115.104,	115.4,     0.60,	2.356),
 'CYS': ( 103.16,	106.7,     0.62,	2.224),
 'GLN': ( 128.16,	145.1,     0.68,	2.722),
 'GLU': ( 129.131,	140.9,     0.66,	2.676),
 'GLY': (  57.0669,	60.6,      0.64,	0.000),
 'HIS': ( 137.157,	152.4,     0.67,	2.798),
 'HSE': ( 137.157,	152.4,     0.67,	2.798), # CHARMM:neutral His, proton on NE2
 'HSP': ( 137.157,	152.4,     0.67,	2.798), # CHARMM:Protonated His
 'HSD': ( 137.157,	152.4,     0.67,	2.798), # CHARMM:neutral HIS, proton on ND1
 'ILE': ( 113.175,	168.9,     0.90,	2.957),
 'LEU': ( 113.175,	168.9,     0.90,	2.957),
 'LYS': ( 128.19,	174.3,     0.82,	3.005),
 'MET': ( 131.215,	163.1,     0.75,	2.903),
 'MSE': ( 178.125,	163.1,     0.75,	2.903),
 'PHE': ( 147.192,	187.9,     0.77,	3.121),
 'PRO': (  97.132,	122.4,     0.76,	2.453),
 'SER': (  87.093,	91.0,      0.63,	1.936),
 'THR': ( 101.12,	117.4,     0.70,	2.385),
 'TRP': ( 186.229,	228.5,     0.74,	3.422),
 'TYR': ( 163.192,	192.1,     0.71,	3.155),
 'VAL': (  99.148,	141.4,     0.86,	2.682)
}

# From voss_jmb_2005.pdf and nadassy_nar_2001.pdf
# Data for nucleoside, = nucleotide - 79.95 for PO3H
# Phosphate for na lacking terminal phosphates
# vbar is 0.60224(V/M)
#          Mass         Volume   vbar      ND
NA_data = {
 'A':  ( 267.25,       272.3,  0.614,   0.00), # Adenosine
 '1MA':  ( 281.25,       299.3,  0.641,   0.00), # 6-HYDRO-1-METHYLADENOSINE
 'MIA':  ( 383.45,       385.1,  0.605,   0.00), # 2-METHYLTHIO-N6-ISOPENTENYL-ADENOSINE
 'C':  ( 243.22,       248.1,  0.614,   0.00), # Cytidine
 '5MC':  ( 257.22,       275.1,  0.644,   0.00), # 5-METHYLCYTIDINE
 'OMC':  ( 257.22,       275.1,  0.644,   0.00), # O2'-METHYLYCYTIDINE
 'G':  ( 283.24,       279.0,  0.593,   0.00), # Guanosine
 '2MG':  ( 297.27,       306.0,  0.620,   0.00), # 2N-METHYLGUANOSINE
 '7MG':  ( 299.27,       306.0,  0.620,   0.00), # 7N-METHYL-8-HYDROGUANOSINE
 'M2G':  ( 311.27,       333.0,  0.644,   0.00), # N2-DIMETHYLGUANOSINE
 'OMG':  ( 297.27,       306.0,  0.620,   0.00), # O2'-METHYLGUANOSINE
 'YYG':  ( 508.54,       500.0,  0.593,   0.00), # MODIFIED GUANOSINE
 'YG':  ( 428.55,       427.1,  0.600,   0.00), # WYBUTOSINE
 'U':  ( 244.20,       243.9,  0.602,   0.00), # Uracil
 '5MU':  ( 258.20,       270.9,  0.632,   0.00), # 5-METHYLURIDINE
 '4SU':  ( 260.35,       270.5,  0.626,   0.00), # 4-THIOURIDINE
 'H2U':  ( 246.20,       243.9,  0.602,   0.00), # 5,6-DIHYDROURIDINE
 'PSU':  ( 244.20,       243.9,  0.602,   0.00), # PSEUDOURIDINE
 'I':  ( 268.23,       273.0,  0.613,   0.00), # Inosine
 'DA':  ( 251.25,       271.0,  0.650,   0.00), # 2'-DEOXYADENOSINE
 'DC':  ( 227.22,       246.8,  0.654,   0.00), # 2'-DEOXYCYTIDINE
 'DG':  ( 267.25,       277.7,  0.626,   0.00), # 2'-DEOXYGUANOSINE
 'DT':  ( 242.23,       264.4,  0.657,   0.00), # 2'-DEOXYTHYMIDINE
 'DI':  ( 252.22,       271.7,  0.649,   0.00), # 2'-DEOXYINOSINE
 'PO2':  (  62.97,	  52.4,  0.501,   0.00)
}

# For default saccharides use average vbar = 0.63
# Schuster, TM & Laue, TM, Eds, Modern Analytical UC, Springer, 2012
# or if found in, use
# Schuck, P. Zhao, H. Brautigam, C.A. and Ghirlando, R.
# Basic Principles of Analytical Ultracentrifugation
# CRC Press, 2016
#
# V = (vbar*M)/0.60224
#          Mass         Volume   vbar      ND
GL_data = {
 'NG6':  ( 301.27,       315.2,  0.630,   0.00), # N-ACETYL-D-GALACTOSAMINE 6-SULFATE
 'NAG':  ( 221.21,       231.4,  0.630,   0.00), # N-ACETYL-D-GLUCOSAMINE
 'BM3':  ( 221.21,       231.4,  0.630,   0.00), # 2-(ACETYLAMINO)-2-DEOXY-ALPHA-D-MANNOPYRANOSE
 'NGA':  ( 221.21,       231.4,  0.684,   0.00), # N-ACETYL-D-GALACTOSAMINE
 'GCU':  ( 194.14,       203.1,  0.630,   0.00), # D-GLUCURONIC ACID
 'IDR':  ( 194.14,       203.1,  0.630,   0.00), # L-IDURONIC ACID
 'BMA':  ( 180.16,       188.4,  0.607,   0.00), # BETA-D-MANNOSE
 'MAN':  ( 180.16,       188.4,  0.607,   0.00), # ALPHA-D-MANNOSE
 'GAL':  ( 180.16,       188.4,  0.622,   0.00), # BETA-D-GALACTOSE
 'GLA':  ( 180.16,       188.4,  0.622,   0.00), # ALPHA D-GALACTOSE
 'GLC':  ( 180.16,       188.4,  0.622,   0.00), # ALPHA-D-GLUCOSE
 'BGC':  ( 180.16,       188.4,  0.622,   0.00), # BETA-D-GLUCOSE
 'AOS':  ( 180.16,       188.4,  0.630,   0.00), # D-ALLOSE
 'GCS':  ( 179.17,       187.4,  0.630,   0.00), # Glucosamine
 'RAM':  ( 164.16,       171.7,  0.630,   0.00), # ALPHA-L-RHAMNOSE
 'FUC':  ( 164.16,       171.7,  0.671,   0.00),# ALPHA-L-FUCOSE
 'SIA':  ( 309.27,       299.9,  0.584,   0.00) # O-SIALIC ACID
}

# Detergents
# From Schuck, P. Zhao, H. Brautigam, C.A. and Ghirlando, R.
# Basic Principles of Analytical Ultracentrifugation
# CRC Press, 2016
#
# DOC from Durschlag, H.
# Specific Volumes of Biological Macromolecules and Some Other
#   Molecules of Biological Interest.
# in Thermodynamic Data for Biochemistry and Biotechnology, ed. H-J Hinz
# Springer-Verlag, 1986
#
# LMT from Suarez et al.
# JBC 259:13791, 1984
#
# FOS (DPC) from Kochendoerfer, et al
# Biochemistry 38:11905, 1999
#  But see Lauterwin, BBA 556:244, 1979 for vbar = 0.937
#
# V = (vbar*M)/0.60224
#          Mass         Volume   vbar      ND
DT_data = {
 'SB3':  ( 335.50,       533.9,  0.957,   0.00), # n-Dodecyl-N,N-dimethyl-3-ammonio-1-propanesulfonate
 'LMT':  ( 510.62,       691.0,  0.820,   0.00), # DODECYL-BETA-D-MALTOSIDE
 'BOG':  ( 417.02,       594.8,  0.859,   0.00), # B-OCTYLGLUCOSIDE
 'LDA':  ( 229.40,       429.7,  1.128,   0.00), # LAURYL DIMETHYLAMINE-N-OXIDE
 'SDS':  ( 266.40,       384.8,  0.880,   0.00), # DODECYL SULFATE
 'DXC':  ( 392.47,       507.0,  0.778,   0.00), # DEOXYCHOLATE
 'FOS':  ( 351.46,       548.6,  0.940,   0.00)  # DODECYLPHOSPHOCHOLINE
}

def openGzip(file_path, mode="rt"):
    try:
        return (
            gzip.open(file_path, mode)
            if mimetypes.guess_type(file_path)[1] == "gzip"
            else open(file_path, mode)
        )
    except Exception as gzip_io_error:
        print("[Error opening mmCIF file]: %s" % gzip_io_error.message)
        return None

def distance(ax,ay,az, bx,by,bz):
    #Euclidean distance between two atoms
    return math.sqrt((ax - bx)**2.0 + (ay - by)**2.0 + (az - bz)**2.0)

def get_coords(model_array):
    # Get x,y,z coords of CA and CB atoms
    coords = []
    for row in range(len(model_array)):
        coords.append((float(model_array[row][5]), float(model_array[row][6]), float(model_array[row][7])))
    if useNumpy:
        return np.array(coords)
    else:
        return coords

def model_from_pdb(file):
    # Makes a reduced atom  model of the pdb file
    # This is the model used to make the convex hull
    # For proteins the CB is displaced along the CA-CB vector a distance
    #   equal to the radius of a sphere equal to the volume of the side chain
    # For nucleic acids only the phosphate, oxygen and nitrogen atoms are used
    # For saccharides only the phosphate, oxygen and nitrogen atoms are used
    # For detergents only the oxygen and nitrogen atoms are used

    # Read the PDB file
    data = open(file, 'r').readlines()

    # Initialize all relevant atom record array
    # Used for Rg
    all_atm_rec = []
    # Initial protein atoms from PDB file
    prot_rec = []
    # Initial nucleic acid atoms from PDB file
    na_rec = []
    # Initial saccharide atoms from PDB file
    gl_rec = []
    # Initial detergent atoms from PDB file
    dt_rec = []
    # Reduced list of atoms
    atom_rec = []
    struc_name = ' '
    num_MG = 0
    num_MN = 0

    # Get all relevant atoms even if in wrong order
    for line in data:
        # Get name of structure
        if (line[:6] == 'HEADER'):
            struc_name = line.strip()

        elif (line[:4] == 'ATOM'):
            resname_pad = line[17:20]
            resname_q = resname_pad.strip()
            resname = resname_q.replace("'","")

            if resname != 'TIP' and resname != 'HOH':
                all_atm_rec.append(line)
                if resname in AA_data and (line[13:15] == 'N ' or \
                  line[13:15] == 'CA' or \
                  line[13:15] == 'C ' or \
                  line[13:15] == 'O ' or \
                  line[13:16] == 'OT1' or \
                  line[13:15] == 'CB') and \
                  (line[16] != 'B') and \
                  (line[:3] != 'END'):
                    prot_rec.append(line)
                # Nucleic acids
                elif resname in NA_data and (line[13:14] == 'N' or \
                  line[13:14] == 'O' or line[13:14] == 'P'):
                    na_rec.append(line)
                # Saccharides
                elif resname in GL_data and (line[13:14] == 'N' or \
                  line[13:14] == 'O' or line[13:14] == 'P'):
                    gl_rec.append(line)
                # Detergents
                elif resname in DT_data and (line[13:14] == 'N' or \
                  line[13:14] == 'O' or line[14:15] == 'O'):
                    dt_rec.append(line)

        elif (line[:6] == 'HETATM'):
            resname_pad = line[17:20]
            resname_q = resname_pad.strip()
            resname = resname_q.replace("'","")
            if resname != 'TIP' and resname != 'HOH':
                for rec in  NA_data.items():
                    if resname == rec[0]:
                        all_atm_rec.append(line)
                for rec in  GL_data.items():
                    if resname == rec[0]:
                        all_atm_rec.append(line)
                for rec in  DT_data.items():
                    if resname == rec[0]:
                        all_atm_rec.append(line)
                if line[12:14] == 'MG':
                    num_MG += 1
                elif line[12:14] == 'MN':
                    num_MN += 1
                # Some nucleic acids are in HETATM (!)
                elif resname in NA_data and (line[13:14] == 'N' or \
                  line[13:14] == 'O' or line[13:14] == 'P'):
                    na_rec.append(line)
                # Some glycosylations are in HETATM (!)
                elif resname in GL_data and (line[13:14] == 'N' or \
                  line[13:14] == 'O' or line[13:14] == 'P'):
                    gl_rec.append(line)
                # Some detergents are in HETATM (!)
                elif resname in DT_data and (line[13:14] == 'N' or \
                  line[13:14] == 'O'):
                    dt_rec.append(line)

    # Get only the N atoms from the protein
    n_rec = []
    for line in prot_rec:
        if (line[13:14] == 'N') and \
           (line[0:3] != 'END'):
            n_rec.append((line))

    # Get only the CA atoms from the protein
    ca_rec = []
    for line in prot_rec:
        if (line[13:15] == 'CA') and \
           (line[0:3] != 'END'):
            ca_rec.append((line))

    # Get only the O atoms from the protein
    o_rec = []
    for line in prot_rec:
        if (line[13:14] == 'O') and \
           (line[0:3] != 'END'):
            o_rec.append((line))

    # Get only the CB atoms from the protein
    cb_rec = []
    for line in prot_rec:
        if (line[13:15] == 'CB') and \
           (line[0:3] != 'END'):
            cb_rec.append((line))

    backbone_atms_check = (2*(len(n_rec))) - len(ca_rec) - len(o_rec)
    if backbone_atms_check != 0:
        print(' Backbone atom missing.')
        print(' Delete residue missing backbone atom(s) and try again.')
        sys.exit()

    # Count the CB atoms and add to atom list in correct order
    num_cb = len(cb_rec)
    line_counter = 0
    cb_counter = 0
    for line in prot_rec:
        if (line[13:15] == 'N ' or \
          line[13:15] == 'CA' or \
          line[13:15] == 'C ' or \
          line[13:15] == 'O ' or \
          line[13:16] == 'OT1'):
            atom_rec.append((line))
            line_counter += 1
            if (math.fmod(line_counter,4) == 0 and line[17:20] != 'GLY' \
                 and cb_counter < num_cb):
                atom_rec.append((cb_rec[cb_counter]))
                cb_counter +=1
            elif (math.fmod(line_counter,4) == 0 and line[17:20] != 'GLY' \
                 and cb_counter >= num_cb):
                print(' Found non-gly residue with no CB atom')
                print(' Change the name of the offending residue to GLY and try again.')
                sys.exit()

    # Put atom info into model array
    num_atoms =  len(atom_rec)
    model_array = [['X' for j in range(8)] for i in range(num_atoms)]
    for row in range(len(model_array)):
        model_array[row][0] = row
        model_array[row][1] = (atom_rec[row][11:16]).strip()
        model_array[row][2] = (atom_rec[row][17:20]).strip()
        if atom_rec[row][20:22] == '  ':
            model_array[row][3] = 'X'
        else:
            model_array[row][3] = (atom_rec[row][20:22]).strip()
        model_array[row][4] = (atom_rec[row][22:26]).strip()
        model_array[row][5] = (atom_rec[row][30:38]).strip()
        model_array[row][6] = (atom_rec[row][38:46]).strip()
        model_array[row][7] = (atom_rec[row][46:54]).strip()

    # Make unified side chain as CB
    for row in range(len(model_array)):
        if (model_array[row][1] == 'CA' and model_array[row][2] != 'GLY'):
            rnam = model_array[row][2]
            # CA coords
            ca_x = float(model_array[row][5])
            ca_y = float(model_array[row][6])
            ca_z = float(model_array[row][7])
            # CB coords (third next coord: C,O,CB)
            cb_x = float(model_array[row+3][5])
            cb_y = float(model_array[row+3][6])
            cb_z = float(model_array[row+3][7])
            # CA-CB distance
            dx = cb_x - ca_x
            dy = cb_y - ca_y
            dz = cb_z - ca_z
            ca_cb_dist = distance(ca_x,ca_y,ca_z,cb_x,cb_y,cb_z)
            # Extend CB outward along CA-CB vector
            aa_mass, vol_aa, vbar_aa, sc_rad_aa = AA_data[rnam]
            new_cb_x = ca_x + (sc_rad_aa * dx/ca_cb_dist)
            new_cb_y = ca_y + (sc_rad_aa * dy/ca_cb_dist)
            new_cb_z = ca_z + (sc_rad_aa * dz/ca_cb_dist)
            model_array[row+3][5] = str(new_cb_x)
            model_array[row+3][6] = str(new_cb_y)
            model_array[row+3][7] = str(new_cb_z)

    # Put nucleic acid atoms into initial array
    num_atoms =  len(na_rec)
    na_array = [['X' for j in range(8)] for i in range(num_atoms)]
    for row in range(len(na_array)):
        na_array[row][0] = row
        na_array[row][1] = (na_rec[row][11:16]).strip()
        na_array[row][2] = (na_rec[row][17:20]).strip()
        na_array[row][3] = (na_rec[row][20:22]).strip()
        na_array[row][4] = (na_rec[row][22:26]).strip()
        na_array[row][5] = (na_rec[row][30:38]).strip()
        na_array[row][6] = (na_rec[row][38:46]).strip()
        na_array[row][7] = (na_rec[row][46:54]).strip()

    # Add nucleic acid atoms to model array
    for row in range(len(na_array)):
        model_array.append(na_array[row])

    # Put saccharide atoms into initial array
    num_atoms =  len(gl_rec)
    gl_array = [['X' for j in range(8)] for i in range(num_atoms)]
    for row in range(len(gl_array)):
        gl_array[row][0] = row
        gl_array[row][1] = (gl_rec[row][11:16]).strip()
        gl_array[row][2] = (gl_rec[row][17:20]).strip()
        gl_array[row][3] = (gl_rec[row][20:22]).strip()
        gl_array[row][4] = (gl_rec[row][22:26]).strip()
        gl_array[row][5] = (gl_rec[row][30:38]).strip()
        gl_array[row][6] = (gl_rec[row][38:46]).strip()
        gl_array[row][7] = (gl_rec[row][46:54]).strip()

    # Add saccharide atoms into model array
    for row in range(len(gl_array)):
        model_array.append(gl_array[row])

    # Put detergent atoms into initial array
    num_atoms =  len(dt_rec)
    dt_array = [['X' for j in range(8)] for i in range(num_atoms)]
    for row in range(len(dt_array)):
        dt_array[row][0] = row
        dt_array[row][1] = (dt_rec[row][11:16]).strip()
        dt_array[row][2] = (dt_rec[row][17:20]).strip()
        dt_array[row][3] = (dt_rec[row][20:22]).strip()
        dt_array[row][4] = (dt_rec[row][22:26]).strip()
        dt_array[row][5] = (dt_rec[row][30:38]).strip()
        dt_array[row][6] = (dt_rec[row][38:46]).strip()
        dt_array[row][7] = (dt_rec[row][46:54]).strip()

    # Add detergent atoms into model array
    for row in range(len(dt_array)):
        model_array.append(dt_array[row])

    # Convert all_atm_rec to multi-item list
    all_atm_array = [['X' for j in range(8)] for i in range(len(all_atm_rec))]
    for row in range(len(all_atm_rec)):
        all_atm_array[row][0] = row					#Atom Index
        all_atm_array[row][1] = (all_atm_rec[row][11:16]).strip()	#Atom Name
        all_atm_array[row][2] = (all_atm_rec[row][17:20]).strip()	#Residue Name
        all_atm_array[row][3] = (all_atm_rec[row][20:22]).strip()	#ChainID
        all_atm_array[row][4] = (all_atm_rec[row][22:26]).strip()	#Residue Number
        all_atm_array[row][5] = (all_atm_rec[row][30:38]).strip()	#x
        all_atm_array[row][6] = (all_atm_rec[row][38:46]).strip()	#y
        all_atm_array[row][7] = (all_atm_rec[row][46:54]).strip()	#z

    return all_atm_array,num_MG,num_MN,model_array


def model_from_chain(file, chainID='all'):
    # Makes a reduced atom  model of a specific chain in the pdb file
    # This is the model used to make the convex hull
    # For proteins the CB is displaced along the CA-CB vector a distance
    #   equal to the radius of a sphere equal to the volume of the side chain
    # For nucleic acids only the phosphate, oxygen and nitrogen atoms are used
    # For saccharides only the phosphate, oxygen and nitrogen atoms are used
    # For detergents only the oxygen and nitrogen atoms are used

    # Read the PDB file
    data = open(file, 'r').readlines()

    # Initialize all relevant atom record array
    # Used for Rg
    all_atm_rec = []
    # Initial protein atoms from PDB file
    prot_rec = []
    # Initial nucleic acid atoms from PDB file
    na_rec = []
    # Initial saccharide atoms from PDB file
    gl_rec = []
    # Initial detergent atoms from PDB file
    dt_rec = []
    # Reduced list of atoms
    atom_rec = []
    struc_name = ' '
    num_MG = 0
    num_MN = 0

    # Get all relevant atoms even if in wrong order
    for line in data:
        # Get name of structure
        if (line[:6] == 'HEADER'):
            struc_name = line.strip()

        # Remove the chains that do not belong
        elif ((chainID != 'all') and (line[21:22] != chainID):
            struc_name = line.strip()

        elif (line[:4] == 'ATOM'):
            resname_pad = line[17:20]
            resname_q = resname_pad.strip()
            resname = resname_q.replace("'","")

            if resname != 'TIP' and resname != 'HOH':
                all_atm_rec.append(line)
                if resname in AA_data and (line[13:15] == 'N ' or \
                  line[13:15] == 'CA' or \
                  line[13:15] == 'C ' or \
                  line[13:15] == 'O ' or \
                  line[13:16] == 'OT1' or \
                  line[13:15] == 'CB') and \
                  (line[16] != 'B') and \
                  (line[:3] != 'END'):
                    prot_rec.append(line)
                # Nucleic acids
                elif resname in NA_data and (line[13:14] == 'N' or \
                  line[13:14] == 'O' or line[13:14] == 'P'):
                    na_rec.append(line)
                # Saccharides
                elif resname in GL_data and (line[13:14] == 'N' or \
                  line[13:14] == 'O' or line[13:14] == 'P'):
                    gl_rec.append(line)
                # Detergents
                elif resname in DT_data and (line[13:14] == 'N' or \
                  line[13:14] == 'O' or line[14:15] == 'O'):
                    dt_rec.append(line)

        elif (line[:6] == 'HETATM'):
            resname_pad = line[17:20]
            resname_q = resname_pad.strip()
            resname = resname_q.replace("'","")
            if resname != 'TIP' and resname != 'HOH':
                for rec in  NA_data.items():
                    if resname == rec[0]:
                        all_atm_rec.append(line)
                for rec in  GL_data.items():
                    if resname == rec[0]:
                        all_atm_rec.append(line)
                for rec in  DT_data.items():
                    if resname == rec[0]:
                        all_atm_rec.append(line)
                if line[12:14] == 'MG':
                    num_MG += 1
                elif line[12:14] == 'MN':
                    num_MN += 1
                # Some nucleic acids are in HETATM (!)
                elif resname in NA_data and (line[13:14] == 'N' or \
                  line[13:14] == 'O' or line[13:14] == 'P'):
                    na_rec.append(line)
                # Some glycosylations are in HETATM (!)
                elif resname in GL_data and (line[13:14] == 'N' or \
                  line[13:14] == 'O' or line[13:14] == 'P'):
                    gl_rec.append(line)
                # Some detergents are in HETATM (!)
                elif resname in DT_data and (line[13:14] == 'N' or \
                  line[13:14] == 'O'):
                    dt_rec.append(line)

    # Get only the N atoms from the protein
    n_rec = []
    for line in prot_rec:
        if (line[13:14] == 'N') and \
           (line[0:3] != 'END'):
            n_rec.append((line))

    # Get only the CA atoms from the protein
    ca_rec = []
    for line in prot_rec:
        if (line[13:15] == 'CA') and \
           (line[0:3] != 'END'):
            ca_rec.append((line))

    # Get only the O atoms from the protein
    o_rec = []
    for line in prot_rec:
        if (line[13:14] == 'O') and \
           (line[0:3] != 'END'):
            o_rec.append((line))

    # Get only the CB atoms from the protein
    cb_rec = []
    for line in prot_rec:
        if (line[13:15] == 'CB') and \
           (line[0:3] != 'END'):
            cb_rec.append((line))

    backbone_atms_check = (2*(len(n_rec))) - len(ca_rec) - len(o_rec)
    if backbone_atms_check != 0:
        print(' Backbone atom missing.')
        print(' Delete residue missing backbone atom(s) and try again.')
        sys.exit()

    # Count the CB atoms and add to atom list in correct order
    num_cb = len(cb_rec)
    line_counter = 0
    cb_counter = 0
    for line in prot_rec:
        if (line[13:15] == 'N ' or \
          line[13:15] == 'CA' or \
          line[13:15] == 'C ' or \
          line[13:15] == 'O ' or \
          line[13:16] == 'OT1'):
            atom_rec.append((line))
            line_counter += 1
            if (math.fmod(line_counter,4) == 0 and line[17:20] != 'GLY' \
                 and cb_counter < num_cb):
                atom_rec.append((cb_rec[cb_counter]))
                cb_counter +=1
            elif (math.fmod(line_counter,4) == 0 and line[17:20] != 'GLY' \
                 and cb_counter >= num_cb):
                print(' Found non-gly residue with no CB atom')
                print(' Change the name of the offending residue to GLY and try again.')
                sys.exit()

    # Put atom info into model array
    num_atoms =  len(atom_rec)
    model_array = [['X' for j in range(8)] for i in range(num_atoms)]
    for row in range(len(model_array)):
        model_array[row][0] = row
        model_array[row][1] = (atom_rec[row][11:16]).strip()
        model_array[row][2] = (atom_rec[row][17:20]).strip()
        if atom_rec[row][20:22] == '  ':
            model_array[row][3] = 'X'
        else:
            model_array[row][3] = (atom_rec[row][20:22]).strip()
        model_array[row][4] = (atom_rec[row][22:26]).strip()
        model_array[row][5] = (atom_rec[row][30:38]).strip()
        model_array[row][6] = (atom_rec[row][38:46]).strip()
        model_array[row][7] = (atom_rec[row][46:54]).strip()

    # Make unified side chain as CB
    for row in range(len(model_array)):
        if (model_array[row][1] == 'CA' and model_array[row][2] != 'GLY'):
            rnam = model_array[row][2]
            # CA coords
            ca_x = float(model_array[row][5])
            ca_y = float(model_array[row][6])
            ca_z = float(model_array[row][7])
            # CB coords (third next coord: C,O,CB)
            cb_x = float(model_array[row+3][5])
            cb_y = float(model_array[row+3][6])
            cb_z = float(model_array[row+3][7])
            # CA-CB distance
            dx = cb_x - ca_x
            dy = cb_y - ca_y
            dz = cb_z - ca_z
            ca_cb_dist = distance(ca_x,ca_y,ca_z,cb_x,cb_y,cb_z)
            # Extend CB outward along CA-CB vector
            aa_mass, vol_aa, vbar_aa, sc_rad_aa = AA_data[rnam]
            new_cb_x = ca_x + (sc_rad_aa * dx/ca_cb_dist)
            new_cb_y = ca_y + (sc_rad_aa * dy/ca_cb_dist)
            new_cb_z = ca_z + (sc_rad_aa * dz/ca_cb_dist)
            model_array[row+3][5] = str(new_cb_x)
            model_array[row+3][6] = str(new_cb_y)
            model_array[row+3][7] = str(new_cb_z)

    # Put nucleic acid atoms into initial array
    num_atoms =  len(na_rec)
    na_array = [['X' for j in range(8)] for i in range(num_atoms)]
    for row in range(len(na_array)):
        na_array[row][0] = row
        na_array[row][1] = (na_rec[row][11:16]).strip()
        na_array[row][2] = (na_rec[row][17:20]).strip()
        na_array[row][3] = (na_rec[row][20:22]).strip()
        na_array[row][4] = (na_rec[row][22:26]).strip()
        na_array[row][5] = (na_rec[row][30:38]).strip()
        na_array[row][6] = (na_rec[row][38:46]).strip()
        na_array[row][7] = (na_rec[row][46:54]).strip()

    # Add nucleic acid atoms to model array
    for row in range(len(na_array)):
        model_array.append(na_array[row])

    # Put saccharide atoms into initial array
    num_atoms =  len(gl_rec)
    gl_array = [['X' for j in range(8)] for i in range(num_atoms)]
    for row in range(len(gl_array)):
        gl_array[row][0] = row
        gl_array[row][1] = (gl_rec[row][11:16]).strip()
        gl_array[row][2] = (gl_rec[row][17:20]).strip()
        gl_array[row][3] = (gl_rec[row][20:22]).strip()
        gl_array[row][4] = (gl_rec[row][22:26]).strip()
        gl_array[row][5] = (gl_rec[row][30:38]).strip()
        gl_array[row][6] = (gl_rec[row][38:46]).strip()
        gl_array[row][7] = (gl_rec[row][46:54]).strip()

    # Add saccharide atoms into model array
    for row in range(len(gl_array)):
        model_array.append(gl_array[row])

    # Put detergent atoms into initial array
    num_atoms =  len(dt_rec)
    dt_array = [['X' for j in range(8)] for i in range(num_atoms)]
    for row in range(len(dt_array)):
        dt_array[row][0] = row
        dt_array[row][1] = (dt_rec[row][11:16]).strip()
        dt_array[row][2] = (dt_rec[row][17:20]).strip()
        dt_array[row][3] = (dt_rec[row][20:22]).strip()
        dt_array[row][4] = (dt_rec[row][22:26]).strip()
        dt_array[row][5] = (dt_rec[row][30:38]).strip()
        dt_array[row][6] = (dt_rec[row][38:46]).strip()
        dt_array[row][7] = (dt_rec[row][46:54]).strip()

    # Add detergent atoms into model array
    for row in range(len(dt_array)):
        model_array.append(dt_array[row])

    # Convert all_atm_rec to multi-item list
    all_atm_array = [['X' for j in range(8)] for i in range(len(all_atm_rec))]
    for row in range(len(all_atm_rec)):
        all_atm_array[row][0] = row					#Atom Index
        all_atm_array[row][1] = (all_atm_rec[row][11:16]).strip()	#Atom Name
        all_atm_array[row][2] = (all_atm_rec[row][17:20]).strip()	#Residue Name
        all_atm_array[row][3] = (all_atm_rec[row][20:22]).strip()	#ChainID
        all_atm_array[row][4] = (all_atm_rec[row][22:26]).strip()	#Residue Number
        all_atm_array[row][5] = (all_atm_rec[row][30:38]).strip()	#x
        all_atm_array[row][6] = (all_atm_rec[row][38:46]).strip()	#y
        all_atm_array[row][7] = (all_atm_rec[row][46:54]).strip()	#z

    return all_atm_array,num_MG,num_MN,model_array

def model_from_cif(file):
    # Initialize all relevant atom record array
    # Used for Rg
    all_atm_rec = []
    # Initial protein atoms from PDB file
    prot_rec = []
    # Initial nucleic acid atoms from PDB file
    na_rec = []
    # Initial saccharide atoms from PDB file
    gl_rec = []
    # Initial detergent atoms from PDB file
    dt_rec = []
    # Reduced list of atoms
    atom_rec = []
    struc_name = ' '
    num_MG = 0
    num_MN = 0

    cif_dict = MMCIF2Dict().parse(file)
    id = list(cif_dict.keys())[0]
    polymer_name = cif_dict[id]['_entity']['pdbx_description'][0]
    records = list(cif_dict[id]['_atom_site'])
    record_types = list(cif_dict[id]['_atom_site']['group_PDB'])
    atms = list(cif_dict[id]['_atom_site']['label_atom_id'])
    alt_rotamer = list(cif_dict[id]['_atom_site']['label_alt_id'])
    res_nams = list(cif_dict[id]['_atom_site']['auth_comp_id'])
    chain_id = list(cif_dict[id]['_atom_site']['auth_asym_id'])
    res_id = list(cif_dict[id]['_atom_site']['auth_seq_id'])
    crdxs = list(cif_dict[id]['_atom_site']['Cartn_x'])
    crdys = list(cif_dict[id]['_atom_site']['Cartn_y'])
    crdzs = list(cif_dict[id]['_atom_site']['Cartn_z'])
    for i in range(len(record_types)):
        if record_types[i] == 'ATOM':
            if res_nams[i] != 'TIP' and res_nams[i] != 'HOH':
                all_atm_rec.append([i,atms[i],res_nams[i],chain_id[i],res_id[i],crdxs[i],crdys[i],crdzs[i]])
            if res_nams[i] in AA_data and (atms[i] == 'N' or \
              atms[i] == 'CA' or \
              atms[i] == 'C' or \
              atms[i] == 'O' or \
              atms[i] == 'OT1' or \
              atms[i] == 'CB' and \
              alt_rotamer[i] != 'B'):
                prot_rec.append([atms[i],res_nams[i],chain_id[i],res_id[i],crdxs[i],crdys[i],crdzs[i]])

            #Nucleic Acids
            elif res_nams[i].replace("'","") in NA_data and (atms[i][:1] == 'N' or \
              atms[i][:1] == 'O' or \
              atms[i][:1] == 'P'):
                na_rec.append([atms[i],res_nams[i],chain_id[i],res_id[i],crdxs[i],crdys[i],crdzs[i]])

            # Saccharides
            elif res_nams[i] in GL_data and (atms[i][:1] == 'N' or \
              atms[i][:1] == 'O' or \
              atms[i][:1] == 'P'):
                na_rec.append([atms[i],res_nams[i],chain_id[i],res_id[i],crdxs[i],crdys[i],crdzs[i]])

            # Detergents
            elif res_nams[i] in DT_data and (atms[i][:1] == 'N' or \
              atms[i][:1] == 'O' or \
              atms[i][1:2] == 'O'):
                dt_rec.append([atms[i],res_nams[i],chain_id[i],res_id[i],crdxs[i],crdys[i],crdzs[i]])
        elif record_types[i] == 'HETATM':
            if res_nams[i] != 'TIP' and res_nams[i] != 'HOH':
                for rec in NA_data.items():
                    if res_nams[i] == rec[0]:
                        all_atm_rec.append([i,atms[i],res_nams[i],chain_id[i],res_id[i],crdxs[i],crdys[i],crdzs[i]])
                for rec in GL_data.items():
                    if res_nams[i] == rec[0]:
                        all_atm_rec.append([i,atms[i],res_nams[i],chain_id[i],res_id[i],crdxs[i],crdys[i],crdzs[i]])
                for rec in DT_data.items():
                    if res_nams[i] == rec[0]:
                        all_atm_rec.append([i,atms[i],res_nams[i],chain_id[i],res_id[i],crdxs[i],crdys[i],crdzs[i]])
                if atms[i] =='MG':
                    num_MG += 1
                elif atms[i] =='MN':
                    num_MG += 1
                # Some nucleic acids are in HETATM (!)
                elif res_nams[i] in NA_data and (atms[i][:1] == 'N' or \
                  atms[i][:1] == 'O' or \
                  atms[i][:1] == 'P'):
                    na_rec.append([atms[i],res_nams[i],chain_id[i],res_id[i],crdxs[i],crdys[i],crdzs[i]])
                # Some glycosylations are in HETATM (!)
                elif res_nams[i] in GL_data and (atms[i][:1] == 'N' or \
                  atms[i][:1] == 'O' or \
                  atms[i][:1] == 'P'):
                    gl_rec.append([atms[i],res_nams[i],chain_id[i],res_id[i],crdxs[i],crdys[i],crdzs[i]])
                # Some detergents are in HETATM (!)
                elif res_nams[i] in DT_data and (atms[i][:1] == 'N' or \
                  atms[i][:1] == 'O' or \
                  atms[i][1:2] == 'O'):
                    dt_rec.append([atms[i],res_nams[i],chain_id[i],res_id[i],crdxs[i],crdys[i],crdzs[i]])
    # Get only the N atoms from the protein
    n_rec = []
    for line in prot_rec:
        if (line[13:14] == 'N') and \
           (line[0:3] != 'END'):
            n_rec.append((line))

    # Get only the CA atoms from the protein
    ca_rec = []
    for line in prot_rec:
        if (line[13:15] == 'CA') and \
           (line[0:3] != 'END'):
            ca_rec.append((line))

    # Get only the O atoms from the protein
    o_rec = []
    for line in prot_rec:
        if (line[13:14] == 'O') and \
           (line[0:3] != 'END'):
            o_rec.append((line))

    # Get only the CB atoms from the protein
    cb_rec = []
    for line in prot_rec:
        if (line[0] == 'CB'):
            cb_rec.append((line))

    backbone_atms_check = (2*(len(n_rec))) - len(ca_rec) - len(o_rec)
    if backbone_atms_check != 0:
        print(' Backbone atom missing.')
        print(' Delete residue missing backbone atom(s) and try again.')
        sys.exit()

    # Count the CB atoms and add to atom list in correct order
    num_cb = len(cb_rec)
    line_counter = 0
    cb_counter = 0
    for line in prot_rec:
        if (line[0] == 'N' or \
          line[0] == 'CA' or \
          line[0] == 'C' or \
          line[0] == 'O' or \
          line[0] == 'OT1'):
            atom_rec.append((line))
            line_counter += 1
            if (math.fmod(line_counter,4) == 0 and line[1] != 'GLY' \
                 and cb_counter < num_cb):
                atom_rec.append((cb_rec[cb_counter]))
                cb_counter +=1
            elif (math.fmod(line_counter,4) == 0 and line[1] != 'GLY' \
                 and cb_counter >= num_cb):
                print(' Found non-gly residue with no CB atom')
                print(' Change the name of the offending residue to GLY and try again.')
                sys.exit()

    # Put atom info into model array
    num_atoms =  len(atom_rec)
    model_array = [['X' for j in range(8)] for i in range(num_atoms)]
    for row in range(len(model_array)):
        model_array[row][0] = row                 #Atom Index
        model_array[row][1] = (atom_rec[row][0])  #Atom Name
        model_array[row][2] = (atom_rec[row][1])  #Residue Name
        if atom_rec[row][2] == '?':
            model_array[row][3] = 'X'
        else:
            model_array[row][3] = (atom_rec[row][2]) #ChainID
        model_array[row][4] = (atom_rec[row][3])  #Residue Number
        model_array[row][5] = (atom_rec[row][4])  #x
        model_array[row][6] = (atom_rec[row][5])  #y
        model_array[row][7] = (atom_rec[row][6])  #z

    # Make unified side chain as CB
    for row in range(len(model_array)):
        if (model_array[row][1] == 'CA' and model_array[row][2] != 'GLY'):
            rnam = model_array[row][2]
            # CA coords
            ca_x = float(model_array[row][5])
            ca_y = float(model_array[row][6])
            ca_z = float(model_array[row][7])
            # CB coords (third next coord: C,O,CB)
            cb_x = float(model_array[row+3][5])
            cb_y = float(model_array[row+3][6])
            cb_z = float(model_array[row+3][7])
            # CA-CB distance
            dx = cb_x - ca_x
            dy = cb_y - ca_y
            dz = cb_z - ca_z
            ca_cb_dist = distance(ca_x,ca_y,ca_z,cb_x,cb_y,cb_z)
            # Extend CB outward along CA-CB vector
            aa_mass, vol_aa, vbar_aa, sc_rad_aa = AA_data[rnam]
            new_cb_x = ca_x + (sc_rad_aa * dx/ca_cb_dist)
            new_cb_y = ca_y + (sc_rad_aa * dy/ca_cb_dist)
            new_cb_z = ca_z + (sc_rad_aa * dz/ca_cb_dist)
            model_array[row+3][5] = str(new_cb_x)
            model_array[row+3][6] = str(new_cb_y)
            model_array[row+3][7] = str(new_cb_z)

    # Put nucleic acid atoms into initial array
    num_atoms =  len(na_rec)
    na_array = [['X' for j in range(8)] for i in range(num_atoms)]
    for row in range(len(na_array)):
        na_array[row][0] = row
        na_array[row][1] = (na_rec[row][0]) #Atom Name
        na_array[row][2] = (na_rec[row][1]) #Residue Name
        na_array[row][3] = (na_rec[row][2]) #ChainID
        na_array[row][4] = (na_rec[row][3]) #Residue Number
        na_array[row][5] = (na_rec[row][4]) #x
        na_array[row][6] = (na_rec[row][5]) #y
        na_array[row][7] = (na_rec[row][6]) #z

    # Add nucleic acid atoms to model array
    for row in range(len(na_array)):
        model_array.append(na_array[row])

    # Put saccharide atoms into initial array
    num_atoms =  len(gl_rec)
    gl_array = [['X' for j in range(8)] for i in range(num_atoms)]
    for row in range(len(gl_array)):
        gl_array[row][0] = row
        gl_array[row][1] = (gl_rec[row][0])
        gl_array[row][2] = (gl_rec[row][1])
        gl_array[row][3] = (gl_rec[row][2])
        gl_array[row][4] = (gl_rec[row][3])
        gl_array[row][5] = (gl_rec[row][4])
        gl_array[row][6] = (gl_rec[row][5])
        gl_array[row][7] = (gl_rec[row][6])

    # Add saccharide atoms into model array
    for row in range(len(gl_array)):
        model_array.append(gl_array[row])

    # Put detergent atoms into initial array
    num_atoms =  len(dt_rec)
    dt_array = [['X' for j in range(8)] for i in range(num_atoms)]
    for row in range(len(dt_array)):
        dt_array[row][0] = row
        dt_array[row][1] = (dt_rec[row][0])
        dt_array[row][2] = (dt_rec[row][1])
        dt_array[row][3] = (dt_rec[row][2])
        dt_array[row][4] = (dt_rec[row][3])
        dt_array[row][5] = (dt_rec[row][4])
        dt_array[row][6] = (dt_rec[row][5])
        dt_array[row][7] = (dt_rec[row][6])

    # Add detergent atoms into model array
    for row in range(len(dt_array)):
        model_array.append(dt_array[row])

    return all_atm_rec,num_MG,num_MN,model_array

def write_pdb(model_array, filename):
    # Write out reduced atom model in PDB format for display
    if type(filename) is type(''):
        ff = open(filename, 'w')
    else:
        ff = filename

    write = ff.write

    pdbfmt = 'ATOM  %5d  %-3s %3s%2s%4d    %8.3f%8.3f%8.3f  0.00  0.00\n'
    for i in range(len(model_array)):
        a = model_array[i]
        atm_name = a[1].strip()
        chainID = a[3].strip()
        write(pdbfmt % (int(a[0]+1),atm_name,a[2],chainID,int(a[4]),float(a[5]),float(a[6]),float(a[7])))
    write('END\n')

    ff.close()

def Sved(all_atm_rec,num_MG,num_MN,model_array):
    #
    # Main function: Does most things and calls HullVolume
    #
    # Use Cohn-Edsall eq. to calc protein specific volume
    # vbar_prot = sum(ni*Mi*vbar_aa) / (ni*Mi) where n is number, M is mass
    #
    prot_mol_mass = 0.0
    numerator = 0.0
    AA = 0
    NA = 0
    GL = 0
    DT = 0

    # Radius of Gyration
    # Calc center of mass
    X = 0.0
    Y = 0.0
    Z = 0.0
    tot = 0.0
    for row in range(len(all_atm_rec)):
        X = X + (float(all_atm_rec[row][5]))
        Y = Y + (float(all_atm_rec[row][6]))
        Z = Z + (float(all_atm_rec[row][7]))
        tot += 1

    com_x = (X/tot)
    com_y = (Y/tot)
    com_z = (Z/tot)
    Rg2  = 0.0
    for row in range(len(all_atm_rec)):
        Rg2 += ((distance(com_x, com_y, com_z, float(all_atm_rec[row][5]),\
               float(all_atm_rec[row][6]), float(all_atm_rec[row][7])))**2)
    Rg = math.sqrt(Rg2/tot)
    # Asphericity
    # 0 for sphere, 1 for rod
    asphr = 0.0
    if useNumpy:
        Ixx,Ixy,Ixz,Iyy,Iyz,Izz = 0,0,0,0,0,0
        for row in range(len(all_atm_rec)):
            Ixx += ((float(all_atm_rec[row][5])) - com_x) * ((float(all_atm_rec[row][5])) - com_x)
            Ixy += ((float(all_atm_rec[row][5])) - com_x) * ((float(all_atm_rec[row][6])) - com_y)
            Ixz += ((float(all_atm_rec[row][5])) - com_x) * ((float(all_atm_rec[row][7])) - com_z)
            Iyy += ((float(all_atm_rec[row][6])) - com_y) * ((float(all_atm_rec[row][6])) - com_y)
            Iyz += ((float(all_atm_rec[row][6])) - com_y) * ((float(all_atm_rec[row][7])) - com_z)
            Izz += ((float(all_atm_rec[row][7])) - com_z) * ((float(all_atm_rec[row][7])) - com_z)

        Ixx= Ixx/row
        Iyy= Iyy/row
        Izz= Izz/row
        Ixy= Ixy/row
        Ixz= Ixz/row
        Iyz= Iyz/row

        gyration_tensor = [[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]]
#       print(gyration_tensor)
        evals, evecs = np.linalg.eig(gyration_tensor)
        L1 = evals[0]
        L2 = evals[1]
        L3 = evals[2]

        asphr = ((L1 - L2)**2 + (L2 - L3)**2 + (L1 - L3)**2)/(2.0*((L1 + L2 + L3)**2))

    # Calc masses for all
    last_rnum = 9999999
    for row in range(len(model_array)):
        if (model_array[row][1].strip() == 'CA'):
            rnam = model_array[row][2]
            aa_mass, vol_aa, vbar_aa, sc_rad_aa = AA_data[rnam]
            numerator += (aa_mass * vbar_aa)
            prot_mol_mass += aa_mass
            AA = AA + 1

        # DNA, RNA, Saccharides, Detergents
        else:
            rnam = model_array[row][2].strip()
            chid = model_array[row][3]
            rnum = model_array[row][4]
            if rnum != last_rnum:
                for rec in  NA_data.items():
                    if rnam == rec[0]:
                        aa_mass, vol_aa, vbar_aa, sc_rad_aa = NA_data[rnam]
                        numerator += (aa_mass * vbar_aa)
                        prot_mol_mass += aa_mass
                        last_rnum = rnum
                        NA = NA +1
                        # Check for nucleotide instead of nucleoside
                        if (model_array[row][1] == "P"):
                            numerator += (62.97 * 0.501)
                            prot_mol_mass += 62.97
                for rec in  GL_data.items():
                    if rnam == rec[0]:
                        aa_mass, vol_aa, vbar_aa, sc_rad_aa = GL_data[rnam]
                        numerator += (aa_mass * vbar_aa)
                        prot_mol_mass += aa_mass
                        last_rnum = rnum
                        GL = GL +1
                for rec in  DT_data.items():
                    if rnam == rec[0]:
                        aa_mass, vol_aa, vbar_aa, sc_rad_aa = DT_data[rnam]
                        numerator += (aa_mass * vbar_aa)
                        prot_mol_mass += aa_mass
                        last_rnum = rnum
                        DT = DT +1

    # Add weight of water for N-term and C-term of protein
    if AA > 0:
        prot_mol_mass += 18.0
    # If nucleic acid add weight of MG and MN
    if NA > 0:
        ion_mass = (num_MG * 24.305) + (num_MN * 54.938)
        prot_mol_mass += ion_mass

    # Correct by -0.0025 because above volumes were measured at 25 deg C and want 20 deg C
    #  as per Svedberg, The Ultracentrifuge, 1940, Oxford Univ. Press
    vbar_prot = (numerator / prot_mol_mass) - 0.0025

    # Get convex hull area and volume
    coords = get_coords(model_array)
    area_hull, vol_hull, Dmax = HullVolume(coords)

    # Calculate shell volume from shell thickness
    # Hydration shell thickness is empirically determined to be optimal
    # This expands each hull plane by hydration thickness
    vol_shell_wat = area_hull * 2.8

    # Calculate total volume of hydrated convex hull (including shell waters)
    vol_hyd_hull = vol_hull + vol_shell_wat

    # Estimate axial ratio
    # Minus 3 Ang because DNA duplex a should be rod length, not diagonal
    #   of hull end vertices (Also necessary to make apoferritin axial ratio = 1)
    a = (Dmax/2.0) - 3.0
    b = math.sqrt((3.0 * vol_hull)/(4.0*math.pi*a))

    # But some very spherical structures may not have the diagonal a full 3 Ang longer.
    # Can't have a < b.
    if a > b:
        # Translational Shape Factor
        _numerator = math.sqrt(1.0 - (b/a)**2)
        _denominator = ((b/a)**0.66666667) * math.log((1 + math.sqrt(1.0 - (b/a)**2))/(b/a))
        Ft = _numerator/_denominator
    else:
        a = b
        Ft = 1.0

    # Weight shape factor
    # Empirically found to work better with expanded volume
    #  (Many combinations tried)
    Ft = math.sqrt(Ft)

    # Axial ratio of prolate ellipsoid of same volume as convex hull
    a_b_ratio = a/b

    # Find radius of sphere of same volume as hydrated convex hull
    factor = 3.0/(4.0 * math.pi)
    Rht = (factor * vol_hyd_hull)**0.333333
    # Include Shape factor to give effective hydrodynamic translational radius
    Rht = Rht * Ft
    # Rht comes as Angstrom
    # need meters in equation below
    Rh_trans = Rht * 1e-8

    # Calculate Svedberg coeff.
    eta = 0.0100194 # poise
    rho = 0.998234 # g/ml water density
    fT = 6.0*math.pi*eta*Rh_trans
    s = prot_mol_mass * (1.0 - (vbar_prot * rho)) / (6.02214e23 * fT)

    # Calculate translational diffusion coeff.
    #R = 8.314e7
    kB = 1.381e-16
    T = 293.15
    Dt = kB*T / fT

    # Calculate fT/fo
    # Vbar is in ml/g, need A^3/Dalton
    vol_prot = prot_mol_mass*vbar_prot/0.60224
    Ro = ((3.0 * vol_prot) / (4.0 * math.pi))**0.333333
    ffo_hyd_P = Rht/Ro

    # Calculate rotational diffusion coeff.
    # Rotation more strongly affected by hydration, Halle & Davidovic, 2003
    # So increase hydration shell thickness
    # Hydration shell thickness empirically found to be optimal
    vol_shell_wat_rot = area_hull * 4.3
    vol_hyd_hull_rot = vol_hull + vol_shell_wat_rot
    Rhr = (factor * vol_hyd_hull_rot)**0.333333

    # Include empirical correction for shape factor to give effective
    #  hydrodynamic rotational radius
    # This obtained by fitting to DNA duplexes
    Fr = Ft**4.0
    Rhr = Rhr * Fr
    # Rhr comes as Angstrom
    # need meters in equation below
    Rh_rot = Rhr * 1e-8
    fR = 8.0*math.pi*eta*(Rh_rot**3.0)
    Dr = kB*T / fR
    tauC = (4.0 * math.pi * eta * (Rh_rot**3.0)) / (3.0 * kB * T)
    tauC = tauC * 1e9

    # Calculate intrinsic viscosity
    # From Einstein viscosity equation,
    # [n] = (2.5*6.02214e23*(4*pi*Rht^3)/3)/prot_mol_mass
    #
    # rearrange to
    #
    #int_vis = (10.0 * math.pi * 0.602214 * ((Rht)**3.0)) / (3.0 * prot_mol_mass)
    ## Deleted.
    ## The above equation is valid only for spherical particles.
    ## Intrinsic viscosity is very sensitive to axial ratio and whether the molecule approximates
    ##   a prolate or oblate ellipsoid.
    ## Decided that int_vis was beyond the scope of HullRad.
    ## Leaving this here for future development.
    int_vis = 0.0

    return s,Dt,Dr,vbar_prot,Rht,ffo_hyd_P,prot_mol_mass,Ro,Rhr,int_vis,a_b_ratio, \
           Ft,Rg,Dmax,tauC,asphr,AA,NA,GL,DT,useNumpy


def HullVolume(coords):
    #
    # Uses qconvex to calculate convex hull
    #
    #   coords are coords of reduced atom model
    #
    if useScipy:
        #Call scipy's ConvexHull using coords
        #This will be an object with attributes that we want
        convex_hull = ConvexHull(coords)

        # Area of convex hull; this is a float
        area_hull = convex_hull.area

        # Volume of convex hull; this is a float
        vol_hull = convex_hull.volume

        # Find max distance between vertices for Dmax calculation
        # Necessary for length of corresponding prolate ellipsoid of revolution
        # Note that the "vertices" from ConvexHull are indices, i.e. pointers to
        # the coordinates in "coords" that serve as the vertices.
        vertices = convex_hull.vertices
        Dmax = 0.0
        for r1 in vertices:
            for r2 in vertices:
                dist = math.sqrt((coords[r1][0] - coords[r2][0])**2.0 + \
                                 (coords[r1][1] - coords[r2][1])**2.0 + \
                                 (coords[r1][2] - coords[r2][2])**2.0 )
                if dist > Dmax:
                    Dmax = dist

        return area_hull, vol_hull, Dmax
    else:
        # Input for qconvex, etc
        dim = 3
        num_points = len(coords)
        instring = " %d\n " % (dim )
        instring += " %d\n " % (num_points )
        for i in range(len(coords)):
            instring += ' %f %f %f\n' % (coords[i][0],coords[i][1],
                                       coords[i][2])
        instring += '\n'

        # Call qconvex with argument, FA (compute total area and volume)
        ### Edit the path to qconvex ###
        p = subprocess.Popen(["/opt/local/bin/qconvex", "FA"],
            universal_newlines=True,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        (input, output) = (p.stdin, p.stdout)
        input.write(instring)
        input.close()
        data = output.read()

        # Area of convex hull
        area_start = data.find("area:")
        area_end = data.find("Total volume:")
        if area_end == -1:
            area_end = data.find("Approximate volume:")
        area_hull = float(data[(area_start + 5) : area_end] )

        # Volume of convex hull
        vol_start = data.find("volume:")
        if vol_start == -1 :
            print(data)
            return 0
        vol_hull = float(data[vol_start + 7 :] )

        # Find max distance between vertices for Dmax calculation
        # Necessary for length of corresponding prolate ellipsoid of revolution
        # Run qconvex with different argument, p (returns vertex coordinates)
        p = subprocess.Popen(["/opt/local/bin/qconvex", "p"],
            universal_newlines=True,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        (input, output) = (p.stdin, p.stdout)
        input.write(instring)
        input.close()
        data = list(map(string.split, output.readlines()))
        Dmax = 0.0
        for i in range(2,len(data)):
            for j in range(3,len(data)):
                dist = math.sqrt((float(data[i][0])-float(data[j][0]))**2.0 + \
                                 (float(data[i][1])-float(data[j][1]))**2.0 + \
                                 (float(data[i][2])-float(data[j][2]))**2.0 )
                if dist > Dmax:
                    Dmax = dist

        return area_hull, vol_hull, Dmax


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(USAGE)
        sys.exit()

    file = sys.argv[1]

    # Print actual file name
    print((' %s' %  (file)))

    str_file = str(file)
    extension = str_file[-3:]
    basename = str_file[:-4]

    # Make the reduced atom model
    if extension == 'pdb' or extension == 'ent':
        all_atm_rec,num_MG,num_MN,model_array = model_from_pdb(file)
    elif extension == 'cif':
        all_atm_rec,num_MG,num_MN,model_array = model_from_cif(file)

    # Write out model
    # Uncomment next line if you want the unified atom SC model written out
#   write_pdb(model_array, basename + '_model.pdb')

    # Call the main function
    s,Dt,Dr,vbar_prot,Rht,ffo_hyd_P,M,Ro,Rhr,int_vis,a_b_ratio,Ft,Rg,Dmax,tauC, \
           asphr,AA,NA,GL,DT,useNumpy  = Sved(all_atm_rec,num_MG,num_MN,model_array)

    # Print coefficients to screen
    print(('  #Amino Acids    :    %9.0f' % (AA)))
    print(('  #Nucleotides    :    %9.0f' % (NA)))
    print(('  #Saccharides    :    %9.0f' % (GL)))
    print(('  #Detergents     :    %9.0f' % (DT)))
    print(('  M               :       %6.0f     g/mol' % (M)))
    print(('  v_bar           :       %6.3f     mL/g'  % (vbar_prot)))
    print(('  Ro(Anhydrous)   :       %6.2f     Angstroms' % (Ro)))
    print(('  Rg(Anhydrous)   :       %6.2f     Angstroms' % (Rg)))
    print(('  Dmax            :       %6.2f     Angstroms' % (Dmax)))
    print(('  Axial Ratio     :       %6.2f' % (a_b_ratio)))
    print(('  f/fo            :       %6.2f'  % (ffo_hyd_P)))
    print(('  Dt              :       %6.2e   cm^2/s' % (Dt)))
    print(('  R(Translation)  :       %6.2f     Angstroms' % (Rht)))
    print(('  s               :       %6.2e   sec' % (s)))
    if a_b_ratio > 2.63:
        print(' ')
        print('  Caution. Axial ratio too large for accurate prediction')
        print('  of following rotational properties by HullRad.')

    print(('  Dr              :       %6.2e   s^-1' % (Dr)))
    print(('  R(Rotation)     :       %6.2f     Angstroms' % (Rhr)))
    print(('  tauC            :       %6.2f     ns (from R_rotation)' % (tauC)))
    if useNumpy:
        print(('  Asphericity     :       %6.2f     (from Gyration Tensor)' % (asphr)))
