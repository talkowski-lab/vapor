from setuptools import setup
from distutils.core import setup,Extension
from Cython.Build import cythonize
import os
import sys
setup(
   name="VaPoR",
   version="0.0.1",
   author="Xuefang Zhao, University of Michigan",
   author_email="xuefzhao@umich.edu",
   description="PacBio long read based genomic structural variants validator.",
   packages=["vapor_vali"],
   scripts=["vapor_vali/vapor"],
   package_data={
       "vapor_vali": [
           "templates/pred_config",
           "templates/truth_config",
           "contrib/SMCScoring.py",
       ],
   },
    # language_level 2: the scoring code relies on Cython 0.29's default (Python 2) semantics
    ext_modules = cythonize("vapor_vali/*.pyx", compiler_directives={"language_level": 2}),
     install_requires=[
     'cython<3', 'numpy','scipy','matplotlib','scikit-learn','pysam'
      ],
    license="Propriety",
    keywords="Long read",
    classifiers=[
       "Development Status :: 1 - Planning",
       "Environment :: Console",
       "Intended Audience :: Information Technology",
       "Intended Audience :: Science/Research",
       "License :: Other/Proprietary License",
       "Natural Language :: English",
       "Programming Language :: Python :: 2",
       "Programming Language :: Python :: 2 :: Only",
       "Topic :: Scientific/Engineering :: Bio-Informatics",
   ],
)

