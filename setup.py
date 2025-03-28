#!/usr/bin/env python3
from setuptools import setup
from setuptools.command.install import install
import os
import sys

__version__ = '0.1.3'

def readme():
    with open('README.md') as f:
        return f.read()

def check_dir_write_permission(directory):
    if os.path.isdir(directory) and not os.access(directory, os.W_OK):
        sys.exit('Error: no write permission for ' + directory + '  ' +
                 'Perhaps you need to use sudo?')

class EnviroAmpDesignerInstall(install):

    def run(self):
        check_dir_write_permission(self.install_lib)
        install.run(self)



setup(name='EnviroAmpDesigner',
      version=__version__,
      description='EnviroAmpDesigner',
      long_description=readme(),
      python_requires='>=3.12',
      classifiers=['Development Status :: Beta',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Intended Audience :: Science/Research'],
      keywords='PCR amplicon primers',
      url='https://github.com/AntonS-bio/EnviroAmpDesigner.git',
      author='Anton Spadar',
      author_email='',
      packages=['scripts'],
      include_package_data=True,
      entry_points={'console_scripts': ['design_primers = design_primers:main']},
      scripts=[
          'scripts/data_classes.py',
          'scripts/design_primers.py',
          'scripts/generate_msa.py',
          'scripts/hierarchy_utils.py',
          'scripts/identify_genotype_snps.py',
          'scripts/identify_species_snps.py',
          'scripts/inputs_validation.py',
          'scripts/load_vcfs.py',
          'scripts/metadata_utils.py',
          'scripts/name_converters.py',
          'scripts/primers_generator.py',
          'scripts/run_blast.py',
          'scripts/snp_optimiser.py'
      ],
      cmdclass={'install': EnviroAmpDesignerInstall}
)
