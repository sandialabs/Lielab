from setuptools import find_packages, setup, Extension, Distribution
import platform
import os
import sys

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

with open('requirements.txt') as f:
    requirements_dev = f.read().splitlines()

long_description = '''A lightweight C++ library for numerical computations.'''

modules = ['lielab.domain',
           'lielab.dynamics',
           'lielab.functions',
           'lielab.kinematics',
           'lielab.optim',
           'lielab.testing',
           'lielab.topos',
           'lielab.transform']

if platform.system() == 'Windows':
  bin_ext = '*.pyd'
elif platform.system() == 'Linux':
  bin_ext = '*.so'

dir_setup = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(dir_setup, 'lielab'))
from cppLielab import __author__, __contact__, __location__, __version__

setup(name="lielab",
      version=__version__,
      description=long_description,
      long_description=long_description,
      author=__author__,
      author_email=__contact__,
      platforms=["any"],
      python_requires='>=3.6',
      url=__location__,
      packages=['lielab'] + modules,
      package_data={'lielab' : [bin_ext]},
      install_requires=requirements,
      extras_require={'dev': requirements_dev}
      )