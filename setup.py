from glob import glob
from os.path import basename, splitext
# from setuptools import find_packages, setup
from setuptools import setup, find_packages
from svnds import __version__


version = __version__
setup(
    name = 'svnds',
    version = version,
    author = 'Eunhee Ko',
    packages = ['svnds', 'svnds/data', 'svnds/data/elcosmos'],
    package_data = {'svnds/data/elcosmos': ['*.csv']},
    include_package_data = True,
    py_modules = [splitext(basename(path))[0] for path in glob('svnds/*.py')]
)


# package_dir = {'': 'svnds'},