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
    packate_dir = {'': 'svnds'},
    package_data = {'': ['data/*']},
    include_package_data = True,
    py_modules = [splitext(basename(path))[0] for path in glob('svnds/*.py')]
)
