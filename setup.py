from glob import glob
from os.path import basename, splitext
# from setuptools import find_packages, setup
from setuptools import setup


setup(
    name = 'svnds',
    version = '0.0.1',
    packages = ['svnds'],
    pakcage_dir = {'svnds': 'svnds/'},
    package_data = {'svnds': ['data/*']},
    py_modules = [splitext(basename(path))[0] for path in glob('svnds/*.py')]
)