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
    packages = ['svnds', 'svnds/data', 'svnds/data/elcosmos', 'svnds/data/filters/EUCLID','svnds/data/filters/LSST','svnds/data/filters/PANSTARRS', 'svnds/data/filters/SDS','svnds/data/filters/SPHEREx', 'svnds/data/filters/VIKING', 'svnds/data/sensitivity', 'svnds/data/systematics'],
    package_data = {'svnds/data/elcosmos': ['*.csv', '*.py'], 'svnds/data/filters/EUCLID': ['*.dat'], 'svnds/data/filters/LSST': ['*.dat'], 'svnds/data/filters/PANSTARRS': ['*.fit'], 'svnds/data/filters/SDS': ['*'], 'svnds/data/filters/SPHEREx': ["*"], 'svnds/data/filters/VIKING': ["*.dat"], 'svnds/data/sensitivity': ['*'], 'svnds/data/systematics': ['*']},
    include_package_data = True,
    py_modules = [splitext(basename(path))[0] for path in glob('svnds/*.py')]
)


# package_dir = {'': 'svnds'},