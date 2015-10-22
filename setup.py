try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Convert MAPS HDF5 files to Scientific Data Exchange format',
    'author': 'David J. Vine',
    'url': 'github.com/djvine/maps_data_exchange',
    'download_url': 'https://github.com/djvine/maps_data_exchange.git',
    'author_email': 'djvine@gmail.com',
    'version': '0.1',
    'packages': ['MAPSToTomoPy'],
    'scripts': [],
    'name': 'MAPSToTomoPy'
}

setup(**config)