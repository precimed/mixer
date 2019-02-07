from __future__ import print_function

from setuptools import setup

setup_kwargs = dict(
    name='precimed',
    version='0.9.0',
    packages=['mixer'],
    install_requires=[],  # dependencies
)

if __name__ == '__main__':
    setup(**setup_kwargs)
