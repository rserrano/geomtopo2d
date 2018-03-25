"""
Copyright 2018 Geomodelr, Inc. 
rserrano at geomodelr.com

This file is part of Geomtopo2d. Geomtopo2d is free software: 
you can redistribute it and/or modify it under the terms of 
the GNU Lesser General Public License as published by the Free
Software Foundation, either version 3 of the License, or
(at your option) any later version.

Geomtopo2d is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License
along with Librtmfp.  If not, see <http://www.gnu.org/licenses/>.
"""

from setuptools import setup
import os

setup(name='geomtopo2d',
    version='0.0.1',
    description='A library with a few geometrical and topological functions.',
    url='https://github.com/rserrano/geomtopo2d',
    author='Ricardo Serrano',
    author_email='rserrano@geomodelr.com',
    license='LGPL',
    packages=['geomtopo2d'],
    install_requires=['numpy', 'scipy'],
    keywords=['geometry', 'topology', '2d'],
    # entry_points = {
    #     'console_scripts': [
    #         'geomodelr=geomodelr.__main__:main'
    #     ]
    # },
    zip_safe=False)

