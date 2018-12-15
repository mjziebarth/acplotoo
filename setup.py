# -*- coding: utf-8 -*-
#
# Setup script for the plot_tools Python module.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

# Imports:
from setuptools import setup
from setuptools.extension import Extension
from numpy import get_include as np_get_include
from Cython.Distutils import build_ext
from Cython.Compiler.Options import _directive_defaults

_directive_defaults['linetrace'] = True
_directive_defaults['binding'] = True

# Extensions:
extensions = []

extensions.append(Extension('acplotoo.geoplot_base.streamplot',
	sources=['acplotoo/geoplot_base/streamplot.pyx',
	         'src/src/grid.cpp',
	         'src/src/streamplot.cpp'],
	include_dirs=[np_get_include(),'src/include'],
	extra_compile_args=['-std=c++17'],
	language='c++'))


# Setup:

setup(
	name='acplotoo',
	version='1.0.0',
	description="Code snippets for plotting with matplotlib.",
	long_description="A collection of code snippets I find useful for working "
	                 "with matplotlib.",
	author='Malte J. Ziebarth',
	author_email='contact@fmvkb.de',
	packages=['acplotoo','acplotoo.euclidean','acplotoo.sphere',
	          'acplotoo.geoplot_base', 'acplotoo.projection'],
	py_modules=['acplotoo'],
	provides=['acplotoo','acplotoo.euclidean','acplotoo.sphere',
	          'acplotoo.projection'],
	cmdclass = {'build_ext': build_ext},
	ext_modules=extensions,
	scripts=[],
	install_requires=['numpy','matplotlib'],
	license='MIT',
)
