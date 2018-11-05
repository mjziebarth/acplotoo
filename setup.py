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

# Setup:

setup(
	name='plot_tools',
	version='1.0.0',
	description="Code snippets for plotting with matplotlib.",
	long_description="A collection of code snippets I find useful for working "
	                 "with matplotlib.",
	author='Malte J. Ziebarth',
	author_email='contact@fmvkb.de',
	packages=['plot_tools','plot_tools.euclidean','plot_tools.sphere',
	          'plot_tools.geoplot_base', 'plot_tools.projection'],
	py_modules=['plot_tools'],
	provides=['plot_tools','plot_tools.euclidean','plot_tools.sphere',
	          'plot_tools.projection'],
	scripts=[],
	license='MIT',
)
