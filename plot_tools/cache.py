# -*- coding: utf-8 -*-
#
# Plot tools caching file.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.

# Try to import joblib:
try:
	from joblib import Memory
	has_joblib = True
	plot_tools_cache =  Memory(".geoplot", verbose=0)
except:
	has_joblib = False
	plot_tools_cache = None

