# -*- coding: utf-8 -*-
#
# Helper methods for the acplotoo.sphereplot Python module.
# Copyright (C) 2018 Malte Ziebarth
# 
# This software is distributed under the MIT license.
# See the LICENSE file in this repository.
import numpy as np

def connect_masked_sequence(data, mask):
	masked = np.argwhere(~mask)
	if mask[0] and mask[-1]:
		indices = [i+masked.flat[-1]+1 for i in
		           range(len(data)-masked.flat[-1]-1)] \
		         + [i for i in range(masked.flat[0])]
	else:
		indices = np.argwhere(mask)
	data = data[indices]
	return data.reshape((data.size,))
