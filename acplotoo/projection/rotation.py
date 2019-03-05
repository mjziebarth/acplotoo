# Rotation wrapper over other projection.

import numpy as np

from .projection import Projection

class Rotation(Projection):
	"""
	This projection is a simple rotation of another projection.
	"""

	def __init__(self, projection, angle):
		# Make sure we have been passed a projection:
		if not isinstance(projection, Projection):
			raise TypeError("Rotation needs to be passed another "
			                "projection!")

		# So far, support only right angles:
		angle = int(angle)
		if not angle in (0, 90, 180, 270, -90):
			raise ValueError("Only right angle rotations (0°, 90°, 180°, 270°/90°) "
			                 "are supported!")

		# Save data:
		self._projection = projection
		self._angle = angle


	def __hash__(self):
		return hash(("Rotation",angle, self._projection))


	### The interface implementation: ###

	def _project(self, lon, lat):
		"""
		Class implementation of _project method.
		"""
		x,y = self._projection._project(lon, lat)
		if self._angle == 0:
			return x,y
		elif self._angle == 90:
			return y,-x
		elif self._angle == 180:
			return -x, -y
		elif self._angle == 270 or self._angle == -90:
			return -y,x


	def _inverse(self, x, y):
		"""
		Class implementation of _inverse method.
		"""
		if self._angle == 0:
			return self._projection._inverse(x,y)
		elif self._angle == 90:
			return self._projection._inverse(-y,x)
		elif self._angle == 180:
			return self._projection._inverse(-x,-y)
		elif self._angle == 270 or self._angle == -90:
			return self._projection._inverse(y, -x)


	def _is_global(self):
		return self._projection._is_global()


	def _identifier(self):
		return ("Rotation", self._angle, self._projection._identifier())


	def __eq__(self, other):
		if isinstance(other, Rotation):
			return self._identifier() == other._identifier()

		return False
