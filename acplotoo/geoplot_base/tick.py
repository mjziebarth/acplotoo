# A class to handle axes ticks.

class Tick:
	"""
	Tick-handling class.
	"""

	def __init__(self, x, y, lon, lat, tick_type='lon', ax=None,
	             renderer=None, label_sign='label', use_latex=True):
		"""
		x,y        : Projected coordinates of tick.
		lon,lat    : Geographic coordinates of tick
		tick_type  : Which coordinate changes at this tick.
		             One of 'lon' or 'lat'.
		ax         : Axes instance.
		renderer   : Where to render.
		label_sign : One of 'label' or 'sign'.
		use_latex  : Whether to create labels using latex syntax.
		"""
		self._label = None
		self._label_text = None
		self._width = None
		self._height = None
		self._ax = ax
		self._renderer = renderer
		self._label_sign = label_sign
		self._use_latex = use_latex

		self._tick_type = tick_type
		if tick_type == 'lon':
			self._value = lon
		elif tick_type == 'lat':
			self._value = lat

		# Save coordinates:
		self.x = x
		self.y = y
		self.lon = lon
		self.lat = lat


	def __repr__(self):
		return "Tick(" + self.label_text() + ")"


	def tick_type(self):
		"""
		Returns the tick type: Which of longitude and latitude
		is represented by this tick.
		"""
		return self._tick_type


	def value(self):
		"""
		Returns the value of this tick in arcdegrees.
		"""
		return self._value


	def label_text(self):
		if self._label_text is None:
			# Create label:
			self._create_label_text()

		return self._label_text


	def invalidate(self):
		"""
		Invalidates the label, i.e. forces the creation
		of a new text.
		"""
		if self._label is not None:
			del self._label
			self._label = None


	def label(self):
		"""
		
		"""
		if self._label is None:
			# Create label:
			self._create_label()

		return self._label


	def set_position(self, x, y):
		del self._label
		self._label = None
		self.label().set_position((x,y))
		self.label().set_zorder(10)


	def width(self):
		"""
		Return width relative to axes width.
		"""
		if self._width is None:
			self._width = self.label().get_window_extent(self._renderer) \
			                  .transformed(self._ax.get_figure()
			                                   .dpi_scale_trans.inverted()).width

		width = self._width

		return width


	def height(self):
		"""
		Return height relative to axes height.
		"""
		if self._height is None:
			self._height = self.label().get_window_extent(self._renderer) \
			            .transformed(self._ax.get_figure()
			                         .dpi_scale_trans.inverted()).height

		return self._height


	def _create_label_text(self, label_sign, use_latex):
		# Calculate degrees, minutes, seconds:
		value = self._value
		sign = 1 if self._value > 0 else -1
		value *= sign
		degrees = round(value)
		value -= degrees
		minutes = round(value / 60.0)
		value -= 60 * minutes
		seconds = round(value / 3600.0)

		# Base label:
		prefix = "$" if use_latex else ""
		degree_symbol = "^\\circ" if use_latex else "°"
		minute_symbol = "'"
		second_symbol = "\""
		space = "\\," if use_latex else " "
		suffix = "$" if use_latex else ""

		if seconds == 0 and minutes == 0:
			label = str(degrees) + degree_symbol
		elif seconds == 0:
			label = str(degrees) + degree_symbol + space + str(minutes) \
			         + minute_symbol
		else:
			label = str(degrees) + degree_symbol + space + str(minutes) \
			          + minute_symbol + space + str(seconds) + second_symbol

		if label_sign == 'label':
			if self._tick_type == 'lon':
				if sign < 0:
					label = prefix + label + space + \
					        ("\\mathrm{W}" if use_latex else "W") + suffix
				else:
					label = prefix + label + space + \
					        ("\\mathrm{E}" if use_latex else "E") + suffix
			else:
				if sign < 0:
					label = prefix + label + space + \
					        ("\\mathrm{S}" if use_latex else "S") + suffix
				else:
					label = prefix + label + space + \
					        ("\\mathrm{N}" if self._use_latex else "N") + suffix
		elif label_sign == 'sign':
			if sign == -1:
				label = prefix + "-" + label + suffix

		self._label_text = label


	def _create_label(self):
		# Make sure an Axes instance exists:
		if self._ax is None:
			raise RuntimeError("Axes instance has to be given to create label!")

		# Make sure label text is created:
		if self._label_text is None:
			self._create_label_text(self._label_sign, self._use_latex)

		# Now create label:
		self._label = self._ax.text(0.5*sum(self._ax.get_xlim()),
		                            0.5*sum(self._ax.get_ylim()),
		                            self._label_text)
