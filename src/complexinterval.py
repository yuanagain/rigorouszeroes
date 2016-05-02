"""
Defines a basic complex interval and implements basic operations
on said intervals.

TODO: implement conjugate operations

Author: 
Yuan Wang
Princeton University
"""

from interval import interval, inf, imath

class ComplexInterval: 
	def __init__(self, a = interval([0, 0]), b = interval([0, 0])):
		self.a = a
		self.b = b

		self.real_max = max(max(self.a))
		self.real_min = min(min(self.a))


		self.im_max = max(max(self.b))
		self.im_min = min(min(self.b))

	def extrema(self):
		topright = ComplexInterval(interval([self.real_max]), interval([self.im_max]))
		topleft = ComplexInterval(interval([self.real_min]), interval([self.im_max]))
		bottomleft = ComplexInterval(interval([self.real_min]), interval([self.im_min]))
		bottomright = ComplexInterval(interval([self.real_max]), interval([self.im_min]))

		return [topright, topleft, bottomleft, bottomright]

	def clone(self):
		"""
		Returns a clone of this complex interval
		"""
		return ComplexInterval(self.a, self.b)

	def max_dim(self):
		"""
		Returns the diameter of the longest dimension of our complex interval
		"""
		return max(self.re_radius(), self.im_radius())

	def re_radius(self):
		"""
		Returns the radius of this complex interval's real part
		"""
		hull = interval.hull((self.a, self.a))
		return max(max(hull)) - min(min(hull))

	def im_radius(self):
		"""
		Returns the radius of this complex interval's imaginary part
		"""
		hull = interval.hull((self.b, self.b))
		return max(max(hull)) - min(min(hull))

	def midpoint(self):
		"""
		Returns the midpoint of the complex interval
		"""
		return ComplexInterval(self.a.midpoint, self.b.midpoint)

	def real(self):
		"""
		Returns real part of complex interval
		"""
		return self.a

	def im(self):
		"""
		Returns imaginary part of complex interval
		"""
		return self.b

	def isZero(self):
		"""
		Returns whether or not this interval is zero
		"""
		zero = interval([0])
		return self.a == zero & self.b == zero

	def __contains__(self, z):
		if (z.a in self.a):
			if (z.b in self.b):
				return True
		return False

	def conjugate(self):
		"""
		Returns the conjugate of the complex interval
		"""
		return ComplexInterval(self.a, -self.b)

	def __add__(self, other):
		"""
		Adds two complex intervals
		"""
		return ComplexInterval(self.a + other.a, self.b + other.b)

	def __mul__(self, other):
		"""
		Multiplies two complex intervals
		"""
		real_1 = 0
		if (self.a == other.a):
			real_1 = self.a**2
		else:
			real_1 = self.a * other.a

		real_2 = 0
		if (self.b == other.b):
			real_2 = self.b**2
		else:
			real_2 = self.b * other.b

		real = real_1 - real_2
		im = self.a * other.b + self.b * other.a
		return ComplexInterval(real, im)

	def __sub__(self, other):
		"""
		Subtracts other from self
		"""
		return self.__add__(-other)

	def __neg__(self):
		"""
		Returns negation of interval
		"""
		return ComplexInterval(-self.a, -self.b)

	def __str__(self):
		"""
		Returns string representation of complex interval
		"""
		return self.a.__str__() + ' + ' + self.b.__str__() + 'i'

	def __div__(self, other):
		"""
		Returns quotients of two complex intervals

		TODO: test more thoroughly
		"""
		d = (other.a**2 + other.b**2)
		a = (self.a * other.a + self.b * other.b) / d
		b = (self.b * other.a - self.a * other.b) / d
		return ComplexInterval(a, b)

	def __pow__(self, n):
		"""
		Returns complex interval to the power n
		"""
		out = _one()

		# TODO: implement powers

		for i in range(n):
			out = out * self
		return out

	def __eq__(self, other):
		"""
		Checks for equality between complex intervals
		"""
		return (self.a == other.a & self.b == other.b)

	def __ne__(self, other):
		"""
		Checks for equality between complex intervals
		"""
		return self.__eq__(other) == False

## Functions for generating copies of identities 
def _one():
	"""
	Returns a complex interval representing one
	"""
	return ComplexInterval(interval([1]), interval([0]))

def _zero():
	"""
	Returns a complex interval representing zero
	"""
	return ComplexInterval(interval([0]), interval([0]))

def _real(x):
	"""
	Returns a complex interval representing x + 0i
	"""
	return ComplexInterval(interval([x]), interval([0]))

def _im(y):
	"""
	Returns a complex interval representing 0 + yi
	"""
	return ComplexInterval(interval([0]), interval([y]))

def _hull(w, z):
	"""
	Returns the hull of two complex intervals
	"""
	return ComplexInterval(interval.hull((w.a, z.a)), interval.hull((w.b, z.b)))

def main():
	print("Testing ComplexInterval")
	xa = interval([1, 2])
	xb = interval([5, 6])
	x = ComplexInterval(xa, xb)

	ya = interval([4, 7])
	yb = interval([2, 3])
	y = ComplexInterval(ya, yb)

	wa = interval([2, 2])
	wb = interval([3, 3])
	w = ComplexInterval(wa, wb)

	za = interval([4, 4])
	zb = interval([5, 5])
	z = ComplexInterval(za, zb)

	print("============================")
	print("Construction")
	print("----------------------------")

	print('x = ' + x.__str__())
	print('y = ' + y.__str__())
	print('w = ' + w.__str__())
	print('z = ' + z.__str__())
	print(min(min(x.a)))
	print("============================")
	print("Radius")
	print("----------------------------")
	print(x.im_radius())
	print(z.im_radius())
	print("============================")
	print("Hull")
	print("----------------------------")
	print(_hull(x, y))

	print("============================")
	print("Extrema")
	print("----------------------------")
	extrema = x.extrema()
	for ex in extrema: 
		print(ex)

if __name__=="__main__":
	main()

