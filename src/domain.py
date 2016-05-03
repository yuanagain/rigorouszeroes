"""
This module provides tools to representing domains
on the complex plane.

Author: 
Yuan Wang
Princeton University
"""
#import matplotlib.pyplot as plt
from interval import interval, inf, imath
from complexinterval import ComplexInterval, _one, _zero, _hull, _im, _real

class RectDomain():
	def __init__(self, x, y):
		## or use interval.hull
		self.a = interval.hull((x.a, y.a))
		self.b = interval.hull((x.b, y.b))
		
		self.real_max = max(max(self.a))
		self.real_min = min(min(self.a))


		self.im_max = max(max(self.b))
		self.im_min = min(min(self.b))

	def toInterval(self):
		"""
		Returns rectangular domain as an equivalent complex interval
		"""
		return ComplexInterval(self.a, self.b)

	def extrema(self):
		return self.toInterval().extrema()

	def __contains__(self, z):
		"""
		Returns true iff complex interval z is contained
		within this complex rectangular domain
		"""
		if (z.a in self.a):
			if (z.b in self.b):
				return True

		return False

	def edges(self):
		"""
		Returns the boundary intervals for the rectangular domain
		"""
		edges = []

		topright = ComplexInterval(interval([self.real_max]), interval([self.im_max]))
		topleft = ComplexInterval(interval([self.real_min]), interval([self.im_max]))
		bottomleft = ComplexInterval(interval([self.real_min]), interval([self.im_min]))
		bottomright = ComplexInterval(interval([self.real_max]), interval([self.im_min]))

		edges.append(_hull(topright, topleft))
		edges.append(_hull(topleft, bottomleft))
		edges.append(_hull(bottomleft, bottomright))
		edges.append(_hull(bottomright, topright))

		return edges

	def midpoint(self):
		"""
		Returns the midpoint of the complex
		rectangular domain.
		"""
		return ComplexInterval(self.a.midpoint, self.b.midpoint)

	def bisect(self, perturbation = 0.0):
		"""
		Bisects this domain along its longest edge, whichever is greater
		"""
		asInterval = self.toInterval()
		if asInterval.im_radius() > asInterval.re_radius():
			return self.bisect_horizontal(perturbation)

		return self.bisect_vertical(perturbation)

	def bisect_horizontal(self, perturbation = 0.0): 
		"""
		Bisect the complex rectangular domain horizontally, 
		i.e. with the bisecting line parallel to the 
		real axis. Perturb with the given 
		perturbation (default = 0). Returns the two 
		resultant rectangular domains.
		"""
		a_1 = interval([self.real_min, max(max(self.a.midpoint)) + perturbation])
		a_2 = interval([max(max(self.a.midpoint)) + perturbation, self.real_max])
		interval_left = ComplexInterval(a_1, self.b)
		interval_right = ComplexInterval(a_2, self.b)
		return RectDomain(interval_left, interval_left), RectDomain(interval_right, interval_right)

	def bisect_vertical(self, perturbation = 0.0): 
		"""
		Bisect the complex rectangular domain vertically, 
		i.e. with the bisecting line parallel to the 
		imaginary axis. Perturb with the given 
		perturbation (default = 0). Returns the two 
		resultant rectangular domains.
		"""
		b_1 = interval([self.im_min, max(max(self.b.midpoint)) + perturbation])
		b_2 = interval([max(max(self.b.midpoint)) + perturbation, self.im_max])
		interval_bottom = ComplexInterval(self.a, b_1)
		interval_top = ComplexInterval(self.a, b_2)
		return RectDomain(interval_bottom, interval_bottom), RectDomain(interval_top, interval_top)

	def show(self):
		"""
		Plots the domain's boundaries
		"""
		print("Not yet implemented")
		# plt.ylabel('imaginary')
		# plt.xlabel('real')
		# plt.plot([2,2],[2,4])

	def __str__(self):
		return "RectDomain( " + self.a.__str__() + ' , ' + self.b.__str__()

def neighborhood(z, diameter = 10**-6):
	"""
	Returns a neighborhood around the midpoint of z of a given diameter
	"""
	mid = z.midpoint()
	topright = z.midpoint() + _real(diameter / 2) + _im(diameter / 2)
	bottomleft = z.midpoint() - _real(diameter / 2) - _im(diameter / 2)
	return RectDomain(topright, bottomleft)

def fromInterval(c_int):
	"""
	Returns a Rect_Domain corresponding to a given complex interval
	"""
	extrema = c_int.extrema()
	return RectDomain(extrema[0], extrema[2])

def main():
	print("Testing Rectangular Domains")
	print("----------------------------")
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
	print("Testing Constructor")
	print("----------------------------")
	domain = RectDomain(x, z)
	print(domain)
	domain_2 = RectDomain(x, y)
	print(domain_2)
	domain_3 = RectDomain(z, w)
	print(domain_3)
	
	print("============================")
	print("Testing Horizontal Bisection")
	print("----------------------------")
	print(domain.bisect_horizontal()[0])
	print(domain.bisect_horizontal()[1])

	print("============================")
	print("Testing Vertical Bisection")
	print("----------------------------")
	print(domain.bisect_vertical()[0])
	print(domain.bisect_vertical()[1])


	print("============================")
	print("Testing Edges")
	print("----------------------------")
	print(domain)
	edges = domain.edges()
	for edge in edges:
		print(edge)

if __name__=="__main__":
	main()