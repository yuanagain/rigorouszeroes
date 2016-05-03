from interval import interval, inf, imath
from complexinterval import ComplexInterval, _one, _zero, _im
from complexpolynomial import ComplexPolynomial
from simpson import Simpson
from complexrational import ComplexRational
from domain import RectDomain, neighborhood, fromInterval
from newton import Newton
import math

class Algo:
	def __init__(self, poly, domain, res=10**-6):
		self.poly = poly
		self.domain = domain
		self.res = res
		self.zeros = []

	def getZeroes(self):
		"""
		Computes all zeroes in this domain. Returns a list of 
		zero, multiplicity tuples.
		"""
		z_count = argument_principle(self.poly, self.domain)
		return self.search(self.poly, self.domain, z_count)

	def search(self, f, domain, z_count):
		"""
		Searches our rectangular domain for zeroes. We are given 
		that there are z_count zeroes of f in this domain.
		"""

		if z_count <= 0:
			return []

		## if at resolution
		if domain.toInterval().max_dim() < self.res:
			return [(domain, z_count)]

		zeroes = []
		print("Zeroes: " + str(z_count))
		if z_count == 1:
			newton = Newton(domain.midpoint(), f)
			if newton.iterate_until(self.res * 10**-6, 50):
				# Verify that it is a zero
				if argument_principle(f, fromInterval(newton.step) ) == 1:
					return [(newton.next, 1)]

				# Expand the ball a little bit to resolution
				# to attempt to capture the zero
				ball = neighborhood(newton.step, self.res)
				if argument_principle(f, fromInterval(ball) ) == 1:
					return [(ball, 1)]

		else:
			# attempt Newton search
			newton = Newton(domain.midpoint(), f)
			if newton.iterate_until(self.res, 50):
				# TODO verify that it is a zero
				if argument_principle(f, fromInterval(newton.step) ) == z_count:
					return [(newton.next, z_count)]

				# Expand the ball a little bit to resolution
				# to attempt to capture the zero
				ball = neighborhood(newton.step, self.res)
				if argument_principle(f, fromInterval(ball) ) == z_count:
					return [(ball, z_count)]

			# Newton failed, bisect and repeat

			domains = domain.bisect()
			for d in domains:
				z_ct = argument_principle(f, d)
				zeroes += self.search(f, d, z_ct)

		return zeroes

def argument_principle(f, domain):
	"""
	Attempts to apply the argument principle to a rectangular domain 
	in order to count the number of zeroes of our polynomial that 
	exist in tje domain. Unexpected behavior likely result of the
	existence of a zero on the domain boundary.
	"""
	return firstInt(integrate(f, domain))

def integrate(f, domain):
	summs = _zero()
	edges = domain.edges()
	g = ComplexRational(f.derive(), f)
	for i in range(len(edges)):

		summand = 0
		if (i >= 2):
			summand = -Simpson(g, edges[i])
		else:
			summand = Simpson(g, edges[i])
		summs += summand
		print(summand)
	summs = summs / _im(math.pi * 2)
	return summs

def firstInt(intvl):
	"""
	Returns real part of the first integer in the interval.
	Returns None if no such integer exists
	"""
	val = int(math.ceil(min(min(intvl.a))))
	if (val not in intvl.a):
		return val # return None once remainders are fixed
	return val

def main():
	print("Testing Algo")
	print("----------------------------")

	xa = interval([-1, -2])
	xb = interval([-5, -6])
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

	domain = RectDomain(x, z)
	#print(domain)
	domain_2 = RectDomain(x, y)
	print("Domain:")
	print(domain_2)
	print("")
	domain_3 = RectDomain(z, w)
	#print(domain_3)

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

	a_0_a = interval([1, 1])
	a_0_b = interval([5, 5])
	a_0 = ComplexInterval(a_0_a, a_0_b)

	a_1_a = interval([1, 1])
	a_1_b = interval([5, 5])
	a_1 = ComplexInterval(a_1_a, a_1_b)

	a_2_a = interval([3, 3])
	a_2_b = interval([2, 2])
	a_2 = ComplexInterval(a_2_a, a_2_b)

	a_3_a = interval([7, 7])
	a_3_b = interval([-4, -4])
	a_3 = ComplexInterval(a_3_a, a_3_b)

	a_4_a = interval([-6, -6])
	a_4_b = interval([1, 1])
	a_4 = ComplexInterval(a_4_a, a_4_b)

	a_5 = ComplexInterval(interval([4]), interval([0]))
	a_6 = ComplexInterval(interval([5]), interval([0]))
	a_7 = ComplexInterval(interval([5]), interval([0]))


	poly_1 = ComplexPolynomial([a_1, a_2, a_3, a_4, a_5, a_6])

	print("Function:")
	print(poly_1)
	#print("Zeroes: " + str(argument_principle(poly_1, domain_2 )))

	print("Function:")
	poly_2 = ComplexPolynomial([_zero(), a_4])
	print(poly_2)
	#print("Zeroes: " + str(argument_principle(poly_2, domain_2 )))

	print("Function:")
	poly_3 = ComplexPolynomial([a_5, a_6, a_3, a_1, a_0, a_3])
	print(poly_3)
	#print("Zeroes: " + str(argument_principle(poly_3, domain_2 )))

	print("Function:")
	poly_4 = ComplexPolynomial([a_5, a_6, a_3])
	print(poly_4)
	#print("Zeroes: " + str(argument_principle(poly_4, domain_2 )))

	print("Function:")
	poly_5 = ComplexPolynomial([a_5])
	print(poly_5)
	#print("Zeroes: " + str(argument_principle(poly_5, domain_2 )))

	# print(domain_2)
	# algo_1 = Algo(poly_1, domain_2)
	# print(algo_1.getZeroes())

	algo_2 = Algo(poly_2, domain_2)
	zeroes_2 = algo_2.getZeroes()
	print(zeroes_2)
	print(len(zeroes_2))


	# algo_3 = Algo(poly_2, domain_2)
	# print(algo_3.getZeroes())

	#algo_5 = Algo(poly_5, domain_2)
	#print(algo_5.getZeroes())

	# print("+++++++++++===========+++++++++++")
	# print(poly_2)
	# h = ComplexRational(poly_2.derive(), poly_2)
	# print(h)
	# print(domain_2)
	# print('side integrals')
	# print(integrate(poly_2, domain_2))

if __name__=="__main__":
	main()

