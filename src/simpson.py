"""
Implement's Simpson's method on the 
complex plane.
"""
from interval import interval, inf, imath
from complexinterval import ComplexInterval, _one, _zero, _real, _im, _hull
from complexpolynomial import ComplexPolynomial
from complexrational import ComplexRational
from domain import RectDomain

import time
import math

def Simpson(g, 
			line = ComplexInterval(interval([0,1]), interval([0])),
			radius = math.pi/2 - 10**-14,
			max_iterates = 10):
	"""
	Returns an interval containing the integral of 
	our polynomial along the line, where the imaginary 
	part does not have radius exceeding the given radius.

	The line is assumed to be a complex interval with a thin
	real or a thin imaginary component.
	"""
	val = _zero()
	extrema = line.extrema()
	nodes = [extrema[0], line.midpoint(), extrema[2]]
	remainder_radius = radius

	G = getG(g)

	j = 0 # temporary, until remainder bounds are fixed
	while ((remainder_radius >= radius) & (j < max_iterates)):
		j += 1 # temporary, until remainder bounds are fixed

		val = _zero()
		#print("len of nodes = " + str(len(nodes)))
		for i in range(0, len(nodes) - 2, 2):
			val = val + approx(g, nodes[i:i+3])  + remainder(G, nodes[i:i+3])
		
		#print("-------")
		#print(val)
		#time.sleep(1)
		new_nodes = []
		for i in range(len(nodes) - 1):
			new_nodes.append(nodes[i])
			new_nodes.append(_hull(nodes[i], nodes[i + 1]).midpoint())
		new_nodes.append(nodes[-1])

		nodes = new_nodes
		#print(remainder_radius)

	return val

def approx(g, nodes):
	"""
	Computes the approximation of g over the nodes for Simpson's method
	"""
	#print(len(nodes))
	factor = (nodes[2] - nodes[0]) / _real(6)
	_sum = g(nodes[0]) + _real(4) * g(nodes[1]) + g(nodes[2])
	return factor * _sum

def getG(g):
	"""
	Returns the 4th derivative of g
	"""
	G = g

	for i in range(4):
		G = G.derive()
		if G.isZero():
			return G

	return G


def remainder(G, nodes):
	"""
	Computes bounds on the remainder term of Simpson's method
	"""

	factor = (nodes[2] - nodes[0])**5 / _real(2880)

	if G.isZero():
		return _real(0)

	intvl = _hull(nodes[0], nodes[2])
	# print(intvl)
	# print(G)
	return G(intvl) * factor


def main():
	print("Testing Simspon's Method")
	print("----------------------------")

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

	line = ComplexInterval(interval([1.88, 1.99]), interval([4.28]))
	print(line)

	poly_1 = ComplexPolynomial([a_1, a_2, a_3, a_4, a_5, a_6])
	print(poly_1)
	poly_2 = ComplexPolynomial([_zero(), a_4])
	print(poly_2)
	poly_3 = ComplexPolynomial([a_5, a_6, a_3, a_1, a_0, a_3])
	print(poly_3)
	poly_4 = ComplexPolynomial([a_5, a_6, a_3])
	print(poly_4)

	print("======")
	print(poly_4)
	print(poly_4.derive())
	print("=====")
	print("g and g'")
	g = ComplexRational(poly_4.derive(), poly_4)
	print(g)
	print(g.derive())
	print("========")
	print(g.evaluate(line))
	print("#######")
	#print(Simpson(g, line))
	line = ComplexInterval(interval([-2, 7]), interval([3]))

	print("**********")
	poly_2 = ComplexPolynomial([_zero(), a_4])
	h = ComplexRational(poly_2.derive(), poly_2 )
	print(h)
	print(Simpson(h, line))


if __name__=="__main__":
	main()