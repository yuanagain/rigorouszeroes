"""
Implements a complex rational function F(z) = P(z) / Q(z),
where P and Q are complex polynomials

Author: 
Yuan Wang
Princeton University
"""
from interval import interval, inf, imath
from complexinterval import ComplexInterval, _one, _zero, _real, _im
from complexpolynomial import ComplexPolynomial, _identity, _constant
from domain import RectDomain, fromInterval
import operator
import numpy as np
import math

class ComplexRational: 
	def __init__(self, p, q):
		self.p = p
		self.q = q

	def isZero(self):
		"""
		Returns whether or not this function's is zero.
		"""
		return (self.p.isZero() & (self.q.isZero() == False))

	def derive(self):
		"""
		Returns a rational function that is the 1st derivative of 
		the current rational function
		"""
		## check if zero

		num = self.q * self.p.derive() - self.p * self.q.derive()
		denom = self.q * self.q
		return ComplexRational(num, denom)

	def evaluate(self, z):
		"""
		Evaluates a complex polynomial at z
		"""
		return self.p(z) / self.q(z)

	def __str__(self):
		"""
		Returns a string representation of the polynomial
		"""
		return 'f(z) = ' + self.p.coeffString() + ' / ' + self.q.coeffString()

	def __call__(self, z):
		return self.evaluate(z)

def main():
	print("Testing Complex Rational Functions")
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

	za = interval([4])
	zb = interval([5])	
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

	coeffs = [a_1, a_2, a_3, a_4]
	poly_1 = ComplexPolynomial(coeffs)
	print(poly_1)
	poly_2 = ComplexPolynomial([_zero(), a_4])
	print(poly_2)
	poly_4 = ComplexPolynomial([a_5, a_6])
	poly_5 = ComplexPolynomial([a_7, a_5])
	print(poly_4)
	print(poly_5)
	print(poly_4 * poly_5)
	print("============================")
	
	print("Construction")
	print("----------------------------")
	rat_1 = ComplexRational(poly_2, poly_4)
	print(rat_1)
	rat_2 = ComplexRational(poly_1, poly_1.derive())
	print(rat_2)
	rat_3 = ComplexRational(poly_5, _identity())
	print(rat_3)


	print("============================")
	print("Evaluation")
	print("----------------------------")

	print(rat_1(a_1))
	print(rat_2(a_2))
	print(rat_3(a_3))


	print("============================")
	print("Derivation")
	print("----------------------------")
	print(rat_3.derive())
	print(rat_3.derive().evaluate(a_2))

	print("============================")

	print("============================")
	print("Testing if Zero")
	print("----------------------------")
	print(rat_3.isZero())
	print(rat_3.derive().isZero())

	print("============================")
	poly_6 = ComplexPolynomial([_zero(), _real(-6) + _im(1)])
	print(poly_6)
	rat_3 = ComplexRational(poly_5, poly_5)
	print(rat_3)
	d_4 = rat_3.derive()
	print(d_4)
	d_5 = d_4.derive()
	print(d_5)
	print(d_4(z))

if __name__=="__main__":
	main()

