"""
Defines a complex polynomial object from coefficients:
f(z) = a_0 + a_1 * z + a_2 * z^2 + ... + a_n * z^n

where for each i, a_i is a complex interval

Author: 
Yuan Wang
Princeton University
"""
from interval import interval, inf, imath
from complexinterval import ComplexInterval, _one, _zero, _real
import operator
import numpy as np
import math

class ComplexPolynomial: 
	def __init__(self, coeffs = [ComplexInterval()]):
		self.coeffs = coeffs

	def coeffs_copy(self):
		"""
		Return a copy of this complex polynomia's coefficients
		"""
		return self.coeffs[:]

	def derive(self):
		"""
		Returns a polynomial that is the first derivative of 
		the current polynomial
		"""
		mycoeffs = self.coeffs
		coeffs = [0] * (len(mycoeffs) - 1)
		for i in range(len(mycoeffs) - 1):
			factor = _real(i + 1)
			coeffs[i] = mycoeffs[i + 1] * factor
		return ComplexPolynomial(coeffs)

	def evaluate(self, z):
		"""
		Evaluates a complex polynomial at z
		"""
		total = _zero()
		coeffs = self.coeffs
		for i in range(len(coeffs)):
			total = total + coeffs[i] * z**i
		return total

	def coeffString(self):
		out = ""
		for i in range(len(self.coeffs)):
			out += '(' + str(min(min(self.coeffs[i].a))) + ' + ' \
			+ str(min(min(self.coeffs[i].b))) + 'i)' + 'z^' + str(i) + ' + '
		return out[:-2]

	def isZero(self):
		"""
		Returns whether or not this polynomial is zero
		"""
		for i in range(len(self.coeffs)):
			if self.coeffs[i].isZero() == False: 
				return False
		return True


	def __str__(self):
		"""
		Returns a string representation of the polynomial
		"""
		
		if (len(self.coeffs) == 0):
			return "f(z) =  0"

		return "f(z) =  " + self.coeffString()

	def __call__(self, z):
		return self.evaluate(z)

	def __add__(self, other):
		coeffs = self.coeffs[:]
		# add coefficientspairwise
		for i in range(min(len(self.coeffs), len(other.coeffs))):
			coeffs[i] = coeffs[i] + other.coeffs[i]

		if (len(other.coeffs) > len(self.coeffs)):
			return ComplexPolynomial(coeffs + other.coeffs[len(other.coeffs):])
		else:
			return ComplexPolynomial(coeffs)

	def __neg__(self):
		coeffs = [-coeff for coeff in self.coeffs]
		return ComplexPolynomial(coeffs)

	def __sub__(self, other):
		return self.__add__(-other)

	def __mul__(self, poly):
		"""
		Multiplies two polynomials.
		"""
		coeffs = []
		if (self.isZero() | poly.isZero()):
			return ComplexPolynomial([_zero()])

		for i in range(len(self.coeffs) + len(poly.coeffs) - 1):
			coeffs.append(_zero())
		P = self.coeffs
		Q = poly.coeffs
		for k in range(len(P)):
			summand = [q_i * P[k] for q_i in Q]
			for j in range(len(Q)):
				coeffs[k + j] = coeffs[k + j] + summand[j]

		return ComplexPolynomial(coeffs)

def _identity():
	"""
	Returns the monomial f(z) = 1
	"""
	return ComplexPolynomial([ComplexInterval(interval([1]), interval([0]))])

def _constant(k):
	"""
	Returns the monomial f(z) = k
	"""
	return ComplexPolynomial([ComplexInterval(interval([k]), interval([0]))])

def main():
	print("Testing Complex Polynomials")
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

	a_5 = ComplexInterval(interval([4]), interval([1]))
	a_6 = ComplexInterval(interval([5]), interval([2]))
	a_7 = ComplexInterval(interval([5]), interval([3]))

	coeffs = [a_1, a_2, a_3, a_4]

	print("Testing Complex Constructor")
	print("----------------------------")
	poly_1 = ComplexPolynomial(coeffs)
	print(poly_1)
	poly_2 = ComplexPolynomial([_zero(), a_4])
	print(poly_2)
	print("============================")


	print("Testing Evaluation")
	print("----------------------------")
	print(poly_1.evaluate(w))
	print(poly_1.evaluate(_one()))
	print(poly_1.evaluate(_zero()))

	print("")

	print(poly_2.evaluate(w))
	print(poly_2.evaluate(_one()))
	print(poly_2.evaluate(_zero()))
	print("============================")

	print("Derivation")
	print("----------------------------")
	print(poly_1.derive())

	print("")

	print(poly_2.derive())
	print("============================")

	print("Multiplication")
	print("----------------------------")
	poly_4 = ComplexPolynomial([a_5, a_6])
	poly_5 = ComplexPolynomial([a_7, a_5])

	print(poly_4)
	print(poly_5)
	print(poly_4 * poly_5)
	print(poly_4 * poly_4)

	print("============================")

	print("Addition")
	print("----------------------------")
	print(poly_4 + poly_2)
	print(poly_1 + poly_5)

	print("============================")

	print("Subtraction")
	print("----------------------------")
	print(poly_4 - poly_5)
	print(poly_4 - poly_2)
	print(poly_1 - poly_5)



	print("============================")

	print(poly_4)
	print(poly_5)

	print(poly_4 * poly_5)
	

if __name__=="__main__":
	main()

