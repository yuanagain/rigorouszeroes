from interval import interval, inf, imath, fpu
from complexinterval import ComplexInterval, _one, _zero
from complexpolynomial import ComplexPolynomial

class Newton: 
	def __init__(self, start, poly):
		self.start = start
		self.poly = poly
		self.iterates = 0
		self.deriv = poly.derive()
		self.step = start

	def iterate(self):
		"""
		Performs one Newton iteration, returns change between values.
		"""
		self.iterates += 1
		x = self.step.midpoint()
		fx = self.poly(x)

		## iterate on derivative
		## self.deriv = self.deriv.derive()

		self.step = x - (fx / self.deriv(x))

		## return the change
		diff = x - self.step
		return diff

	def iterate_until(self, res = 10**-6, max_iterates = 20):
		"""
		Iterates until at resolution or until maximum number
		of iterations has been reached. Returns True if convergence
		achieved, returns False otherwise.
		"""
		res_box = ComplexInterval(interval([res, -res]), interval([res, -res]))

		while (self.iterates < max_iterates - 1):
			if self.iterate() in res_box:
				return True

		if self.iterate() in res_box:
				return True

		return False

	def __str__(self):
		"""
		Returns string representation
		"""
		return "Newton's Iterator\n" + "Start: " + str(self.start) + "\nFunction: " + str(self.poly)

def main():
	print("Testing Newton")
	
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

	a_5 = ComplexInterval(interval([2]), interval([0]))
	a_6 = ComplexInterval(interval([2]), interval([0]))

	coeffs = [a_0, a_1, a_2, a_3, a_4, a_5, a_6]

	print("Testing Complex Constructor")
	print("----------------------------")
	poly_1 = ComplexPolynomial(coeffs)
	print(poly_1)
	poly_2 = ComplexPolynomial([_zero(), a_4])
	print(poly_2)
	poly_3 = ComplexPolynomial([a_5, a_6, a_3, a_1, a_0])
	print(poly_3)
	print("============================")


	print("Testing Evaluation")
	print("----------------------------")
	print(poly_1(w))
	print(poly_1(_one()))
	print(poly_1(_zero()))

	print("")

	print(poly_2(w))
	print(poly_2(_one()))
	print(poly_2(_zero()))

	print("")

	print(poly_3(w))
	print(poly_3(_one()))
	print(poly_3(_zero()))
	print("============================")

	print("Derivation")
	print("----------------------------")
	print(poly_1.derive())
	print(poly_1.derive().derive())
	print(poly_1.derive().derive().derive())
	print("")

	print(poly_2.derive())
	print(poly_2.derive().derive())
	print("")

	print(poly_3.derive())
	print(poly_3.derive().derive())
	print("============================")

	print("Newton's Method Constructor")
	print("----------------------------")
	start1 = ComplexInterval(interval([0]), interval([0]))
	start2 = ComplexInterval(interval([1]), interval([1]))
	start3 = ComplexInterval(interval([0]), interval([0]))
	n_1 = Newton(start1, poly_1)
	n_2 = Newton(start2, poly_2)
	n_3 = Newton(start3, poly_3)
	print(n_1)
	print("")
	print(n_2)
	print("")
	print(n_3)
	print("")
	print("============================")

	print("Testing Iteration")
	print("----------------------------")
	for i in range(10):
		print(n_1.iterate())
	print("----------------------------")
	for i in range(10):
		print(n_2.iterate())
	print("----------------------------")
	for i in range(10):
		print(n_3.iterate())
	# print(fpu.isnan(n_3.iterate().a))
	print("============================")

	print("Testing convergence")
	print("----------------------------")
	print(n_1.iterate_until())
	print("----------------------------")
	print(n_2.iterate_until())
	print("----------------------------")
	print(n_3.iterate_until())
	# print(fpu.isnan(n_3.iterate().a))
	print("============================")



if __name__=="__main__":
	main()