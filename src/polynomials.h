#ifndef POLYNOMIAL_HEADER
#define POLYNOMIAL_HEADER

int solve_cubic_equation(double c3, double c2, double c1, double c0,
			 double roots[3]);

int solve_quartic_equation(double d4, double d3,
			   double d2, double d1, double d0,
			   double roots[4]);

#endif /* POLYNOMIAL_HEADER */
