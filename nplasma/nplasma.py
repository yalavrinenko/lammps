import math
class Unit:
	def __init__(self, si=0.0, cgs=0.0):
		self.si = si
		self.cgs = cgs

class Contants:
	kB = Unit(1.3807e-23, 1.3807e-16)
	e = Unit(1.6022e-19, 4.8032e-10)
	me = Unit(9.1094e-31, 9.1094e-28)
	mp = Unit(1.6726e-27, 1.6726e-24)

	nu0 = Unit(4.0 * 3.14159265 * 1e-7, 1)
	eps0 = Unit(8.8542e-12, 1)


def T_eV(Tk):
	kT = Contants.kB.si * Tk / Contants.e.si
	return kT


def T_K(TeV):
	return TeV * Contants.e.si / Contants.kB.si


def w_pl(ne):
	cgs = (4.0 * math.pi * ne * Contants.e.cgs**2 / Contants.me.cgs)**0.5
	return cgs


def Gamma(ne, Tk):
	c1 = (4.0 * math.pi / 3.0) ** (1.0 / 3.0)
	c2 = Contants.e.cgs**2 * ne ** (1.0 / 3.0)

	return c1 * c2 / (Contants.kB.cgs * Tk)


def G_el_density(G, Tk):
	c1 = (4.0 * math.pi / 3.0) ** (1.0 / 3.0)
	c2 = G * (Contants.kB.cgs * Tk) / c1
	c3 = c2 / (Contants.e.cgs**2)
	return c3**3.0
