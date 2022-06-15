#!/usr/bin/env python
## Codes for Information Manifold.
from sympy import *
## Symbols & variables
x = Symbol('x')
y = Symbol('y')
##################################################################
E = Function('E')(x, y)
F = Function('F')(x, y)
G = Function('G')(x, y)
# Gammaijk := \Gamma^{i}{}_{jk}, for i, j, k \in \{1, 2\}.
Gamma_111 = Function('Gamma_111')(x, y)
Gamma_112 = Function('Gamma_112')(x, y)
Gamma_122 = Function('Gamma_122')(x, y)
Gamma_211 = Function('Gamma_211')(x, y)
Gamma_212 = Function('Gamma_212')(x, y)
Gamma_222 = Function('Gamma_222')(x, y)
##################################################### General nabla
def nabla(VEC1, VEC2, VVEC1, VVEC2):
	'''
	\nabla_{\partial i}{\partial j} = \Gamma^{k}{}_{ij} \partial_k.
	Assume that X = VEC1\partial_1 + VEC2\partial_2, Y = VVEC1\partial_1 + VVEC2\partial_2 
	and f \in C^{\infty}(M). This function, computes (\nabla_{X}Y)(f).
	'''
	_nabla1 = VEC1 * VVEC1 * Gamma_111 + VEC1 * diff(VVEC1, x) + VEC1 * VVEC2 * Gamma_112 + VEC2 * VVEC1 * Gamma_112 + VEC2 * diff(VVEC1, y) + VEC2 * VVEC2 * Gamma_122
	_nabla2 = VEC1 * VVEC1 * Gamma_211 + VEC1 * VVEC2 * Gamma_212 + VEC1 * diff(VVEC2, x) + VEC2 * VVEC1 * Gamma_212 + VEC2 * VVEC2 * Gamma_222 + VEC2 * diff(VVEC2, y)
	_nabla = [_nabla1, _nabla2]
	return _nabla
##################################################### Metric computations
'''
Here, we compute g(X, Y); where X = VEC1\partial_1 + X2\partial_2 and Y = Y1\partial_1 + Y2\partial_2.
'''
def inner_g(VEC1, VEC2, VVEC1, VVEC2):
	return VEC1 * VVEC1 * E + (VEC1 * VVEC2 + VEC2 * VVEC1) * F + VEC2 * VVEC2 * G
#####################################################
##################################################### Action of section on C^{\infty}-function
'''
If X = X1\partial_1 + X2\partial_2 and f \in C^{\infty}(M), then
X.f = X(f) = X^i \frac{\partial f}{\partial x^i}.
'''
def dir_der(VEC1, VEC2, input_func):
	return VEC1 * diff(input_func, x) + VEC2 * diff(input_func, y)
##################################################### Test
#################################################### Metric compability
def nabla_g(X1, X2, Y1, Y2, Z1, Z2):
	##X.g(Y, Z)
	gYZ = inner_g(Y1, Y2, Z1, Z2)
	a = dir_der(X1, X2, gYZ)

	##g(\nabla_X Y, Z)
	nablaXY = nabla(X1, X2, Y1, Y2)
	b = inner_g(nablaXY[0], nablaXY[1], Z1, Z2)

	#g(Y, \nabla_X Z)
	nablaXZ = nabla(X1, X2, Z1, Z2)
	c = inner_g(Y1, Y2, nablaXZ[0], nablaXZ[1])
	
	ex = a - b - c
	return ex
#################################################### Test



## Codes for Levi-Civita connection induced by a Riemannian metric

## Symbols & variables
x = Symbol('x')
y = Symbol('y')
##
#input_func = Function('input_func')(x, y)
#X1 = Function('X1')(x, y)
#X2 = Function('X2')(x, y)
#Y1 = Function('Y1')(x, y)
#Y2 = Function('Y2')(x, y)

##
##################################################################
'''
	These are the components of a Riemannian mtric on \( \mathbb{R}^{3}\) that is induced on a surface.
	E = g_11, F = g_12 = g_21, G = g_22
'''
## 
E_x = diff(E, x)
E_y = diff(E, y)
F_x = diff(F, x)
F_y = diff(F, y)
G_x = diff(G, x)
G_y = diff(G, y)
'''
	These are Christoffel components of the Riemannian metric.
'''
# Main Denominator
MainD = 2 * (E * G - F ** 2)
# Gammaijk := \Gamma^{i}{}_{jk}, for i, j, k \in \{1, 2\}.
Gamma_111 = (G * E_x - 2 * F * F_x + F * E_y) / MainD
Gamma_112 = (G * E_y - F * G_x) / MainD
Gamma_122 = (2 * G * F_y - G * G_x - F * G_y) / MainD
Gamma_211 = (2 * E * F_x - E * E_y - F * E_x) / MainD
Gamma_212 = (E * G_x - F * E_y) / MainD
Gamma_222 = (E * G_y - 2 * F * F_y + F * G_x) / MainD
##
GammaList = [Gamma_111, Gamma_112, Gamma_122, Gamma_211, Gamma_212, Gamma_222]
##################################################### Riemannian tensor and its inverse
g = Matrix(2,2, [E, F, F, G])
g_inverse = g ** (-1)
#####################################################
##################################################### General nabla
def nabla(VEC1, VEC2, VVEC1, VVEC2):
	'''
	\nabla_{\partial i}{\partial j} = \Gamma^{k}{}_{ij} \partial_k.
	Assume that X = VEC1\partial_1 + VEC2\partial_2, Y = VVEC1\partial_1 + VVEC2\partial_2 
	and f \in C^{\infty}(M). This function, computes (\nabla_{X}Y)(f).
	'''
	_nabla1 = VEC1 * VVEC1 * Gamma_111 + VEC1 * diff(VVEC1, x) + VEC1 * VVEC2 * Gamma_112 + VEC2 * VVEC1 * Gamma_112 + VEC2 * diff(VVEC1, y) + VEC2 * VVEC2 * Gamma_122
	_nabla2 = VEC1 * VVEC1 * Gamma_211 + VEC1 * VVEC2 * Gamma_212 + VEC1 * diff(VVEC2, x) + VEC2 * VVEC1 * Gamma_212 + VEC2 * VVEC2 * Gamma_222 + VEC2 * diff(VVEC2, y)
	_nabla = [_nabla1, _nabla2]
	return _nabla
#####################################################
##################################################### Metric computations
'''
Here, we compute g(X, Y); where X = VEC1\partial_1 + X2\partial_2 and Y = Y1\partial_1 + Y2\partial_2.
'''
def inner_g(VEC1, VEC2, VVEC1, VVEC2):
	return VEC1 * VVEC1 * E + (VEC1 * VVEC2 + VEC2 * VVEC1) * F + VEC2 * VVEC2 * G
#####################################################
##################################################### Action of section on C^{\infty}-function
'''
If X = X1\partial_1 + X2\partial_2 and f \in C^{\infty}(M), then
X.f = X(f) = X^i \frac{\partial f}{\partial x^i}.
'''
def dir_der(VEC1, VEC2, input_func):
	return VEC1 * diff(input_func, x) + VEC2 * diff(input_func, y)
#####################################################
##################################################### Curvature
'''
For a 2-dim Riemannian manifold, R_{1212} is enough. 
'''
def R1212():
	return (1 / 2) * ( 2 * diff(F, x, y) - diff(E, y, y) - diff(G, x, x)) + E * (Gamma_112 ** 2 - Gamma_122 * Gamma_111) + F * (2 * Gamma_112 * Gamma_212 - Gamma_122 * Gamma_211 - Gamma_222 * Gamma_111) + G * (Gamma_212 ** 2 - Gamma_222 * Gamma_211)
#R1212().replace(E, x ** 2).replace(F, y).replace(G, x * y).doit()
#################################################### Ricci in local form generally
#################################################### Ricci in dimension 2 in orthogonal frame
def R11():
	return diff(Gamma_211, y) - diff(Gamma_212, x) + Gamma_111 * Gamma_212 + Gamma_211 * Gamma_222 - Gamma_112 * Gamma_211 - Gamma_212 * Gamma_212

def R22():
	return diff(Gamma_122, x) - diff(Gamma_112, y) + Gamma_122 * Gamma_111 + Gamma_222 * Gamma_112 - Gamma_112 * Gamma_112 - Gamma_122 * Gamma_212
def R12():
	return diff(Gamma_112, x) - diff(Gamma_111, y) + Gamma_112 * Gamma_212 + Gamma_222 * Gamma_212 - Gamma_211 * Gamma_122 - Gamma_212 * Gamma_222
def R21():
	return diff(Gamma_212, y) - diff(Gamma_222, x) + Gamma_112 * Gamma_111 + Gamma_112 * Gamma_212 - Gamma_112 * Gamma_111 - Gamma_122 * Gamma_211
#################################################### Scalar curvature
'''
	R_{ij} = R^{1}{}_{i1j} + R^{2}{}_{i2j}
	Trace of the Ricci is the scalar curvature. We shall show it by \( S \).
'''
def sectional_cur():
	return g_inverse[0] * R11() + g_inverse[1] * (R12()  + R21()) + g_inverse[3] * R22()




## Geodesic computations...
## Symbols & variables
x = Symbol('x')
y = Symbol('y')

##  Initial point: p_0 = (p01, p02)
# p01 = 
# p02 = 
## Initial rate: A  = (a^1, a^2)
# a1 = 
# a2 =  
## Functions
'''
	These are the components of a Riemannian mtric on \( \mathbb{R}^{3}\) that is induced on a surface.
'''
	These are Christoffel components of the Riemannian metric.
##################################################### Code Functions
def D2(TwoVarFunc, number):
	
	'''
		This function acts on a two-variable function and a number.
		Then lists all of its partial derivatives until number.
		=================
		Example:
		D2(sin(x * y), 2) 
		=================
	'''
	
	AnswerList = []
	for i in range(0, number):
		a = []
		a.append(i)
		a.append(number - i)
		AnswerList.append(a)
	AnswerList.append([number, 0])
	Diff = []
	DerList = []
	for List in AnswerList:
		SmallDiffer = "TwoVarFunc_X%sY%s"  % (List[0], List[1])
		DerList.append(diff(TwoVarFunc, x, List[0], y, List[1]))
		Diff.append(SmallDiffer)
	#return(Diff , DerList)
	ExactEqual = []
	for i in range(0, len(Diff)):
		S1 = Diff[i]
		S2 = str(DerList[i])
		S = S1 + " = " + S2
		ExactEqual.append(S)
	return(ExactEqual)	
# End of D2 function

#####################################################
#####################################################
##################################################### End of Code Functions
print(D2(x**2 * y**2, 3))





## run a test
n1g12 = nabla_g(1, 0, 1, 0, 0, 1)
n1g21 = nabla_g(1, 0, 0, 1, 1, 0)
Eq1 = n1g12 - n1g21

n1g22 = nabla_g(1, 0, 0, 1, 0, 1)
n2g21 = nabla_g(0, 1, 0, 1, 1, 0)
Eq2 = n1g22 - n2g21

n1g21 = nabla_g(1, 0, 0, 1, 1, 0)
n2g11 = nabla_g(0, 1, 1, 0, 1, 0)
Eq3 = n1g21 - n2g11

n2g11 = nabla_g(0, 1, 1, 0, 1, 0)
n1g12 = nabla_g(1, 0, 1, 0, 0, 1)
Eq4 = n2g11 - n1g12

n2g12 = nabla_g(0, 1, 1, 0, 0, 1)
n1g22 = nabla_g(1, 0, 0, 1, 0, 1) 
Eq5 = n2g12 - n1g22

n2g21 = nabla_g(0, 1, 0, 1, 1, 0)
n2g12 = nabla_g(0, 1, 1, 0, 0, 1)
Eq6 = n2g21 - n2g12
##

print(Eq1)
print('-----')
print(Eq2)
print('-----')
print(Eq3)
print('-----')
print(Eq4)
print('-----')
print(Eq5)
print('-----')
print(Eq6)
print('-----')



x1 = Symbol('x1') # first component of direction vector
y1 = Symbol('y1') # second component of direction vector
G1 = Function('G1')(x, y, x1, y1)
G2 = Function('G2')(x, y, x1, y1)
## Spray
def spray(G1, G2):
	'''
		We define spray by formula 
		S = x1 \par_x + y1 \par_y - 2G1\par_x1 -2G\par_y1,
		as a 4-list. 
	'''
	return [x1, y1, - 2 * G1, - 2 * G2] 

print(spray(2 * x1, y1 - y + x))

