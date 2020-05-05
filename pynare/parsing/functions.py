""" 
module for recognizing, parsing, and interpreting the built-in functions 
supported by Pynare 
"""

import numpy as np
from scipy.special import erf 

import pynare.parsing.ast as ast
import pynare.parsing.base as base

import pynare.parsing.algebra as alg

from pynare.parsing.tokens import Token
from pynare.parsing.generic import ABCVisitor


FUNCTION 	= 'FUNCTION'

EXP 	= 'exp'
LOG 	= 'log'
LN		= 'ln'
LOG10 	= 'log10'
SQRT 	= 'sqrt'
ABS 	= 'abs'
SIGN 	= 'sign'
SIN 	= 'sin'
COS 	= 'cos'
TAN 	= 'tan'
ASIN 	= 'asin'
ACOS 	= 'acos'
ATAN 	= 'atan'
MAX 	= 'max'
MIN 	= 'min'
NORMCDF = 'normcdf'
NORMPDF = 'normpdf'
ERF 	= 'erf'


RESERVED_KEYWORDS = {
	EXP: 		Token(FUNCTION, EXP),
	LOG: 		Token(FUNCTION, LOG),
	LN: 		Token(FUNCTION, LN),
	LOG10: 		Token(FUNCTION, LOG10),
	SQRT: 		Token(FUNCTION, SQRT),
	ABS: 		Token(FUNCTION, ABS),
	SIGN: 		Token(FUNCTION, SIGN),
	SIN: 		Token(FUNCTION, SIN),
	COS: 		Token(FUNCTION, COS),
	TAN: 		Token(FUNCTION, TAN),
	ASIN: 		Token(FUNCTION, ASIN),
	ACOS: 		Token(FUNCTION, ACOS),
	ATAN: 		Token(FUNCTION, ATAN),
	MAX: 		Token(FUNCTION, MAX),
	MIN: 		Token(FUNCTION, MIN),
	NORMCDF: 	Token(FUNCTION, NORMCDF),
	NORMPDF: 	Token(FUNCTION, NORMPDF),
	ERF: 		Token(FUNCTION, ERF)
}


class Parser(alg.Parser):
	""" 
	The FunctionParser adds a call() method inbetween the exponent() and atom() 
	methods of the super-class AlgebraParser
	"""

	def __init__(self, lexer):
		super().__init__(lexer)

	def call(self):
		node = self.atom()

		while self.current_token.type == FUNCTION:
			token = self.current_token
			self.eat(FUNCTION)
			self.eat(base.LPARE)

			node = ast.Function(token=token,
								expr=self.expr())
			self.eat(base.RPARE)

		return node

	def exponent(self):
		node = self.call()

		while self.current_token.type == base.POWER:
			token = self.current_token
			self.eat(base.POWER)

			node = ast.BinaryOp(left=node,
								op=token,
								right=self.exponent())
		return node




class Interpreter(alg.Interpreter):

	def visit_Function(self, node):
		f = node.value

		if f == EXP:
			return np.exp(self.visit(node.expr))

		if (f == LOG) or (f == LN):
			return np.log(self.visit(node.expr))

		if f == LOG10:
			return np.log10(self.visit(node.expr))

		if f == SQRT:
			return np.sqrt(self.visit(node.expr))

		if f == ABS:
			return np.abs(self.visit(node.expr))

		if f == SIGN:
			return np.sign(self.visit(node.expr))

		if f == SIN:
			return np.sin(self.visit(node.expr))

		if f == COS:
			return np.cos(self.visit(node.expr))

		if f == TAN:
			return np.tan(self.visit(node.expr))

		if f == ASIN:
			return np.arcsin(self.visit(node.expr))

		if f == ACOS:
			return np.arccos(self.visit(node.expr))

		if f == ATAN:
			return np.arctan(self.visit(node.expr))

		if f == MAX:
			raise NotImplementedError(MAX)

		if f == MIN:
			raise NotImplementedError(MIN)	

		if f == NORMCDF:
			raise NotImplementedError(NORMCDF)

		if f == NORMPDF:
			raise NotImplementedError(NORMPDF)	

		if f == ERF:
			return erf(self.visit(node.expr))
























