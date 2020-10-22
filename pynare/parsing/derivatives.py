
import copy
import numpy as np

import pynare.parsing.ast as ast
import pynare.parsing.base as base


def number(value):
	return ast.Num(base.Token(base.NUMBER, value))

def add(left, right):
	return ast.BinaryOp(left, base.Token(base.PLUS, '+'), right)

def minus(left, right):
	return ast.BinaryOp(left, base.Token(base.MINUS, '-'), right)

def mul(left, right):
	return ast.BinaryOp(left, base.Token(base.MUL, '*'), right)

def div(left, right):
	return ast.BinaryOp(left, base.Token(base.DIV, '/'), right)

def power(left, right):
	return ast.BinaryOp(left, base.Token(base.POWER, '**'), right)


class FunctionDerivative(object):

	def _copy(self, tree):
		return copy.deepcopy(tree)

	def function_derivative(self, node):
		method = 'diff_' + node.value
		fprime = getattr(self, method, self._notimplemented_diff)
		return fprime(node)

	def _notimplemented_diff(self, node):
		msg = 'Derivative for {} not implemented.'
		raise NotImplementedError(msg.format(repr(node.value)))

	def diff_exp(self, node):
		""" e^x -> e^x """
		return self._copy(node)

	def diff_log(self, node):
		""" log(x) -> x^(-1) """
		return power(
			left=self._copy(node.expr),
			right=number(-1)
		)

	def diff_ln(self, node):
		""" log(x) -> x^(-1) """
		return self.diff_log(node)

	def diff_log10(self, node):
		""" log10(x) -> (x*ln(10))^(-1) """
		denom = mul(
			left=self._copy(node.expr),
			right=number(np.log(10))
		)
		return power(
			left=denom,
			right=number(-1)
		)

	def diff_sqrt(self, node):
		""" sqrt(x) -> (2*x^(1/2))^(-1) """
		radical = power(
			left=self._copy(node.expr),
			right=number(-1/2)
		)
		return div(
			left=radical,
			right=number(2)
		)

	def diff_abs(self, node):
		""" abs(x) -> sign(x) """
		return ast.Function(
			token=base.Token(base.FUNCTION, base.SIGN),
			expr=self._copy(node.expr)
		)

	def diff_sign(self, node):
		""" sign(x) -> 0 , ignoring undefined value at 0 """
		return number(0)

	def diff_sin(self, node):
		""" cos(x) -> sin(x) """
		return ast.Function(
			token=base.Token(base.FUNCTION, base.COS),
			expr=self._copy(node.expr)
		)

	def diff_cos(self, node):
		""" cos(x) -> -sin(x) """
		sine = ast.Function(
			token=base.Token(base.FUNCTION, base.IN),
			expr=self._copy(node.expr)
		)
		return ast.UnaryOp(
			op=base.Token(base.MINUS, '-'),
			expr=sine
		)

	def diff_tan(self, node):
		""" tan(x) -> (cos(x))^(-2) """
		cosine = ast.Function(
			token=base.Token(base.FUNCTION, base.COS),
			expr=self._copy(node.expr)
		)
		return power(
			left=cosine,
			right=number(-2)
		)

	def diff_asin(self, node):
		""" asin(x) -> (1 - x^2)^(-1/2) """
		xsq = power(
			left=self._copy(node.expr),
			right=number(2)
		)
		one_minus_xsq = minus(
			left=number(1),
			right=xsq
		)
		return power(
			left=one_minus_xsq,
			right=number(-1/2)
		)

	def diff_acos(self, node):
		""" acos(x) -> -(1 - x^2)^(-1/2) """
		radical = self.diff_asin(node)
		return ast.UnaryOp(
			op=base.Token(base.MINUS, '-'),
			expr=radical
		)

	def diff_atan(self, node):
		""" atan(x) -> (1 + x^2)^(-1) """
		xsq = power(
			left=self._copy(node.expr),
			right=number(2)
		)
		one_plus_xsq = add(
			left=number(1),
			right=xsq
		)
		return power(
			left=one_plus_xsq,
			right=number(-1)
		)

	def diff_floor(self, node):
		""" floor(x) -> 0, ignoring discontinuities where x is integer """
		return number(0)

	def diff_ceil(self, node):
		""" ceil(x) -> 0, ignoring discontinuities where x is integer """
		return number(0)

	def diff_erf(self, node):
		""" erf(x) -> (2/sqrt(pi))*exp(-x^2) """
		sqrtpi = ast.Function(
			token=base.Token(base.FUNCTION, base.SQRT),
			expr=ast.Num(base.Token(base.NUMBER, np.pi))
		)
		twopi = div(
			left=number(2),
			right=sqrtpi
		)

		xsq = power(
			left=self._copy(node.expr),
			right=number(2)
		)
		nxsq = ast.UnaryOp(
			op=base.Token(base.MINUS, '-'),
			expr=xsq
		)
		exp = ast.Function(
			token=base.Token(base.FUNCTION, base.EXP),
			expr=nxsq
		)
		return mul(
			left=twopi,
			right=exp
		)





class NodeFinder(base.ABCVisitor):
	"""
	boolean AST walker that checks if a node type is present in a tree. To be
	specific, for each node class it checks the:
		BinaryOp Node 			- operation
		UnaryOp Node 			- operation
		Num Node 				- value
		Var Node 				- variable name
		Function Node 			- function name
		MultivarFunction Node 	- function name
		Arg Node 				- checkings arg expression
	"""

	def __init__(self, tree):
		super().__init__(tree)
		self.value = None
		self.has_type = False

	def contains(self, value):
		self.has_type = False
		self.value = value
		self.visit(self.tree)
		return self.has_type

	def check_match(self, node_value):
		if node_value == self.value:
			self.has_type = True

	def visit_BinaryOp(self, node):
		self.check_match(node.op.value)
		self.visit(node.left)
		self.visit(node.right)

	def visit_UnaryOp(self, node):
		self.check_match(node.op.value)
		self.visit(node.expr)

	def visit_Num(self, node):
		self.check_match(node.value)

	def visit_Var(self, node):
		self.check_match(node.value)

	def visit_Function(self, node):
		self.check_match(node.value)
		self.visit(node.expr)

	# def visit_MultivarFunction(self, node):
	# 	self.check_match(node.value)
	# 	for arg in node.arguments:
	# 		self.visit(arg)

	# def visit_Arg(self, node):
	# 	self.visit(arg.expr)





def is_constant(tree, wrt):
	"""
	checks that an (full or partial) AST is constant with respect to some variable 
	"""
	checker = NodeFinder(tree)
	return not checker.contains(wrt)





class _Jacobian(base.ABCVisitor, FunctionDerivative):

	def __init__(
		self,
		tree: ast.AST,
		diff_var: str
	):
		self.tree = self._copy(tree)
		self.diff_var = diff_var


	def differentiate(self):
		self.visit(self.tree)
		return self.tree

	def _is_constant(self, node):
		return is_constant(node, self.diff_var)


	def _copy(self, tree):
		return copy.deepcopy(tree)


	def maybe_singleton(self, node, attr):
		"""
		Checks if node's attr attribute if it is a singleton - a Var, Num or 
		Function node. If so, take the derivative of it. If not, walk the 
		attribute as normal
		"""
		attribute = getattr(node, attr)
		if isinstance(attribute, ast.Var):
			if attribute.value == self.diff_var:
				setattr(node, attr, number(1))
			else:
				setattr(node, attr, number(0))
		elif isinstance(attribute, ast.Num):
			setattr(node, attr, number(0))
		elif isinstance(attribute, ast.Function):
			setattr(node, attr, self._chain_rule(attribute))
		elif isinstance(attribute, ast.MultivarFunction):
			raise NotImplementedError
		else:
			self.visit(attribute)


	def _chain_rule(self, node):
		"""
		Use function_derivative method of the FunctionDerivative superclass to 
		compute the derivative of the functions. That method doesn't use the 
		Chain Rule though; i.e. it handles f(g(x)) as if it's f(x). So, we 
		multiply by g'(x) here (if necessary)
		"""
		gprime = self._copy(node)
		self.maybe_singleton(gprime, 'expr')
		fprime = self.function_derivative(node)

		return mul(gprime.expr, fprime)


	def _product_rule(self, node):
		"""
		(f*g)' = f'*g + f*g' 
		"""

		# The pairing of left_copy with maybe_singleton's 'left' argument, and 
		# right_copy with the 'right' argument is not important for the product 
		# rule here. But, because the _quotient_rule method makes use of this
		# method, those pairings can't be switched
		left_copy = self._copy(node)
		right_copy = self._copy(node)

		# making the f'*g term
		node.op = base.Token(base.PLUS, '+')
		self.maybe_singleton(left_copy, 'left')
		node.left = mul(left_copy.left, left_copy.right)

		# making the f*g' term
		self.maybe_singleton(right_copy, 'right')
		node.right = mul(right_copy.left, right_copy.right)


	def _quotient_rule(self, node):
		"""
		(f/g)' = (g*f' - f*g')/(g^2)
		"""

		# use the product rule method to evaluate the numerator, then replace the
		#	plus sign with a minus
		num_copy = self._copy(node)
		self._product_rule(num_copy)
		num_copy.op = base.Token(base.MINUS, '-')

		# squaring the denominator function with _new_binary
		node.right = power(self._copy(node.right), number(2))		
		node.op = base.Token(base.DIV, '/')
		node.left = num_copy


	def _power_rule(self, node):
		"""
		This method covers more than just the basic (x^a)' = a*x^(a-1). There
		are four cases. If:
			Base & Exponent Constant 	- left and right nodes set to 0
			B Constant, E Not 			- (a^f)' = ln(a)*f'*a^f
			B Not, E Constant 			- basic power rule
			B Not, E Not 				- (f^g)' = (f^g)*((g*f')/f + g'*ln(f))

		Here's the algebra for the fourth case. Define y = f^g. Then
				ln(y) = g*ln(f)
			 [ln(y)]' = g * (f'/f) + g'*ln(f)
			 (y'/y)   = g * (f'/f) + g'*ln(f)
			  y' 	  = y * ((g*f')/f + g'*ln(f))
		"""

		base_const = self._is_constant(node.left)
		expo_const = self._is_constant(node.right)

		if base_const and expo_const:
			node.left = number(0)
			node.op = base.Token(base.MUL, '*')
			node.right = number(0)

		elif base_const and (not expo_const):
			# first copy is used to evaluate the ln(a) and f' coefficient terms;
			# 	the second copy is used to track the full remaining a^f term
			coef_copy = self._copy(node)
			full_copy = self._copy(node)

			# evaluating the f' term
			self.maybe_singleton(coef_copy, 'right')

			# the ln(a) coefficient
			node.left = ast.Function(
				token=base.Token(base.FUNCTION, base.LN),
				expr=coef_copy.left
			)
			node.op = base.Token(base.MUL, '*')
			node.right = mul(coef_copy.right, full_copy)

		elif (not base_const) and expo_const:
			node_copy = self._copy(node)

			node.left = node_copy.right
			node.op = base.Token(base.MUL, '*')

			exp_minus_one = minus(node_copy.right, number(1))

			node.left = node_copy.left
			node.op = base.Token(base.POWER, '^')
			node.right = exp_minus_one

		elif (not base_const) and (not expo_const):

			node_copy = self._copy(node)
			prime_copy = self._copy(node)

			self.maybe_singleton(prime_copy, 'left')
			self.maybe_singleton(prime_copy, 'right')

			f, g = node_copy.left, node_copy.right
			fp, gp = prime_copy.left, node_copy.right

			lnf = ast.Function(
				token=base.Token(base.FUNCTION, base.LN),
				expr=self._copy(f)
			)
			gfp = mul(g, fp)
			gfpf = div(gfp, f)
			gplnf = mul(gp, lnf)

			node.left = node_copy
			node.op = base.Token(base.MUL, '*')
			node.right = add(gfpf, gplnf)


	def visit_BinaryOp(self, node):
		"""
		for products, divisions, and exponents, use their respective methods. For
		addition and subtraction, check to see if the 'left' and 'right' nodes
		are singleton variables before visiting their children ASTs
		"""
		if node.op.type == base.PLUS:
			self.maybe_singleton(node, 'left')
			self.maybe_singleton(node, 'right')

		elif node.op.type == base.MINUS:
			self.maybe_singleton(node, 'left')
			self.maybe_singleton(node, 'right')

		elif node.op.type == base.MUL:
			self._product_rule(node)

		elif node.op.type == base.DIV:
			self._quotient_rule(node)

		elif node.op.type == base.POWER:
			self._power_rule(node)


	def visit_UnaryOp(self, node):
		""" Check if singleton variable; visit child nodes if not """
		self.maybe_singleton(node, 'expr')


	def visit_Num(self, node):
		""" 
		if this node is reached it means the entire tree is just the one 
		number node. In larger trees, the maybe_singleton method handles
		numbers within binary & unary operations, and expressions
		"""
		self.tree = number(0)


	def visit_Var(self, node):
		""" 
		if this node is reached it means the entire tree is just the one 
		variable node. In larger trees, the maybe_singleton method handles
		variables within binary & unary operations, and expressions
		"""
		if node.value == self.diff_var:
			self.tree = number(1)
		else:
			self.tree = number(0)


	def visit_Function(self, node):
		""" 
		if this node is reached it means the entire tree is just the one 
		function node. In larger trees, the maybe_singleton method handles
		functions within binary & unary operations, and expressions
		"""
		self.tree = self._chain_rule(node)