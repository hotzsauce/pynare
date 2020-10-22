"""
module for computing Jacobian & higher-order derivative matrices of the model

For now we compute the derivative by walking the ASTs of the model expressions.
This is pretty slow (the BasicNK example model has 9 expressions and 15 
endogenous & exogenous variables, and takes around 0.4 seconds to compute the
Jacobian), and I have worries about the curse of dimensionality. Automatic
differentiation might be a better solution, and the `autograd' library has 
potential (the associated JAX library looks like overkill). autograd isn't 
actively being developed, however, although I don't know how big of a 
problem that is 
"""

from __future__ import annotations

import numpy as np 

import pynare.parsing.ast as ast
import pynare.parsing.base as base

from pynare.parsing.derivatives import _Jacobian


def compute_jacobian(model: ABCModel):

	_ml = _ModelLinearizer(model)
	return _ml.jac



class _ModelLinearizer(object):

	def __init__(
		self,
		model: ABCModel
	):
		# the scope (model parameters & locally-declared variables) will be
		#	substituted before Jacobian is computed
		self.scope = model.model_scope

		# the steady-state of endogenous variables & exogenous variables 
		#	correspond to the columns of the Jacobian
		self.endo = dict(zip(model.endogenous.name, model.steady_state))
		self.exo = dict(zip(model.shocks.name, np.zeros(len(model.shocks))))
		self.steady = {**self.endo, **self.exo}

		self.idx_mgr = model._index_mgr
		self.model_exprs = model._model_exprs_asts
		self.mexpr_vars = model._model_exprs_vars

		n_exprs = len(self.model_exprs)
		self.jac = np.zeros((n_exprs, self.idx_mgr.n_endo_cols + self.idx_mgr.n_exos))

		for expr_idx in range(n_exprs):
			self._linearize_model_expression(expr_idx)

		
	def _linearize_model_expression(
		self, 
		expr_idx: int
	):
		# get the endogenous variables that are in this model expression. Each
		#	element of the set 'var' is a tuple of the form of (var_name, offset)
		#	where offset is either -1, 0, 1
		model_var = self.mexpr_vars[expr_idx]

		# substitute all the parameter values into the model expression ast
		model_ast = self.model_exprs[expr_idx]
		model_ast = base.ASTSubstitution(model_ast, self.scope)

		for jac_var in self.steady.keys():

			# column of variable in the model jacobian, and the time periods
			#	it appears in 
			jac_cols = self.idx_mgr.jacobian_column(jac_var)
			offset_periods = self.idx_mgr.lead_lag(jac_var)

			for col_idx, offset in zip(jac_cols, offset_periods):

				# we only want to compute the Jacobian for variables that appear
				#	in this model expression
				if (jac_var, offset) in model_var or jac_var in self.exo:

					model_jac = _DynamicJacobian(
						tree=model_ast,
						diff_var=jac_var,
						offset=offset
					)

					j = JacobianEvaluator(model_jac.differentiate(), self.steady)
					self.jac[expr_idx, col_idx] = j



class _DynamicJacobian(_Jacobian):
	"""
	calculates the derivative of an AST with respect to one variable. Basically
	the same as _Jacobian, but we need to add functionality for PeriodVar nodes
	"""

	def __init__(
		self, 
		tree: AST,
		diff_var: str,
		offset: int
	):
		super().__init__(tree, diff_var)

		self.offset = offset


	def maybe_singleton(self, node, attr):
		"""
		Checks if node's attr attribute if it is a singleton - a Var, Num or 
		Function node. If so, take the derivative of it. If not, walk the 
		attribute as normal
		"""
		attribute = getattr(node, attr)

		if isinstance(attribute, ast.Var):
			var_name = attribute.value 
			if (var_name == self.diff_var) & (self.offset == 0):
				setattr(node, attr, ast.Num(base.Token(base.NUMBER, 1)))
			else:
				setattr(node, attr, ast.Num(base.Token(base.NUMBER, 0)))

		elif isinstance(attribute, ast.PeriodVar):
			var_name = attribute.value 
			offset = attribute.period_offset
			if (var_name == self.diff_var) & (self.offset == offset):
				setattr(node, attr, ast.Num(base.Token(base.NUMBER, 1)))
			else:
				setattr(node, attr, ast.Num(base.Token(base.NUMBER, 0)))

		elif isinstance(attribute, ast.Num):
			setattr(node, attr, ast.Num(base.Token(base.NUMBER, 0)))

		elif isinstance(attribute, ast.Function):
			setattr(node, attr, self._chain_rule(attribute))

		else:
			self.visit(attribute)


	def visit_Var(self, node):
		""" 
		if this node is reached it means the entire tree is just the one 
		variable node. In larger trees, the maybe_singleton method handles
		variables within binary & unary operations, and expressions
		"""
		if (node.value == self.diff_var) & (self.offset == 0):
			self.tree = ast.Num(base.Token(base.NUMBER, 1))
		else:
			self.tree = ast.Num(base.Token(base.NUMBER, 0))


	def visit_PeriodVar(self, node):
		""" 
		if this node is reached it means the entire tree is just the one 
		period variable node. In larger trees, the maybe_singleton method handles
		variables within binary & unary operations, and expressions
		"""
		if (node.value == self.diff_var) & (self.offset == node.period_offset):
			self.tree = ast.Num(base.Token(base.NUMBER, 1))
		else:
			self.tree = ast.Num(base.Token(base.NUMBER, 0))



class JacobianEvaluator(object):

	def __new__(
		cls,
		tree,
		scope
	):
		evaluator = _JacobianEvaluator(tree, scope)
		return evaluator.evaluate()



class _JacobianEvaluator(base._BaseEvaluator):
	"""
	adds a visit_PeriodVar method to the _BaseEvaluator class
	"""

	def visit_PeriodVar(self, node):
		try:
			try:
				return self.scope.lookup(node.value)
			except AttributeError:
				return self.scope.get(node.value)

		except KeyError:
			raise SyntaxError(f'symbol {node.value} was not assigned value before use')