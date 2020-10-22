"""
module for interpreting the ASTs of the model definition
"""

from __future__ import annotations

import numpy as np

from typing import Iterable, Callable, Union


from pynare.parsing.base import (
	_BaseEvaluator,
	ASTSubstitution
)


class SteadyStateFunction(_BaseEvaluator):
	"""
	the callable that corresponds to the steady-state of individual model
	expressions. When called, it is assumed that the only inputs are endogenous
	variables for the moment
	"""
	def __init__(
		self,
		tree: AST,
		scope: Union[ScopedMemory, dict],
		model_vars: tuple
	):
		simplified_tree = ASTSubstitution(tree, scope)
		super().__init__(simplified_tree, dict())

		self.model_vars = model_vars

	def __call__(self, args):
		self.args = args
		return self.visit(self.tree)

	def visit_Var(self, node):
		try:
			idx = self.model_vars.index(node.value)
			return self.args[idx]
		except ValueError:
			# an exogenous variable
			return 0

	def visit_PeriodVar(self, node):
		try:
			idx = self.model_vars.index(node.value)
			return self.args[idx]
		except ValueError:
			# an exogenous variable
			return 0


def _vectorize_model_functions(
	funcs: Iterable[SteadyStateFunction]
) -> Callable:
	"""
	This functions is originally intended to take the LHS or RHS of the model
	equations in a model definition and return a single function that returns
	the output of those equations in a single vector

	Parameters
	----------
	funcs : an iterable of 'n' functions

	Returns
	-------
	a single function that evaluates to a numpy vector of length n
	"""
	def _vectorized_func(*args):
		return np.array([f(*args) for f in funcs])

	return np.vectorize(
		_vectorized_func,
		signature='(n)->(m)'
	)