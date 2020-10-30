"""
module for translating the ASTs of the model definition into the steady-state
and dynamic functions represented by each model expression
"""

from __future__ import annotations

import inspect 
import functools
import numpy as np

import pynare.parsing.ast as ast
import pynare.parsing.base as base




# ---------- Steady-State Functions ----------

def steady_func(
	node: AST, 
	endo: tuple, 
	params: dict
):
	"""
	Creates & returns a steady-state function based on a single model 
	expression in the definition of a model. 

	Parameters
	----------
	node : AST
		the abstract syntax tree of the model expression
	endo : tuple
		a tuple of the endogenous model names in declaration order. these 
		are the steady state names (no 'm1' or 'p1' PeriodVar suffixes)
	params : dict
		dictionary of parameter_names: parameter_values. This dictionary
		must provide the locally-defined model variables, too

	Returns
	-------
	steady-state function that corresponds to the model expression. The 
	parameters of the function are the variable names in `endo`, with 
	order preserved. The `__name__` attribute is set to 'steady_mexpr',
	and an `__endo__` attribute is a set of the endogenous variables
	that appear in the function
	"""
	sub = _SteadyStateSubstitution(node, params)
	sub_tree, expr_endo = sub.tree, sub.not_provided

	ast_walker = SteadyFunc(sub_tree)

	@model_signature(endo)
	def _steady_func(ss):
		ss_vals = {e: ss[i] for i, e in enumerate(endo) if e in expr_endo}
		return ast_walker(ss_vals)

	_steady_func.__name__ = 'steady_mexpr'
	_steady_func.__endo__ = expr_endo

	return _steady_func


class SteadyFunc(base._BaseEvaluator):
	"""
	the callable that corresponds to the steady-state of individual model
	expressions. When called, it is assumed that the only inputs are endogenous
	variables for the moment
	"""
	def __init__(
		self,
		tree: AST,
	):
		super().__init__(tree, dict())

	def __call__(self, ss_vals):
		self.ss_vals = ss_vals
		return self.visit(self.tree)

	def visit_Var(self, node):
		try:
			return self.ss_vals[node.value]
		except KeyError:
			# an exogenous variable
			return 0

	def visit_PeriodVar(self, node):
		try:
			return self.ss_vals[node.value]
		except KeyError:
			# an exogenous variable
			return 0




# ---------- Dynamic Functions ----------
def dynamic_func(
	node: AST,
	endo: tuple,
	exo: tuple, 
	llx: np.ndarray,
	params: dict
):
	"""
	Creates & returns a dynamic function based on a single model expression 
	in the definition of a model. 

	Parameters
	----------
	node : AST
		the abstract syntax tree of the model expression
	endo : tuple
		a tuple of the endogenous model names in declaration order
	exo : tuple
		a tuple of the exogenous model names in declaration order
	llx : np.ndarray
		the lead-lag incidence matrix of the model
	params : dict
		dictionary of parameter_names: parameter_values. This dictionary
		must provide the locally-defined model variables, too

	Returns
	-------
	steady-state function that corresponds to the model expression. The 
	parameters of the function are the variable names in `endo`, with 
	order preserved. The `__name__` attribute is set to 'steady_mexpr',
	and an `__endo__` attribute is a set of the endogenous variables
	that appear in the function
	"""
	sub = _DynamicSubsitution(node, params)
	sub_tree, expr_endo = sub.tree, sub.not_provided

	# `llx` has `len(endo)` columns by construction, so we make use of that
	#	to get the `endos` names in the correct order
	endo_suffix = np.array(
		[[f'{v}{s}' for v in endo] for s in ('m1', '', 'p1')], dtype=object
	)
	dyn_vars = list(endo_suffix[~np.isnan(llx)])
	dyn_vars.extend(exo) # `exo` in declaration order by construction
	
	ast_walker = DynamicFunc(sub_tree)

	@model_signature(dyn_vars)
	def _dynamic_func(dyn):
		dyn_vals = {e: dyn[i] for i, e in enumerate(dyn_vars) if e in expr_endo}
		return ast_walker(dyn_vals)

	_dynamic_func.__name__ = 'dynamic_mexpr'
	_dynamic_func.__endo__ = expr_endo

	return _dynamic_func


class DynamicFunc(base._BaseEvaluator):
	"""
	the callable that corresponds to the dynamic version of individual model
	expressions
	"""
	def __init__(
		self,
		tree: AST,
	):
		super().__init__(tree, dict())

	def __call__(self, dyn_vals):
		self.dyn_vals = dyn_vals
		return self.visit(self.tree)

	def visit_Var(self, node):
		return self.dyn_vals[node.value]
		
	def visit_PeriodVar(self, node):
		var_name = node.value
		offset = node.period_offset
		name_offset = f'{var_name}p1' if offset > 0 else f'{var_name}m1'
		return self.dyn_vals[name_offset]
		


# ---------- Miscellaneous Functions ----------


class _SteadyStateSubstitution(base._ASTSubstitution):
	"""
	Modify the maybe_substitute of the base class to treat PeriodVar and 
	Var ASTs identically. Also keeps a record of variables that are not
	in the provided `scope` dictionary, which are the endogenous variables
	in the `tree` parameter
	"""	

	def __init__(
		self, 
		tree: AST,
		scope: dict
	):
		super().__init__(tree, scope) # initializes self.tree attribute too
		self.not_provided = set()

		self.visit(self.tree)

	def maybe_substitute(self, node, attr):
		attribute = getattr(node, attr)
		if isinstance(attribute, (ast.Var, ast.PeriodVar)):
			var_name = attribute.value
			try:
				var_value = self.scope[var_name]
				setattr(node, attr, self.substitute(attribute, var_value))
			except KeyError:
				self.not_provided.add(var_name)


class _DynamicSubsitution(base._ASTSubstitution):
	"""
	Modify the maybe_substitute of the base class to differentiate between 
	PeriodVar and Var ASTs. Also keeps a record of variables that are not
	in the provided `scope` dictionary, which are the endogenous variables
	in the `tree`. Var variables are recorded with no name modifications,
	so variable `x` will be added to `self.not_provided` as `x`. However,
	PeriodVars are altered in accordance with the naming scheme used in 
	dynamic_func:
		PeriodVar(name='x', period_offset=-1) -> 'xm1'
		PeriodVar(name='y', period_offset=1)  -> 'ym1'
	"""	

	def __init__(
		self, 
		tree: AST,
		scope: dict
	):
		super().__init__(tree, scope) # initializes self.tree attribute too
		self.not_provided = set()

		self.visit(self.tree)

	def maybe_substitute(self, node, attr):
		attribute = getattr(node, attr)
		if isinstance(attribute, ast.Var):
			var_name = attribute.value
			try:
				var_value = self.scope[var_name]
				setattr(node, attr, self.substitute(attribute, var_value))
			except KeyError:
				self.not_provided.add(var_name)

		elif isinstance(attribute, ast.PeriodVar):
			var_name = attribute.value
			offset = attribute.period_offset
			name_offset = f'{var_name}p1' if offset > 0 else f'{var_name}m1'
			self.not_provided.add(name_offset)



def model_signature(
	mvars: Iterable[str]
):
	"""
	Decorator that redefines a steady-state (or dynamic) model expression's 
	signature into one made up only of positional steady-state (or dynamic)
	model variable names. Although the existing SteadyStateFunction callable
	class (commit 60e5bc594b6a08a5bbff61e5baa320ba3378af06) can figure out 
	which steady-state values are assigned to each variable based just on 
	an `*args` argument in the `__call__` method, the `autograd` library needs
	each endogenous variable to be its own positional argument in order for
	its differential operators to work properly.

	Parameters
	----------
	mvars : Iterable[str]
		iterable of model variable names

	Returns
	-------
	model expression function with `mvars` as its signature

	Based on:
	https://stackoverflow.com/questions/1409295/set-function-signature-in-python
	"""
	def sig_decorator(func):
		@functools.wraps(func)
		def sig_wrapper(*args, **kwargs):
			return func(*args, **kwargs)

		sig_wrapper.__signature__ = inspect.Signature(ensure_parameter(mvars))
		return sig_wrapper

	return sig_decorator


def ensure_parameter(
	names: Iterable[str], 
	kind: str = 'POSITIONAL_ONLY'
) -> Tuple[inspect.Parameter]:
	"""
	Translates an iterable of strings into a tuple of inspect.Parameter objects
	of type `kind`

	Parameters
	----------
	names : Iterable
		an iterable of strings representing the parameter names
	kind : str
		the kind of Parameter to create. Acceptable values are:
			POSITIONAL_ONLY
			POSITIONAL_OR_KEYWORD
			VAR_POSITIONAL
			KEYWORD_ONLY
			VAR_KEYWORD

	Returns
	-------
	tuple of inspect.Parameter objects
	"""
	# translates an iterable of strings into a tuple of inspect.Parameter
	#	of type `kind`
	param_kind = getattr(inspect.Parameter, kind)
	return tuple([inspect.Parameter(n, kind=param_kind) for n in names])


def vectorize_functions(
	funcs: Iterable[Callable]
) -> Callable:
	"""
	Vectorizes an iterable of functions and returns a single function that 
	returns the output of those equations in a single vector

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