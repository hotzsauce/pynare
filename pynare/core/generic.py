from __future__ import annotations

from pynare.core.expressions import steady_func, dynamic_func

from pynare.core.variables import (
	EndogenousVariableList,
	StochasticVariableList,
	VariableIndexManager
)

class ABCModel(object):

	def __init__(
		self,
		outline: ModelOutline,
		language: str,
		*args,
		**kwargs
	):
		self.outline = outline
		self.language = language
		self.is_altered = False

		# using a copy of the outline to set the various attributes allows the 
		#	original model outline to be preserved
		_copy = outline.__deepcopy__()

		# parameters & locally declared variables are stored in dicts of the 
		#	form 'var_name': 'var_value'
		self._parameters = _copy._parameters
		self._local_model_variables = _copy._local_model_variables

		# endogenous and exogenous variables are stored as tuples of variable
		#	names in the order they were declared
		self._endogenous = EndogenousVariableList.from_outline(_copy)
		self._stoch_shocks = StochasticVariableList.from_outline(_copy)
		# self._determ_shocks = _copy._deterministic_exogenous

		# manages indices of variables in declaration & decision-rule order
		self._index_mgr = VariableIndexManager.from_outline(_copy)

		self._model_exprs_asts = _copy._model_expression_asts


	def _init_steady_state_funcs(self):
		"""
		Define the steady state functions based on the ASTs in the model's 
		definition
		"""
		self._steady_state_funcs = list()

		endo = self.endogenous
		scope = self.model_scope

		for n_expr, ast in enumerate(self._model_exprs_asts):
			# define and rename function to 'steady_mexpr{n_expr}'
			ss_func = steady_func(ast, endo, scope)
			ss_func.__name__ = f'{ss_func.__name__}{n_expr}'
			
			self._steady_state_funcs.append(ss_func)

	def _init_dynamic_funcs(self):
		"""
		Define the dynamic functions based on the ASTs in the model's definition
		"""
		self._dynamic_funcs = list()

		endo = self.endogenous
		exo = self.exogenous
		llx = self._index_mgr.llx
		scope = self.model_scope

		for n_expr, ast in enumerate(self._model_exprs_asts):
			# define and rename function to 'dynamic_mexpr{n_expr}'
			dyn_func = dynamic_func(ast, endo, exo, llx, scope)
			dyn_func.__name__ = f'{dyn_func.__name__}{n_expr}'

			self._dynamic_funcs.append(dyn_func)


	@property
	def parameters(self):
		return self._parameters

	@property
	def model_scope(self):
		return {**self._parameters, **self._local_model_variables}

	@property
	def endogenous(self):
		"""
		when accessing public-facing endogenous attribute, just return tuple 
		of the variable names
		"""
		return tuple((v.name for v in self._endogenous))

	@property
	def exogenous(self):
		"""
		when accessing public-facing exogenous attribute, just return tuple 
		of the variable names
		""" 
		return tuple((v.name for v in self._stoch_shocks))