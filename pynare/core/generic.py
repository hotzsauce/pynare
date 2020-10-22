from __future__ import annotations

from typing import Union, Any

import numpy as np 

from pynare.core.expressions import (
	_vectorize_model_functions,
	SteadyStateFunction
)

from pynare.core.variables import (
	ParameterDict,
	VariableList,
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

		# using a copy of the outline to set the various attributes allows the 
		#	original model outline to be preserved
		_copy = outline.__deepcopy__()

		self._parameters = ParameterDict(_copy._parameters)
		self._endogenous = VariableList(_copy._endogenous, 'endogenous')
		
		self._endogenous.set_initial_value(_copy._initial_values)
		self._endogenous.set_terminal_value(_copy._terminal_values)
		self._endogenous.set_historical_value(_copy._historical_values)
		self._endogenous.set_lead_lag(_copy._endo_lead_lags)

		self._index_mgr = VariableIndexManager(
			_copy._endo_lead_lags,
			_copy._stochastic_exogenous
		)


		# at the moment the ModelFactory doesn't handle deterministic shocks, and
		#	it will only handle 'var' and 'stderr' stochastic shocks, which is
		#	stored in the '_shocks' dict as the variance regardless of how its
		#	declared
		self._stochastic_shocks = \
			VariableList(_copy._stochastic_exogenous, 'stochastic')
		# self._deterministic_shocks = \
		# 	VariableList(_copy._deterministic_exogenous, 'deterministic')

		self._stochastic_shocks.set_variance(_copy._shocks)



		self._local_model_variables = ParameterDict(_copy._local_model_variables)

		# we save the original ASTs in case a parameter or model expression is 
		#	updated
		self._model_exprs_asts = _copy._model_expression_asts
		self._model_exprs_vars = _copy._model_expression_vars

		self._init_steady_state_exprs()


		
	def _init_steady_state_exprs(self):
		"""
		initialize the steady state expressions of the model. 
		"""
		self._steady_state_exprs = _vectorize_model_functions(
			[
				SteadyStateFunction(
					tree=ast,
					scope=self.model_scope,
					model_vars=self.outline._endogenous
				)
				for ast in self._model_exprs_asts
			]
		)

	@property
	def state_variables(self):
		return [v for v in self._endogenous if v.is_state]
	
	@property
	def model_scope(self):
		return {**self._parameters, **self._local_model_variables}
	
	@property
	def parameters(self):
		return self._parameters

	@property
	def parameter_values(self):
		return self._parameters

	@property
	def endogenous(self):
		return self._endogenous

	@property
	def shocks(self):
		return self._stochastic_shocks # + self._deterministic_shocks