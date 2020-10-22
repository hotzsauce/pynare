
from __future__ import annotations

from typing import Union


from pynare.core.generic import ABCModel
from pynare.core.steady import _SteadyStateManager
from pynare.core.deriv.linearizing import compute_jacobian

from pynare.parsing.api import read_model

from pynare.utils.io import (
	read_file,
	read_example
)


class Model(ABCModel):

	def __init__(
		self, 
		outline,
		language: str = 'dynare',
		*args, 
		**kwargs
	):
		super().__init__(outline, language)

	@classmethod
	def from_path(
		cls, 
		filepath: Union[str, Path],
		language: str = 'dynare'
	) -> Model:
		"""

		"""
		try:
			model_definition = read_file(filepath)
		except FileNotFoundError:
			try:
				model_definition = read_example(filepath)
			except ValueError:
				raise FileNotFoundError(filepath)

		model_outline = read_model(model_definition, language)
		return Model(model_outline, language)


	def compute_steady_state(
		self,
		method: Union[str, int, Any] = None,
		**kwargs
	):
		"""
		Computes and saves the steady state of the model. Stores the results of 
		the optimization algorithm in '_steady_state_results', in case objective
		function value, number of iterations, or number of function evaluations
		are desired. The steady state vector is stored as a numpy array in 
		'_steady_state_values'. The steady state vector is also returned

		Parameters
		----------
		method : int | str
			method specification. Whether to use an integer or string depends on 
			what language the model was defined in 
		kwargs : dict
			optional keyword arguments to pass to the steady-state algorithm. Valid
			keywords depend on what algorithm is chosen

		Returns
		-------
		numpy 0-dim array (vector) of steady state values
		"""
		ss_algo = _SteadyStateManager(
			language=self.language,
			method=method
		)

		self._steady_state_results = ss_algo(self, **kwargs)
		self._steady_state_values = self._steady_state_results.x

		self.endogenous.set_steady_state(self._steady_state_values)
		return self._steady_state_values


	def change_parameter(
		self, 
		param_name: str, 
		value: float
	):
		"""
		Change the value of a model parameter

		Parameters
		----------
		param_name : str
			the name of the parameter that will be changed
		value : float
			the new value of the parameter
		"""
		self._parameters[param_name] = value
		self.update_model()


	def update_model(self):
		"""
		After updating the model in some way - be it changing a parameter, changing
		a model expression, or some other alteration - this resets the steady state
		of the model
		"""
		self._init_steady_state_exprs()
		self.compute_steady_state()

	@property
	def jacobian(self):
		# returns the dynamic Jacobian of the model
		if hasattr(self, '_jacobian'):
			return self._jacobian

		else:
			self._jacobian = compute_jacobian(self)
			return self._jacobian

	@property
	def steady_state(self):
		"""
		Access the steady-state values of the model. If necessary, compute them
		"""
		if hasattr(self, '_steady_state_values'):
			pass
		else:
			_ = self.compute_steady_state()
		return self._steady_state_values