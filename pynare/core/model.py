
from __future__ import annotations

from pynare.core.generic import ABCModel
from pynare.core.steady import _SteadyStateManager
from pynare.core.deriv.linearizing import compute_jacobian

from pynare.parsing import read_model

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
		Create a model from a pathname. If file is not found, assume it's a
		pre-made example model in pynare/examples

		Parameters
		----------
		filepath : str | Path
			path to .mod or .txt file that defines the model
		language : str
			the language the model is written in

		Returns
		-------
		Model
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

		try:
			self._steady_state_results = ss_algo(self, **kwargs)
		except AttributeError:
			self.update_model()
			self._steady_state_results = ss_algo(self, **kwargs)


		self._steady_state_values = self._steady_state_results.x
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
		self.is_altered = True
		self._parameters[param_name] = value


	def update_model(self):
		"""
		After updating the model in some way - be it changing a parameter, changing
		a model expression, or some other alteration - this resets the steady state
		of the model
		"""
		self.is_altered = False
		self._init_steady_state_funcs()
		self._init_dynamic_funcs()

	@property
	def steady_state_model(self):
		if self.is_altered:
			self.update_model()
		return self._steady_state_funcs

	@property
	def dynamic_model(self):
		if self.is_altered:
			self.update_model()
		return self._dynamic_funcs

	@property
	def jacobian(self):
		# returns the dynamic Jacobian of the model
		if hasattr(self, '_jacobian'):
			if self.is_altered:
				self.update_model()
			self._jacobian = compute_jacobian(self)

		else:
			self._jacobian = compute_jacobian(self)

		return self._jacobian

	@property
	def steady_state(self):
		"""
		Access the steady-state values of the model. If necessary, compute them
		"""
		if hasattr(self, '_steady_state_values'):
			if self.is_altered:
				self.update_model()
				self.compute_steady_state()
		else:
			_ = self.compute_steady_state()

		return self._steady_state_values