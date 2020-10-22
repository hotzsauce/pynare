
from collections import (
	MutableMapping,
	Sequence
)

import numpy as np 


class ABCModelVariable(object):

	def __init__(
		self,
		name: str,
		vtype: str
	):
		self.name = name
		self.vtype = vtype


class EndogenousModelVariable(ABCModelVariable):

	def __init__(
		self,
		name
	):
		super().__init__(name, 'endogenous')
		self.lead_lag = None
		self.initial = self.terminal = self.historical = self.steady = None

	def __repr__(self):
		bounds = dict()
		for b in ('initial', 'terminal', 'historical'):
			v = getattr(self, b)
			if v is not None:
				bounds[b] = v

		return f'<EndoVar({self.name}, {bounds}, periods:{self.lead_lag})>'

	def __str__(self):
		return f'<EndoVar({self.name}, ss={self.steady})>'

	@property
	def is_static(self):
		# appears only at time t
		return self.lead_lag == {0}

	@property
	def is_pforward(self):
		# purely forward variables must appear at t+1 and cannot appear at t-1
		return (1 in self.lead_lag) & (-1 not in self.lead_lag)

	@property
	def is_pbackward(self):
		# purely backward variables must appear at t-1 and cannot appear at t+1
		return (-1 in self.lead_lag) & (1 not in self.lead_lag)

	@property
	def is_mixed(self):
		# appears at times t-1, t, and t+1
		return {-1, 1}.issubset(self.lead_lag)

	@property
	def is_state(self):
		# state variables are purely backward variables and mixed variables
		return self.is_pbackward | self.is_mixed

	def set_lead_lag(self, ll):
		self.lead_lag = ll

	def set_initial_value(self, value):
		self.initial = value

	def set_terminal_value(self, value):
		self.terminal = value

	def set_historical_value(self, value):
		self.historical = value

	def set_steady_state(self, value):
		self.steady = value


class StochasticModelVariable(ABCModelVariable):

	def __init__(
		self,
		name
	):
		super().__init__(name, 'stochastic')
		self.variance = None

	def __repr__(self):
		return f'<StochVar({self.name}, var={self.variance})>'

	def set_variance(self, value):
		self.variance = value


class DeterministicModelVariable(ABCModelVariable):

	def __init__(
		self,
		name
	):
		super().__init__(name, 'deterministic')





"""
The VariableList and _CallableVarList are pretty closely intertwined - when an 
attribute of VariableList is called that's a method name of the component
ABCModelVariables, we add all those methods to a _CallableVarList instance, then
assuming the number of 'args' is the same length as 'methods', evaluate each
method with the single argument

The 'if isinstance(args, dict)' branch in '__call__' ensures that the arguments
are sent to the right method, as self.method and self.variables are always 
guaranteed to be in the same order
"""
class _CallableVarList(object):

	def __init__(self, methods, variables):
		self.methods = methods
		self.variables = variables

	def __call__(self, args):
		if len(args) > 0:
			if isinstance(args, dict):
				var_names = [v.name for v in self.variables]
				return [
					method(args[v]) for method, v in zip(self.methods, var_names)
				]

			return [method(a) for method, a in zip(self.methods, args)]


class VariableList(Sequence):

	def __init__(
		self,
		variables: list,
		vtype: str
	):
		if vtype == 'endogenous':
			vclass = EndogenousModelVariable
		elif vtype == 'stochastic':
			vclass = StochasticModelVariable
		elif vtype == 'deterministic':
			vclass = DeterministicModelVariable

		self._variables = [vclass(n) for n in variables]
		self._vtype = vtype

	def __getattribute__(self, attr):
		return super(VariableList, self).__getattribute__(attr)

	def __getattr__(self, attr):
		is_callable = callable(getattr(self._variables[0], attr))
		if is_callable:
			variable_methods = [getattr(v, attr) for v in self._variables]
			return _CallableVarList(variable_methods, self._variables)

		return [getattr(v, attr) for v in self._variables]

	def __getitem__(self, key):
		return self._variables[key]

	def __len__(self):
		return len(self._variables)

	def __str__(self):
		return f'<VariableList{self._variables}>'






class ModelParameter(ABCModelVariable):

	def __init__(
		self,
		name,
		value
	):
		super().__init__(name, 'parameter')
		self.value = value

	def __repr__(self):
		return f'<ModelParameter({self.name}: {self.value})>'


class ParameterDict(MutableMapping):

	def __init__(
		self,
		params: dict
	):
		self._params = {n: ModelParameter(n, v) for n, v in params.items()}

	def __getitem__(self, key):
		return self._params[key].value

	def __setitem__(self, key, value):
		mp = ModelParameter(key, value)
		self._params[key] = mp

	def __delitem__(self, key):
		del self._params[key]

	def __iter__(self):
		for p in self._params.keys():
			yield p
		return

	def __len__(self):
		return len(self._params)

	def __repr__(self):
		return f'<ParameterDict{list(self._params.values())}>'




class VariableIndexManager(object):
	"""
	manager class for keeping track of indices of model variables in the Jacobian
	and Decision Rule spaces
	"""


	def __init__(
		self,
		lead_lags: dict,
		stoch_shocks: list
	):	

		"""
		The values of the `lead_lags' dict are sets with elements -1, 0, 1, 
		depending on which of the t-1, t, t+1 time periods the endogenous variable
		appears in. For consistency across the variables, we cast these to a
		sorted tuple
		"""
		self.lead_lags = {k: tuple(sorted(v)) for k, v in lead_lags.items()}

		"""
		Initialize the lead lag incidence matrix. There will always be three 
		rows: the first corresponds to time t-1 appearances; the second, time
		t appearances; the third, time t+1 appearances. Each column corresponds
		to the endogenous variables (which by construction of the ModelOutline's
		_endo_lead_lags attribute, are in declaration order). The values in the 
		matrix give that variable's column index in the Model's jacobian
		"""
		llx = np.empty((3, len(lead_lags)))
		llx[:] = np.nan

		for endo_idx, ll in enumerate(lead_lags.values()):
			idx = np.array(list(ll)) + 1
			llx[idx, endo_idx] = 1

		# dictionary of endo_names: col idx of lead-lag incidence
		self.var_cols = dict(zip(lead_lags.keys(), range(len(lead_lags))))

		# number of endogenous variables at all time periods
		self.n_endo_cols = int(np.nansum(llx))

		# filling in the spots top-to-bottom left-to-right
		llx[llx == 1] = np.arange(self.n_endo_cols)
		self.llx = llx

		self.stoch_shocks = stoch_shocks
		self.n_exos = len(stoch_shocks)


	def jacobian_column(
		self, 
		key: str
	):
		"""
		Given the variable name, returns the column indices of the Jacobian
		corresponding to that variable. 

		Parameters
		----------
		key : str
			variable name

		Returns
		-------
		list of column indices
		"""
		try:
			# endogenous variables
			cols = self.llx[:, self.var_cols[key]]
			return cols[~np.isnan(cols)].astype(int).tolist()
			# return self.llx[:, self.var_cols[key]]

		except KeyError:
			# exogenous variables
			return [self.n_endo_cols + self.stoch_shocks.index(key)]


	def lead_lag(
		self,
		key: str
	):
		"""
		Given the variable name, returns the t-1, t, t+1 time periods
		that variable appears in 

		Parameters
		----------
		key : str
			variable name

		Returns
		-------
		list of -1, 0, 1
		"""	

		# default value of {0} for exogenous varaibles
		return self.lead_lags.get(key, {0})

	@property
	def n_static(self):
		return len(self.static_jacobian)

	@property
	def n_dynamic(self):
		return len(self.dynamic_jacobian)

	@property
	def n_mixed(self):
		llx = self.llx
		_mj = llx[0, ~(np.isnan(llx[0, :]) | np.isnan(llx[2, :]))]
		return len(_mj[~np.isnan(_mj)])
	
	@property 
	def static_jacobian(self):
		# static variables in jacobian matrix are those in the lead lag incidence
		#	matrix that have non-nan values only in the middle row
		llx = self.llx
		_sj = llx[1, np.isnan(llx[0, :]) & np.isnan(llx[2, :])]
		return _sj[~np.isnan(_sj)].astype(int)

	@property
	def dynamic_jacobian(self):
		# dynamic variables in the jacobian matrix are the time t variables that
		#	also appear in periods t-1 or t+1
		llx = self.llx
		_dj = llx[1, ~(np.isnan(llx[0, :]) & np.isnan(llx[2, :]))]
		return _dj[~np.isnan(_dj)].astype(int)

	@property
	def forward_jacobian(self):
		# forward variables in jacobian matrix appear at time t+1
		return self.llx[2, ~np.isnan(self.llx[2, :])].astype(int)

	@property
	def backward_jacobian(self):
		# backward variables in jacobian matrix appear at time t-1
		return self.llx[0, ~np.isnan(self.llx[0, :])].astype(int)

	@property
	def exogenous_jacobian(self):
		# exogenous variables are the right-most columns of the Jacobian
		return np.arange(self.n_endo_cols, self.n_endo_cols+self.n_exos, dtype=int)
	
	@property
	def cont_forward_jacobian(self):
		# contemporaneous variables that also appear at time t+1
		llx = self.llx
		_fj = llx[1, ~np.isnan(llx[2, :])]
		return _fj[~np.isnan(_fj)].astype(int)

	@property
	def cont_pure_forward_jacobian(self):
		# contemporaneous instance of variables that appear at time t+1, but
		#	do not appear at time t-1
		llx = self.llx
		_fj = llx[1, np.isnan(llx[0, :]) & ~np.isnan(llx[2, :])]
		return _fj[~np.isnan(_fj)].astype(int)

	@property
	def cont_backward_jacobian(self):
		# contemporaneous variables that also appear at time t-1
		llx = self.llx
		_fj = llx[1, ~np.isnan(llx[0, :])]
		return _fj[~np.isnan(_fj)].astype(int)
	
	@property
	def pure_backward_jacobian(self):
		# purely backward variables only appear at time t-1
		llx = self.llx
		_bj = llx[0, np.isnan(llx[1, :]) & np.isnan(llx[2, :])]
		return _bj[~np.isnan(_bj)].astype(int)
	
	


















