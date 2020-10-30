from __future__ import annotations

from collections import Sequence

import numpy as np 




class ABCModelVariable(object):

	def __init__(self, name: str, vtype: str):
		self.name = name
		self.vtype = vtype


class VariableList(Sequence):

	def __init__(self, vtype: str, variables: Iterable[ABCModelVariable]):
		self.vtype = vtype
		self._variables = variables

	def __getitem__(self, key):
		return self._variables[key]

	def __len__(self):
		return len(self._variables)



# ---------- Endogenous Variables ----------

endo_bounds = ('initial', 'terminal', 'historical')
class EndogenousVariable(ABCModelVariable):

	initial = None
	terminal = None
	historical = None
	steady = None

	def __init__(self, name: str, lead_lag=None, *args, **kwargs):
		super().__init__(name, 'endogenous')
		self.lead_lag = lead_lag

		for k, v in kwargs.items():
			if k in endo_bounds:
				setattr(self, k, v)

	def __repr__(self):
		bounds = {
			b: getattr(self, b) for b in endo_bounds
			if getattr(self, b) is not None
		}
		return f'<EndoVar({self.name}, {bounds}, periods:{self.lead_lag})>'


class EndogenousVariableList(VariableList):

	@classmethod
	def from_outline(self, outline: ModelOutline):
		bounds = {
			'initial': '_initial_values',
			'terminal': '_terminal_values',
			'historical': '_historical_values'
		}
		bound_dicts = {
			k: getattr(outline, v) for k, v in bounds.items() 
			if getattr(outline, v)
		}

		# although we don't know which of the bounds will be given, the lead lags
		# 	dict will always be nonempty
		lead_lags = outline._endo_lead_lags
		var_dict = {}
		for var_name in lead_lags.keys():
			var_dict[var_name] = {b: d[var_name] for b, d in bound_dicts.items()}

		v = [EndogenousVariable(n, ll, **var_dict[n]) for n, ll in lead_lags.items()]

		return EndogenousVariableList('endogenous', v)




# ---------- Stochastic Variables ----------

class StochasticVariable(ABCModelVariable):

	def __init__(self, name: str, variance=None):
		super().__init__(name, 'stochastic')
		self.variance = variance


class StochasticVariableList(VariableList):

	@classmethod
	def from_outline(cls, outline: ModelOutline):
		v = [StochasticVariable(n, var) for n, var in outline._shocks.items()]
		return StochasticVariableList('stochastic', v)




# ---------- Managing Variable Indices ----------

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

	@classmethod
	def from_outline(cls, outline: ModelOutline):
		# kind of unnecesary, but keeps it consistent with *VariableList classes
		return VariableIndexManager(
			outline._endo_lead_lags,
			outline._stochastic_exogenous
		)

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