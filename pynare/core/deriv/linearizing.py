"""
module for computing Jacobian & higher-order derivative matrices of the model
"""

from __future__ import annotations

import numpy as np 

from pynare.core.deriv.autograd_wrapper import grad_array


def compute_jacobian(
	model: Model
):
	
	# `llx` has `len(steady_state)` columns by construction, and each column of 
	#	the two arrays correspond to the same endo variable. We make use of that
	#	when constructing the steady-state vector with variables ordered in the 
	#	dynamic (t-1, t, t+1) declaration order
	llx = model._index_mgr.llx
	ss_block = np.tile(model.steady_state, (3, 1))
	ss_endos = ss_block[~np.isnan(llx)].flatten()

	# all shocks are zero in steady state
	zero_shocks = np.zeros(len(model.exogenous))
	ss = np.concatenate((ss_endos, zero_shocks))

	# list of the functions representing the dynamic model
	dyn_funcs = model.dynamic_model

	n_exprs, n_cols = len(dyn_funcs), len(ss)
	jac = np.zeros((n_exprs, n_cols))
	
	# grad_array wraps autograd.grad & returns an array instead of tup-of-tups
	for i, func in enumerate(dyn_funcs):
		grad = grad_array(func)
		jac[i, :] = grad(ss)
	
	return jac