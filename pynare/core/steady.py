"""
module with functions for finding the steady state of a model
"""

from __future__ import annotations

import numpy as np
import scipy.optimize as optim

from pynare.core.expressions import vectorize_functions

def _define_simple_least_squares(M: Model):
	"""
	Create a simple least squares function based on the steady state equations
	of a model. 

	Parameters
	----------
	M : Model
		the pynare model with a list of steady state functions in the 
		`_steady_state_funcs` attribute

	Returns
	-------
	least squares function
	"""

	ss_expr = vectorize_functions(M.steady_state_model)

	def _least_squares(x):
		return np.sum(np.power(ss_expr(x), 2))

	return _least_squares


def minimize_dogleg(
	M: Model,
	**kwargs
):
	"""
	Finds the steady state of a model using the dog-leg trust-region algorithm.
	This corresponds to the steady-state command
		'steady(solve_algo=0)'
	in dynare.

	Parameters
	----------
	initial_trust_radius : float
		initial trust-region radius
	max_trust_radius : float
		maximum value of the trust-region radius. No steps that are longer than
		this value will be proposed
	eta : float
		trust region related acceptance stringency for proposed steps
	gtol : float
		gradient norm must be less than 'gtol' before successful termination
	"""
	raise NotImplementedError('dog-leg algorithm')


def minimize_conjugate_gradient(
	M: Model,
	**kwargs
):
	"""
	Finds the steady state of a model using the conjugate gradient algorithm. 
	This corresponds most closely to the steady-state command 
		'steady(solve_algo=8)'
	in dynare. The correspondence between the dynare algorithm and the one 
	implemented here is not exact, as the dynare algorithm used is the Newton
	algorithm with a Stabilized Bi-Conjugate Gradient solver at each iteration

	Parameters
	----------
	M: Model
		the pynare model with valid steady state equations

	Returns 
	-------
	vector of steady-state values for endogenous variables in the order they
	were declared
	"""
	obj_func = _define_simple_least_squares(M)
	x0 = tuple([v.initial for v in M._endogenous])

	return optim.minimize(
		fun=obj_func,
		x0=x0,
		method='cg',
		**kwargs
	)


def minimize_nelder_mead(
	M: Model,
	**kwargs
):
	"""
	Finds the steady state of a model using the Nelder-Mead algorithm. This 
	doesn't correspond with any of dynare's steady-state solvers, but it's 
	the default algorithm for scipy's minimize() function, so I thought it 
	would be good to include it

	Parameters
	----------
	M : Model
		the pynare model with valid steady state equations

	Returns
	-------
	vector of steady-state values for endogenous variables in the order they
	were declared
	"""
	obj_func = _define_simple_least_squares(M)
	x0 = tuple([v.initial for v in M._endogenous])

	return optim.minimize(
		fun=obj_func,
		x0=x0,
		method='nelder-mead',
		**kwargs
	)






class _SteadyStateManager(object):

	# _dynare_default = 4
	_dynare_default = 8
	_dynare_mappings = {
		0: minimize_dogleg,
		1: 'dynares own nonlinear solver',
		2: 'split model into recurve blocks and perform (1)',
		3: 'chris sims solver',
		4: 'split model into blocks and use trust-region',
		5: 'newton algo with sparse Gaussian elimination',
		6: 'newton algo with sparse LU solver',
		7: 'newton algo with generalized minimal residual',
		8: minimize_conjugate_gradient,
		9: 'trust-region algo on entire model',
		10: 'levenberg-marquardt mixed with complementatiry problem',
		11: 'Ferris and Munson PATH solver'
	}


	_nonlang_mappings = {
		'nelder-mead': minimize_nelder_mead
	}

	def __new__(
		cls,
		language, 
		method
	):
		try:
			if language == 'dynare':
	
				# dynare allows selection with '0' int key, so we can't use a simple
				#	'select = method if method else cls._dynare_default'
				if isinstance(method, str):
					if method == '':
						select = cls._dynare_default
					else:
						raise KeyError('dynare steady state keys must be integers')
	
				elif isinstance(method, (float, int)):
					select = int(method)
	
				elif not method:
					select = cls._dynare_default
	
				return cls._dynare_mappings[select]

		except KeyError:
			return cls._nonlang_mappings[method]