"""
first-order solver
"""

from __future__ import annotations

from typing import Tuple

from scipy import linalg

import numpy as np 
import pandas as pd 

import pynare.utils.numpy as unp

from pynare.errors import ModelIdentificationError



def solve_first_order(
	M: Model
) -> Tuple[np.ndarray]:
	"""
	Computes the first-order linear solution to a model. Notation follows
	Villemot (2011). A nonlinear model can be written as 
				f(y_tp1^p, y_t, y_tm1^m, u_t) = 0,
	where 
		y_t 	- vector of endogenous variables
		y_t^p 	- subset of y_t that appears with a lead
		y_t^m 	- subset of y_t that appears with a lag
		u_t 	- vector of exogenous shocks. 
	If the model is expressed as above, the solution g solves the functional 
	equation
				f(g(g(y_tm1, u_t), u_tp1), g(y_tm1, u_t), y_tm1, u_t) = 0.
	This function solves for the first-order solution, i.e.
				g(y_tm1, u_t) = ybar + gy*yhat_tm1 + gu*uhat_t
	where the 'hat' suffix denotes percent deviations from the model's steady
	state, ybar.

	Parameters
	----------
	M : Model
		a pynare model

	Returns
	-------
	tuple of numpy ndarrays
		the tuple of vectors gy, gu described above
	"""

	idx_mgr = M._index_mgr

	J = M.jacobian
	S = J[:, idx_mgr.static_jacobian]

	Q, R = np.linalg.qr(S, mode='complete')


	# for model to be properly, identified, rank(R) must equal number of static 
	_r = np.linalg.matrix_rank(R)
	if _r != idx_mgr.n_static:
		raise ModelIdentificationError(_r, idx_mgr.n_static)

	A = np.dot(np.transpose(Q), J)

	# solving for backward- and forward-looking dynamic variables & merging 
	# 	the arrays. by construction the bottom n_mixed rows of gym & top 
	# n_mixed rows of gyp are the same
	gym, gyp = solve_nonstatic_endo(A, idx_mgr)
	gyd = unp.array_shingle(gym, gyp, idx_mgr.n_mixed)


	# solving for static variables, then stacking the static and dynamic arrays
	gys = solve_static_endo(A, idx_mgr, gym, gyp)
	gy = np.vstack((gys, gyd))

	fu = J[:, idx_mgr.exogenous_jacobian]
	gu = solve_exo(A, idx_mgr, gyp, fu)

	return gy, gu
	


def solve_nonstatic_endo(
	A: np.ndarray,
	idx_mgr: VariableIndexManager
):
	"""
	Solving for the non-static endogenous variables of the model. This corresponds
	to Section 4.1 of Villemot (2011)

	Parameters
	----------
	A: np.ndarray
		the rotated Jacobian matrix of the model
	idx_mgr : VariableIndexManager
		the index manager of the model

	Returns
	-------
	2-tuple of np.ndarray
		respectively, the state response to backward- and forward-looking dynamic
		endogenous variables
	"""
	nd = idx_mgr.n_dynamic
	Ap_tilde = A[-nd:, idx_mgr.forward_jacobian]
	Am_tilde = A[-nd:, idx_mgr.backward_jacobian]
	A0p_tilde = A[-nd:, idx_mgr.cont_forward_jacobian]

	# model might not have any purely backward variables
	A0m_tilde = A[-nd:, idx_mgr.pure_backward_jacobian]
	if not A0m_tilde:
		A0m_tilde = np.zeros((nd, 1))

	# setting up structural state space representation
	nm = idx_mgr.n_mixed
	nb = len(idx_mgr.backward_jacobian)
	nf = len(idx_mgr.forward_jacobian)

	Im = np.eye(nm, nb)
	Ip = np.eye(nm, nf)

	D = np.block([[A0m_tilde, Ap_tilde], [Im, np.zeros((nm, nf))]])
	E = np.block([[-Am_tilde, -A0p_tilde], [np.zeros((nm, nb)), Ip]])

	# solving the generalized Schur decomposition
	eigs, Z_11, _, Z_21, Z_22, T_11, _, _, S_11, _, _, = unp.partitioned_qz(D, E)

	# checking for the Blanchard-Kahn (1980) order and rank conditions 
	n_explosive = np.sum((eigs > 1) | np.isinf(eigs))
	n_forward = len(idx_mgr.forward_jacobian)
	if n_explosive != n_forward:
		raise ModelIdentificationError('order', n_explosive, n_forward)

	if np.linalg.matrix_rank(Z_22) != Z_22.shape[0]:
		raise ModelIdentificationError('rank')

	# solution for forward-looking variables
	gyp = - np.linalg.solve(Z_22, Z_21)

	# solution for backward-looking variables
	Z_11t = np.transpose(Z_11)
	gym = np.dot(
		np.dot(Z_11t, np.linalg.inv(T_11)),
		np.dot(S_11, np.linalg.inv(Z_11t))
	)

	return gym, gyp


def solve_static_endo(
	A: np.ndarray,
	idx_mgr: VariableIndexManager,
	gym: np.ndarray,
	gyp: np.ndarray
):
	"""
	Solving for the static endogenous variables of the model. This corresponds
	to Section 4.2 of Villemot (2011)

	Parameters
	----------
	A: numpy ndarray
		the rotated Jacobian matrix of the model
	idx_mgr : VariableIndexManager
		the index manager of the model
	gym: numpy ndarray
		the state response to backward-looking dynamic endogenous variables
	gyp: numpy ndarray
		the state response to forward-looking dynamic endogenous variables

	Returns
	-------
	np.ndarray
		the state response to static endogenous variables
	"""
	ns = idx_mgr.n_static

	Ap_cup = A[:ns, idx_mgr.forward_jacobian]
	Am_cup = A[:ns, idx_mgr.backward_jacobian]
	A0s_cup = A[:ns, idx_mgr.static_jacobian]
	A0d_cup = A[:ns, idx_mgr.dynamic_jacobian]

	# merging forward- and backward-looking dynamic response
	gyd = unp.array_shingle(gym, gyp, idx_mgr.n_mixed)

	Agp = np.dot(np.dot(Ap_cup, gyp), gym)
	Agd = np.dot(A0d_cup, gyd)

	return -np.dot(np.linalg.inv(A0s_cup), Agp + Agd + Am_cup)


def solve_exo(
	A: np.ndarray,
	idx_mgr: VariableIndexManager,
	gyp: np.ndarray,
	fu: np.ndarray
):
	"""
	Solving for the response to exogenous variables of the model. This corresponds
	to Section  5 of Villemot (2011)

	Parameters
	----------
	A: numpy ndarray
		the rotated Jacobian matrix of the model
	idx_mgr : VariableIndexManager
		the index manager of the model
	gyp: numpy ndarray
		the state response to forward-looking dynamic endogenous variables
	fu : numpy ndarray
		the (unrotated) Jacobian column(s) corresponding to exogenous responses

	Returns
	-------
	np.ndarray
		the state response to exogenous variables
	"""
	cont_static = A[:, idx_mgr.static_jacobian]
	cont_backward = A[:, idx_mgr.cont_backward_jacobian]
	cont_forward = A[:, idx_mgr.cont_pure_forward_jacobian]

	Ap = A[:, idx_mgr.forward_jacobian]
	implied_backward = np.dot(Ap, gyp) + cont_backward

	U = np.hstack((cont_static, implied_backward, cont_forward))
	return - np.linalg.solve(U, fu)