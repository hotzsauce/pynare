"""
utility functions for numpy- and scipy-related operations
"""

from __future__ import annotations

import scipy
import numpy as np 

def find_first(
	vec: np.ndarray
):
	# very fast numpy-pure method for finding first index wehre true
	idx = vec.view(bool).argmax() // vec.itemsize
	return idx if vec[idx] else -1 


def partitioned_qz(
	A: np.ndarray,
	B: np.ndarray
):
	"""
	Given two (n x n) matrices, computes the QZ (aka Generalized Schur)
	decomposition & generalized eigenvectors. Then the resulting right unitary
	matrix and triangular matrices are partitioned according to the magnitude
	of the eigenvalues. The eigenvalues & partitioned matrices are then returned

	More specifically, given the eigenvalue problem A*x = s*B*x, the matrices 
	A and B are decomposed as
		A = Q*T*Z
		B = Q*S*Z
	where
		T is upper triangular,
		S is quasi upper triangular,
		Q and Z are unitary matrices

	For matrices X = T, S, Z, they are all partitioned as 
		X = [[X_11 X_12]
			  X_21 X_22]]
	according to whether the eigenvalues are greater than or less than unity in
	absolute value, with those associated with eigenvalues less than 1 in the
	upper left submatrix

	Parameters
	----------
	A : numpy 2d-array
		an (n x n) square matrix
	B : numpy 2d-array
		an (n x n) square matrix

	Returns
	-------
	eigs : numpy 1d-array
		n-vector of generalized eigenvalues
	Z_11 : numpy 2d-array
		upper left submatrix of unitary matrix Z
	Z_12 : numpy 2d-array
		upper right submatrix of unitary matrix Z
	Z_21 : numpy 2d-array
		lower left submatrix of unitary matrix Z
	Z_22 : numpy 2d-array
		lower right submatrix of unitary matrix Z
	T_11 : numpy 2d-array
		upper left submatrix of triangular matrix T
	T_12 : numpy 2d-array
		upper right submatrix of triangular matrix T
	T_22 : numpy 2d-array
		lower right submatrix of triangular matrix T
	S_11 : numpy 2d-array
		upper left submatrix of triangular matrix S
	S_12 : numpy 2d-array
		upper right submatrix of triangular matrix S
	S_22 : numpy 2d-array
		lower right submatrix of triangular matrix S
	"""

	S, T, alph, bet, Q, Zt = scipy.linalg.ordqz(B, A, sort='iuc')

	with np.errstate(divide='ignore', invalid='ignore'):
		# locally silence RuntimeWarnings about dividing by zero or infinity
		eigs = alph/bet

	# explosive eigenvalues
	expl = find_first((eigs > 1) | np.isinf(eigs))

	# unitary matrix first
	Z = np.transpose(Zt)
	Z_11 = Z[:expl, :expl]
	Z_12 = Z[:expl, expl:]
	Z_21 = Z[expl:, :expl]
	Z_22 = Z[expl:, expl:]

	# triangular matrix associated with matrix A
	T_11 = T[:expl, :expl]
	T_12 = T[:expl, expl:]
	T_22 = T[expl:, expl:]

	# triangular matrix associated with matrix B
	S_11 = S[:expl, :expl]
	S_12 = S[:expl, expl:]
	S_22 = S[expl:, expl:]	

	return eigs, Z_11, Z_12, Z_21, Z_22, T_11, T_12, T_22, S_11, S_12, S_22


def array_shingle(
	a0: np.ndarray,
	a1: np.ndarray, 
	overlap: int
):
	"""
	Merges two arrays along their first dimension by `covering up' the first 
	`overlap' rows of the second array with the last `overlap' rows of the 
	first array. If the two arrays are 2-dimensional, the analogy to roofing
	shingles is more clearly seen. All dimensions (except the first) of both 
	arrays must be equal, or else a ValueError is thrown.

	Parameters
	----------
	a0 : numpy ndarray
		the array on top that will cover up the bottom array
	a1 : numpy ndarray
		the bottom array that will be covered up by a0
	overlap : int
		the number of rows of the bottom array to be covered up

	Returns
	-------
	shingle : numpy ndarray
		the conjoined array

	Example
	-------
	>>> a = np.random.rand(5, 3)
	[[0.87599626 0.88090669 0.24360141]
 	 [0.21553097 0.57258704 0.75483751]
 	 [0.05710367 0.85748782 0.86286493]
 	 [0.92451993 0.63332986 0.96809019]
 	 [0.00214015 0.95635431 0.25054558]]

	>>> b = np.random.rand(7, 3)
	[[0.46366265 0.34500481 0.25021755]
 	 [0.2536837  0.02746545 0.92360626]
 	 [0.78506297 0.1625407  0.47006583]
 	 [0.80043353 0.33777573 0.1051361 ]
 	 [0.84375664 0.48336254 0.41655925]
 	 [0.88435366 0.76542689 0.77689997]
 	 [0.13782912 0.7792985  0.17253379]]

	>>> array_shingle(a, b, 2)
	[[0.87599626 0.88090669 0.24360141]
 	 [0.21553097 0.57258704 0.75483751]
 	 [0.05710367 0.85748782 0.86286493]
 	 [0.92451993 0.63332986 0.96809019]
 	 [0.00214015 0.95635431 0.25054558]
 	 [0.84375664 0.48336254 0.41655925]
 	 [0.88435366 0.76542689 0.77689997]
 	 [0.13782912 0.7792985  0.17253379]]
	"""
	
	s0, s1 = a0.shape, a1.shape

	if s0[1:] != s1[1:]:
		raise ValueError(
			'arrays must be have identical shapes, excluding first dim.'
		)
	
	shingle = np.zeros((s0[0] + s1[0] - overlap, *s0[1:]))
	shingle[:s0[0], ...] = a0
	shingle[s0[0]:, ...] = a1[overlap:, ...]

	return shingle 