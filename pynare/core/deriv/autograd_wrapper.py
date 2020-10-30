import numpy as np
import functools

from autograd import grad

def grad_array(func):

	_grad = grad(func)

	@functools.wraps(func)
	def _grad_array(x):
		return np.array(_grad(x))

	return _grad_array