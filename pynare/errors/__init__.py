""" pynare errors """

class PynareError(Exception):
	pass


class PynareWarning(object):
	def __init__(self):
		pass

	def issue_warning(self, msg):
		issued = ': '.join([type(self).__name__, msg])
		print(issued)


class PynareSyntaxError(PynareError):
	def __init__(self, line, pos):
		self.line = line
		self.pos = pos
		super().__init__()