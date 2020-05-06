

from pynare.errors import PynareSyntaxError

class ABCParser(object):
	""" Abstract Base Class for parsers throughout Pynare """

	def __init__(self, lexer):
		self.lexer = lexer
		self.current_token = self.lexer.get_next_token()

	def error(self):
		raise PynareSyntaxError(self.lexer.text, 20)

	def type_check(self, *args):
		if self.current_token.type not in args:
			raise TypeError(self.current_token.type)

	def eat(self, token_type):
		"""
		compare the current token type wit hthe passed token type and if they match
		eat the current token and assign the next token to the self.current_token.
		Otherwise, throw an error
		"""
		if self.current_token.type == token_type:
			self.current_token = self.lexer.get_next_token()
		else:
			self.error()

	def peek(self):
		return self.lexer.see_next_token()

	def peek_type(self):
		# just because I don't like seeing 'self.peek().type ==' in higher-level
		#	parsers
		return self.peek().type



class ABCVisitor(object):
	""" Abstract Base Class for AST node visitors """

	def __init__(self, tree):
		self.tree = tree

	def visit(self, node):
		method = 'visit_' + type(node).__name__
		visitor = getattr(self, method, self._nonexistent_node)
		return visitor(node)

	def _nonexistent_node(self, node):
		raise Exception('No visit_{} method'.format(type(node).__name__))