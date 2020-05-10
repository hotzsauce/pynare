

from pynare.parsing.lexer import Lexer
from pynare.parsing.parser import Parser
from pynare.parsing.interpreter import Interpreter
from pynare.parsing.semantics import SemanticAnalyzer



class Tester(object):
	
	def __init__(self, text):
		self.lexer = Lexer(text)
		self.parser = Parser(self.lexer)


	def read(self):
		self.check_semantics()

		# self.interpreter = Interpreter(self.tree)
		# return self.interpreter.interpret()


	def check_semantics(self):
		self.sem_analyzer = SemanticAnalyzer()
		self.sem_analyzer.visit(self.tree)


	def summarize_tokens(self):
		for t in self.lexer.token_stream:
			print(t)

	def show_parse_tree(self):
		tree = self.tree
		tree.describe()

	@property
	def tree(self):
		if hasattr(self, '_tree'):
			return self._tree
		else:
			self._tree = self.parser.parse()
			return self._tree