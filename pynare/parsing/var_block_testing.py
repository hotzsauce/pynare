
import ast
import base 

import model as mdl
import algebra as alg
import variables as vbl 
import functions as fnc 

from tokens import Token
from errors import PynareSyntaxError
from semantics import SemanticAnalyzer





RESERVED_KEYWORDS = {
	**vbl.RESERVED_KEYWORDS, 
	**mdl.RESERVED_KEYWORDS,
	**fnc.RESERVED_KEYWORDS
}



class Lexer(object):

	def __init__(self, text):
		self.text = text

		self.token_pos = self.pos = 0
		self.current_char = self.text[self.pos]
		self.token_stream = list()

		self.tokenize()


	def error(self, c):
		msg = 'Invalid character: {}'.format(repr(c))
		raise Exception(msg)

	def skip_whitespace(self):
		while self.current_char is not None and self.current_char.isspace():
			self.advance()

	def skip_comment(self):
		self.advance()
		while self.current_char is not '\n':
			self.advance()
		self.advance()
		self.skip_whitespace()

	def skip_multiline_comment(self):
		self.advance() # skip '/'
		self.advance() # skip '*'
		while (self.current_char != '*') and (self.peek() != '/'):
			self.advance()
		self.advance() # skip '*'
		self.advance() # skip '/' 
		self.skip_whitespace() # just in case

	def _id(self):
		result = ''
		while base.is_valid_varchar(self.current_char):
			result += self.current_char
			self.advance()
		token = RESERVED_KEYWORDS.get(result, Token(base.ID, result))
		return token

	def number(self):
		result = ''
		while self.current_char is not None and base.is_valid_numchar(self.current_char):
			result += self.current_char
			self.advance()

			if self.current_char == '.':
				result += self.current_char
				self.advance()

				while (
					self.current_char is not None and
					self.current_char.isdigit()
				):
					result += self.current_char
					self.advance()

		return Token(base.NUMBER, float(result))

	def string(self, quote):
		self.advance()
		result = ''
		while self.current_char is not quote:
			result += self.current_char
			self.advance() 
		self.advance()
		return Token(base.STRING, result)

	def advance(self):
		""" Advance the position and set the current character """
		self.pos += 1
		try:
			self.current_char = self.text[self.pos]
		except IndexError:
			self.current_char = None

	def peek(self):
		""" returns the next character, without incrementing the curreot one """
		peek_pos = self.pos + 1
		try:
			return self.text[peek_pos]
		except IndexError:
			return None

	def next_token(self):

		while self.current_char is not None:

			if self.current_char.isspace():
				self.skip_whitespace()

			if self.current_char == '%':
				self.skip_comment()

			if (
				self.current_char == '/'
				and self.peek() == '*'
			):
				self.skip_multiline_comment()

			if base.is_valid_varname_firstchar(self.current_char):
				return self._id()

			if base.is_valid_numchar(self.current_char):
				return self.number()

			if self.current_char == '\'':
				return self.string(quote='\'')

			if self.current_char == '\"':
				return self.string(quote='\"')

			if self.current_char == ',':
				self.advance()
				return Token(base.COMMA, ',')

			if self.current_char == ';':
				self.advance()
				return Token(base.SEMI, ';')

			if self.current_char == '+':
				self.advance()
				return Token(base.PLUS, '+')

			if self.current_char == '-':
				self.advance()
				return Token(base.MINUS, '-')

			if self.current_char == '*':
				self.advance()
				return Token(base.MUL, '*')

			if self.current_char == '/':
				self.advance()
				return Token(base.DIV, '/')

			if self.current_char == '^':
				self.advance()
				return Token(base.POWER, '^')

			if self.current_char == '(':
				self.advance()
				return Token(base.LPARE, '(')

			if self.current_char == ')':
				self.advance()
				return Token(base.RPARE, ')')

			if self.current_char == '[':
				self.advance()
				return Token(base.LBRACKET, '[')

			if self.current_char == ']':
				self.advance()
				return Token(base.RBRACKET, ']')

			if self.current_char == '#':
				self.advance()
				return Token(base.POUND, '#')

			if self.current_char == '=':
				self.advance()
				return Token(base.EQUALS, '=')

			if self.current_char is not None:
				self.error(self.current_char)

	def tokenize(self):

		while self.current_char is not None:
			t = self.next_token()
			# 'if' only necessary b/c self.next_token() returns None when Lexer 
			#	reaches the end of the source code. Probably indicates this 
			#	part needs to be refactored sometime
			if isinstance(t, Token):
				self.token_stream.append(t)
		self.token_stream.append(Token(base.EOF, None))

	def get_next_token(self):
		"""
		Returns the token at the current token position (it's the next token from
		the Parser's perspective, however, hence the name), and iterates the token
		position by one
		"""
		token = self.token_stream[self.token_pos]
		self.token_pos += 1 
		return token

	def see_future_token(self, k):
		future_pos = (self.token_pos - 1) + k
		if future_pos > len(self.token_stream) - 1:
			raise SyntaxError('Trying to look too far into the future.')
		else:
			return self.token_stream[future_pos]

	def see_next_token(self):
		"""
		Returns the token at the current token position (it's the next token from
		the Parser's perspective, however, hence the name), but unlike 
		self.get_next_token, the position is not incremented
		"""
		return self.see_future_token(k=0)



###############################################################################
#
# 	PARSER
#
###############################################################################

class Parser(fnc.FunctionParser):

	def __init__(self, lexer):
		super().__init__(lexer)


	def parameter_list(self):
		self.eat(base.LPARE)

		# must have at least one parameter
		paras = [ast.Param(self.current_token)]
		self.eat(base.ID)

		while self.current_token.type == base.COMMA:
			self.eat(base.COMMA)
			paras.append(ast.Param(self.current_token))
			self.eat(base.ID)	
		self.eat(base.RPARE)

		return paras

	def optional_parameter_list(self):
		if self.peek_type() == base.LPARE:
			return self.parameter_list()
		else:
			return list()


	def model_block(self):
		"""
		model_block : MODEL SEMI model END SEMI
					| MODEL LPARE parameter_list RPARE SEMI model END SEMI
		"""
		
		self.eat(mdl.MODEL)
		parameters = self.optional_parameter_list()
		self.eat(base.SEMI)

		model = self.model()
		self.eat(mdl.END)
		self.eat(base.SEMI)

		return ast.ModelBlock(parameters, model)

	def model(self):
		"""
		model : model_statement SEMI (model_statement SEMI)*
		"""

		root = ast.Model()
		root.children.append(self.model_statement())

		# i would like to have this condition based on the current token, rather
		# 	than one it's expecting, but it throws a syntax error either way
		while self.current_token.type is not mdl.END: 
			self.eat(base.SEMI)
			state = self.model_statement()
			# having the while condition based on an expected token necessitates
			# 	this 'if', however:
			if state.left:
				root.children.append(state)
			
		return root


	def model_statement(self):

		if self.current_token.type == base.POUND:
			return self.local_declaration()
		elif self.current_token.type == base.LBRACKET:
			tag = self.tag()
			node = self.model_expression()
			node.tag = tag
			return node
		else:
			return self.model_expression()

	def model_expression(self):
		"""
		model_expression : mexpr
						 | mexpr EQUALS mexpr
		"""
		left = self.mexpr()
		if self.peek_type() == base.EQUALS:
			self.eat(base.EQUALS)
			right = self.mexpr()
		else:
			# if there's no equals sign, it's assumed this is a homogenous equation
			zero = Token(base.NUMBER, 0)
			right = ast.Num(zero)

		return ast.ModelExpression(left, right)

	def local_declaration(self):
		"""
		local_declaration : POUND ID EQUALS mexpr
		"""
		self.eat(base.POUND)

		left = ast.Var(self.current_token)
		self.eat(base.ID)
		token = self.current_token
		self.eat(base.EQUALS)
		right = self.mexpr()

		return ast.VarAssignment(left, token, right)

	def tag(self):
		"""
		tag : LBRACKET tag_list RBRACKET
		"""
		self.eat(base.LBRACKET)
		root = ast.Tag()
		root.children.append(self.tag_pair())

		# tag_list doesn't have it's own method; its done right here
		while self.current_token.type == base.COMMA:
			self.eat(base.COMMA)
			root.children.append(self.tag_pair())
		self.eat(base.RBRACKET)

		return root

	def tag_pair(self):
		"""
		tag_pair : ID EQUALS STRING
		"""
		key = self.current_token
		self.eat(base.ID)
		self.eat(base.EQUALS)
		value = self.current_token
		self.eat(base.STRING)
		return ast.TagPair(key, value)

	

	# MODEL
	def matom(self):
		token = self.current_token

		if token.type == base.ID:
			var = self.variable()
			if self.peek().type == base.LPARE:
				pvar = self.period_variable(var)
				# return pvar
			return var
		elif token.type == base.NUMBER:
			self.eat(base.NUMBER)
			return ast.Num(token)
		elif token.type == base.LPARE:
			self.eat(base.LPARE)
			node = self.mexpr()
			self.eat(base.RPARE)
			return node
		elif token.type in (base.MINUS, base.PLUS):
			self.eat(token.type)
			node = ast.UnaryOp(op=token,
							expr=self.mterm())
			return node

	def mexponent(self):
		node = self.matom()

		while self.current_token.type == base.POWER:
			token = self.current_token
			self.eat(base.POWER)

			node = ast.BinaryOp(left=node,
							op=token,
							right=self.mexponent())
		return node

	def mterm(self):
		node = self.mexponent()

		while self.current_token.type in (base.MUL, base.DIV):
			token = self.current_token
			self.eat(token.type)

			node = ast.BinaryOp(left=node,
							op=token,
							right=self.mexponent())
		return node

	def mexpr(self):
		node = self.mterm()

		while self.current_token.type in (base.PLUS, base.MINUS):
			token = self.current_token
			self.eat(token.type)

			node = ast.BinaryOp(left=node,
							op=token,
							right=self.mterm())
		return node

	# MODFILE 
	def period_variable(self, var):
		self.eat(base.LPARE)
		if self.current_token.type in (base.PLUS, base.MINUS):
			direction = self.current_token
			self.eat(direction.type)
		else:
			self.error()

		periods = self.current_token
		period_str = direction.value + str(periods.value)
		self.eat(base.NUMBER)
		self.eat(base.RPARE)

		# int() can't read the periods.value b/c it is a string float,
		# 	so first we cast to float
		offset = int(float(period_str))
		return ast.PeriodVar(var.token, offset)

	def variable(self):
		node = ast.Var(self.current_token)
		self.eat(base.ID)
		return node

	def variable_declaration(
		self,
		vtype
	):
		# if a vtype keyword has been encountered, there must be at least one variable
		var_nodes = [ast.Var(self.current_token)]
		self.eat(base.ID)

		# constructing nodes for variables after the first
		while self.current_token.type in (base.COMMA, base.ID):
			
			if self.current_token.type == base.COMMA:
				self.eat(base.COMMA)
				var_nodes.append(ast.Var(self.current_token))
				self.eat(base.ID)

			# commas between variables are optional
			elif self.current_token.type == base.ID:
				var_nodes.append(ast.Var(self.current_token))
				self.eat(base.ID)

		return [ast.VarDeclaration(v, vtype) for v in var_nodes]

	def _declaration(self, vtypes):

		var_decl = list()
		while self.current_token.type in vtypes:
			# save the variable type to pass to variable_declaration
			vtype = self.current_token.type 

			# save the list of variables of type 'vtype'
			self.eat(self.current_token.type)
			vds = self.variable_declaration(vtype=vtype)
			var_decl.extend(vds)
			self.eat(base.SEMI)

		return var_decl

	def declaration(self):

		var_decl = list()
		var_decl.extend(self._declaration(vbl.REQUIRED_DECLARATIONS))
		var_decl.extend(self._declaration(vbl.OPTIONAL_DECLARATIONS))
		return var_decl

	def variable_assignment(self):

		left = self.variable()
		token = self.current_token
		self.eat(base.EQUALS)
		right = self.expr()
		node = ast.VarAssignment(left, token, right)
		self.eat(base.SEMI)
		return node

	def assignment(self):

		root = ast.Compound()
		while self.current_token.type == base.ID:
			node = self.variable_assignment()
			root.children.append(node)

		return root

	def modfile(self):
		decl_node = self.declaration()
		assn_node = self.assignment()
		model_node = self.model_block()
		# model_node = self.model_declaration()
		# model_node = None
		return ast.ModFile(decl_node, assn_node, model_node)
		
	def parse(self):
		return self.modfile()



###############################################################################
#
# 	INTERPRETER
#
###############################################################################

class Interpreter(fnc.FunctionInterpreter):

	def __init__(self, tree):

		self.tree = tree
		self._GLOBAL_SCOPE = dict()


	# VARIABLES STUFF
	def summarize_variables(self):
		print('\n')
		max_name_len = max([len(var) for var in self._GLOBAL_SCOPE])
		max_value_len = max([len(str(val)) for val in self._GLOBAL_SCOPE.values()])
		for k, v in self._GLOBAL_SCOPE.items():
			srm = '{var:<{var_width}} : {val:>{val_width}}'.format(
					var=k, var_width=max_name_len, 
					val=v, val_width=max_value_len
			)
			print(srm)


	def visit_VarAssignment(self, node):
		var_name = node.left.value
		self._GLOBAL_SCOPE[var_name] = self.visit(node.right) # to ALGEBRA method, or visit_Var

	def visit_Compound(self, node):
		for child in node.children:
			self.visit(child) 		# visit_VarAssignment

	def visit_VarDeclaration(self, node):
		pass

	def visit_ModFile(self, node):
		for vd in node.declaration:
			self.visit(vd) 			# visit_VarDeclaration
		self.visit(node.assignment) # visit_Compound

	def interpret(self):
		return self.visit(self.tree)



class Tester(object):

	def __init__(self, text):
		self.lexer = Lexer(text)
		self.parser = Parser(self.lexer)

	def read(self):
		tree = self.tree

		self.sem_analyzer = SemanticAnalyzer()
		self.sem_analyzer.visit(tree)

		self.interpreter = Interpreter(tree)
		result = self.interpreter.interpret()

		self.interpreter.summarize_variables()

		# print(self.sem_analyzer.scope)

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

