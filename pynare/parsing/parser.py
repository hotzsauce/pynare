""" 
module for the most-extensive pynare Parser - "most-extensive" here meaning it
is the parser that is able to parse the largest subset of the Dynare/Matlab 
language. 
	There are other, smaller, Parsers in other pynare.parsing modules that
handle much simpler things, like the parser in the algebra module that can handle
expressions with '+', '-', '*', '/', '^', for example. The parser here inherits 
the methods of those simpler parsers, emphasizing the distinctions between 
different parts of a .mod file, and also exemplifying the inheritance and 
polymorphism of Python classes
"""

import pynare.parsing.algebra as alg


import pynare.parsing.ast as ast
import pynare.parsing.base as base

import pynare.parsing.model as mdl 
import pynare.parsing.variables as vbl 
import pynare.parsing.functions as fnc 
import pynare.parsing.simulation as sim

from pynare.parsing.tokens import Token 


class Parser(fnc.Parser):

	def __init__(self, lexer):
		super().__init__(lexer)

	# MODEL EXPRESSION METHODS - ONLY DIFFERENCE BETWEEN THESE ANALOGOUS ALBEGRA
	#	AND FUNCTION PARSER METHODS IS THE matom() METHOD ALLOWS FOR OFFSERT VARS.
	#	IN THE FUTURE, mcall() WILL BE ABLE TO RECOGNIZE EXPECTATION AND 
	#	STEADYSTATE OPERATORS
	def matom(self):
		token = self.current_token

		if token.type == base.ID:
			var = self.maybe_offset_variable()
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

	def mcall(self):
		while self.current_token.type == fnc.FUNCTION:
			token = self.current_token
			self.eat(FUNCTION)
			self.eat(base.LPARE)
			node = ast.Function(token=token,
								expr=self.mexpr())
			self.eat(base.RPARE)

		try:
			return node
		except UnboundLocalError:
			return self.matom()

	def mexponent(self):
		node = self.mcall()

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


	# OUT-OF- and IN-PERIOD VARIABLES
	def period_offset(self):
		"""
		period_offset : LPARE (PLUS | MINUS) NUMBER RPARE
		"""
		self.eat(base.LPARE)
		if self.current_token.type in (base.PLUS, base.MINUS):
			direction = self.current_token
			self.eat(direction.type)
		else:
			self.error()
		periods = self.current_token
		self.eat(base.NUMBER)
		self.eat(base.RPARE)

		# int() can't read the periods.value b/c it is a string float, so 
		# 	we first cast to float
		period_str = direction.value + str(periods.value)
		return int(float(period_str))

	def maybe_offset_variable(self):
		"""
		maybe_offset_variable : ID
							  | ID period_offset
		returns ast.OffsetVar if offset given, ast.Var node otherwise
		"""
		var = self.variable()
		if self.peek_type() == base.LPARE:
			return ast.OffsetVar(var=var, offset=self.period_offset())
		return var

	def variable(self):
		node = ast.Var(self.current_token)
		self.eat(base.ID)
		return node


	# BEFORE MODEL DECLARATION
	def variable_declaration(self, vtype):
		"""
		variable_declaration : ID (COMMA) variable_declaration
							 | empty
		"""
		# if a vtype keyword has been encountered, there must be at least one var
		var_nodes = [ast.Var(self.current_token)]
		self.eat(base.ID)

		# constructing nodes for any vars after the first
		while self.current_token.type in (base.COMMA, base.ID):

			# commas between variable names are optional
			if self.current_token.type == base.COMMA:
				self.eat(base.COMMA)
				var_nodes.append(self.variable())

			elif self.current_token.type == base.ID:
				var_nodes.append(self.variable())

		return [ast.VarDeclaration(v, vtype) for v in var_nodes]

	def vtype_declaration(self, vtypes):

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
		root = ast.Compound()
		root.children.extend(self.vtype_declaration(vbl.REQUIRED_DECLARATIONS))
		root.children.extend(self.vtype_declaration(vbl.OPTIONAL_DECLARATIONS))
		return root

	def variable_assignment(self):
		"""
		variable_assignment : ID EQUALS expr SEMI
		"""
		left = self.variable()
		token = self.current_token
		self.eat(base.EQUALS)
		right = self.expr()
		node = ast.VarAssignment(left, token, right)
		self.eat(base.SEMI)
		return node

	def assignment(self):
		"""
		assignment : (variable_assignment)+
		"""

		root = ast.Compound()
		while self.current_token.type == base.ID:
			node = self.variable_assignment()
			root.children.append(node)
		return root


	# MODEL DECLARATION
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
		#	than one it's expecting, but it throws a syntax error either way for now
		while self.current_token.type is not mdl.END:
			self.eat(base.SEMI)
			state = self.model_statement()
			# having the while conditino based on an expected token necessitates 
			#	this 'if', however
			if state.left:
				root.children.append(state)
		return root

	def model_statement(self):
		"""
		model_statement : local_declaration
						| tag model_expression
						| model_expression
		"""

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


	# AUXILLIARY FUNCTIONALITY IN MODEL STATEMENTS: LOCAL VARIABLES, TAGS
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


	# SETTING INITIAL AND TERMINAL VALUES FOR SOLVING & SIMULATING
	def model_condition_block(self):

		unused_types = sim.CONDITION_TYPES.copy()
		root = ast.ModelConditionBlock()

		while self.current_token.type in sim.CONDITION_TYPES:

			# ensuring each type of block is called exactly once
			if self.current_token.type not in unused_types:
				self.error()
			unused_types.remove(self.current_token.type)

			if self.current_token.type == sim.INITVAL:
				self.eat(sim.INITVAL)
				root.initial = self.model_condition_values()
			if self.current_token.type == sim.ENDVAL:
				self.eat(sim.ENDVAL)
				root.terminal = self.model_condition_values()
			if self.current_token.type == sim.HISTVAL:
				self.eat(sim.HISTVAL)
				root.historical = self.model_condition_values()

		return root
		
	def model_condition_values(self):
		param = self.optional_parameter_list()
		self.eat(base.SEMI)

		root = ast.ModelConditionValues(parameters=param)
		while self.current_token.type == base.ID:
			node = self.variable_assignment()
			root.children.append(node)
		self.eat(mdl.END)
		self.eat(base.SEMI)
		return root


	# PARSING THE BLOCK OF SHOCKS
	def shock_block(self):
		self.eat(sim.SHOCKS)
		self.eat(base.SEMI)
		node = self.shock_list()
		self.eat(mdl.END)
		self.eat(base.SEMI)
		return node
		
	def shock_list(self):

		root = ast.ShockBlock()
		root.children.append(self.shock())
		self.eat(base.SEMI)
		while self.current_token.type in (sim.CORR, vbl.VAR):
			root.children.append(self.shock())
			self.eat(base.SEMI)
		return root

	def shock(self):

		if self.current_token.type == sim.CORR:
			return self.correlated_shock()
		elif (
			self.current_token.type == vbl.VAR
			and self.peek_type(k=2) == base.COMMA
		):
			return self.covarying_shock() 
		elif (
			self.current_token.type == vbl.VAR
			and self.peek_type(k=2) == base.EQUALS
		):
			return self.variance_shock()
		elif (
			self.current_token.type == vbl.VAR
			and self.peek_type(k=2) == base.SEMI
			and self.peek_type(k=3) == sim.STDERR
		):
			return self.stderr_shock()

		else:
			return self.deterministic_shock()

	def correlated_shock(self):
		"""
		correlated_shock : CORR ID COMMA ID EQUALS expr
		"""	
		self.eat(sim.CORR)
		x = self.variable()
		self.eat(base.COMMA)
		y = self.variable()
		self.eat(base.EQUALS)
		expr = self.expr()

		node = ast.CorrelatedShock(
					left=x,
					right=y,
					expr=expr)
		return node

	def covarying_shock(self):
		"""
		covarying_shock : VAR ID COMMA ID EQUALS expr
		"""
		self.eat(vbl.VAR)
		x = self.variable()
		self.eat(base.COMMA)
		y = self.variable()
		self.eat(base.EQUALS)
		expr = self.expr()

		node = ast.CovaryingShock(
					left=x,
					right=y,
					expr=expr)
		return node

	def variance_shock(self):
		"""
		variance_shock : VAR ID EQUALS expr
		"""
		self.eat(vbl.VAR)
		token = self.variable()
		self.eat(base.EQUALS)
		expr = self.expr()

		node = ast.VarianceShock(
					token=token,
					expr=expr)
		return node

	def stderr_shock(self):
		"""
		stderr_shock : VAR ID SEMI STDERR expr
		"""
		self.eat(vbl.VAR)
		token = self.variable()
		self.eat(base.SEMI)
		self.eat(sim.STDERR)
		expr = self.expr()

		node = ast.StdErrShock(
					token=token,
					expr=expr)
		return node

	def deterministic_shock(self):
		"""
		deterministic_shock : VAR ID SEMI PERIODS periods_list SEMI VALUES values_list
		"""
		self.eat(vbl.VAR)
		token = self.variable()
		self.eat(base.SEMI)
		self.eat(sim.PERIODS)
		periods = self.shock_periods_list()
		self.eat(base.SEMI)
		self.eat(sim.VALUES)
		values = self.shock_values_list()

		node = ast.DeterministicShock(
					token=token,
					periods=periods,
					values=values)
		return node

	def shock_periods_list(self):
		"""
		shock_periods_list : shock_period COMMA shock_periods_list
						   | shock_period shock_periods_list
						   | shock_period
		"""
		root = ast.Compound()
		root.children.append(self.shock_period())
		while self.current_token.type in (base.COMMA, base.NUMBER):

			if self.current_token.type == base.COMMA:
				self.eat(base.COMMA)
				root.children.append(self.shock_period())

			elif self.current_token.type == base.NUMBER:
				root.children.append(self.shock_period())

		return root

	def shock_period(self):
		"""
		shock_period : NUMBER COLON NUMBER
					 | NUMBER
		"""
		start = self.current_token
		self.eat(base.NUMBER)
		if self.current_token.type == base.COLON:
			self.eat(base.COLON)
			end = self.current_token
			self.eat(base.NUMBER)
			return ast.ShockPeriod(start=start, end=end)
		return ast.ShockPeriod(start=start, end=start)

	def shock_values_list(self):
		"""
		shock_values_list : shock_value COMMA shock_values_list
						  | shock_value shock_values_list
						  | shock_value
		"""
		root = ast.Compound()
		root.children.append(self.shock_value())
		while self.current_token.type in (base.COMMA, base.NUMBER, base.LPARE):

			if self.current_token.type == base.COMMA:
				self.eat(base.COMMA)
				root.children.append(self.shock_value())

			else:
				root.children.append(self.shock_value())

		return root

	def shock_value(self):
		"""
		shock_value : NUM
					| LPARE expr RPARE
		"""
		if self.current_token.type == base.NUMBER:
			token = self.current_token
			self.eat(base.NUMBER)
			return ast.ShockValue(token=token)
		self.eat(base.LPARE)
		expr = self.expr()
		self.eat(base.RPARE)
		return ast.ShockValue(token=expr)
	
	
	# PARSING LISTS OF PARAMETERS - NOT EXCLUSIVE TO MODEL BLOCK, BUT IT'S 
	#	THE FIRST TIME IT CAME UP
	def parameter_list(self):
		"""
		parameter_list : LPARE ID (COMMA ID)* RPARE
		"""
		self.eat(base.LPARE)

		root = ast.Compound()
		root.children.append(ast.Param(self.current_token))
		self.eat(base.ID)

		while self.current_token.type == base.COMMA:
			self.eat(base.COMMA)
			root.children.append(ast.Param(self.current_token))
			self.eat(base.ID)
		self.eat(base.RPARE)

		return root

	def optional_parameter_list(self):
		""" auxillary function to check for parameter list """
		if self.peek_type() == base.LPARE:
			return self.parameter_list()
		else:
			return ast.Compound()


	def modfile(self):
		decl_node = self.declaration()
		assn_node = self.assignment()

		model_block_node = self.model_block()
		cond_node = self.model_condition_block()

		shock_node = self.shock_block()

		return ast.ModFile(declaration=decl_node, 
							assignment=assn_node, 
							model_block=model_block_node,
							simconditions=cond_node
		)
		
	def parse(self):
		return self.modfile()