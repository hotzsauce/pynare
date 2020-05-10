

from pynare.parsing.variables import (
	VAR,
	VAREXO,
	VAREXO_DET,
	PARAMETERS,
	PREDETERMINED,
	TREND_VAR,
	LOG_TREND_VAR
)

from pynare.parsing.generic import ABCVisitor



class Symbol(object):
	def __init__(self, name, stype=None):
		self.name = name
		self.stype = stype


class BuiltinTypeSymbol(Symbol):
	def __init__(self, name):
		super().__init__(name)

	def __str__(self):
		return self.name

	def __repr__(self):
		return '<{class_name}({name})>'.format(
				class_name=self.__class__.__name__,
				name=self.name
		)


class VarSymbol(Symbol):
	def __init__(self, name, stype):
		super().__init__(name, stype)

	def __str__(self):
		return self.__repr__()

	def __repr__(self):
		return '<{class_name}({name}: {stype})>'.format(
				class_name=self.__class__.__name__,
				name=self.name,
				stype=self.stype
		)


class ScopedSymbolTable(object):

	def __init__(self, scope_name, scope_level):
		self._symbols = {}
		self.scope_name = scope_name
		self.scope_level = scope_level
		self._init_builtins()

	def _init_builtins(self):
		self.insert(BuiltinTypeSymbol(VAR))
		self.insert(BuiltinTypeSymbol(VAREXO))
		self.insert(BuiltinTypeSymbol(VAREXO_DET))
		self.insert(BuiltinTypeSymbol(PARAMETERS))
		self.insert(BuiltinTypeSymbol(PREDETERMINED))
		self.insert(BuiltinTypeSymbol(TREND_VAR))
		self.insert(BuiltinTypeSymbol(LOG_TREND_VAR))

	def insert(self, symbol):
		self._symbols[symbol.name] = symbol

	def lookup(self, name):
		return self._symbols.get(name)

	def __str__(self):

		n_delims = 50
		h1 = 'SCOPE (Scoped Symbol Table)'
		lines = ['\n', h1, '='*n_delims]
		for header_name, header_value in (
			('Scope Name', self.scope_name),
			('Scope Level', self.scope_level)
		):
			lines.append('{:<15}: {}'.format(header_name, header_value))

		h2 = 'Scope (Scoped Symbol Table) Contents'
		lines.extend([h2, '-'*n_delims])
		lines.extend(
			(('{:>15}: {}').format(k, v)
			for k, v in self._symbols.items()
			if not isinstance(v, BuiltinTypeSymbol))
		)
		lines.append('\n')
		return '\n'.join(lines)


class SemanticAnalyzer(ABCVisitor):

	def __init__(self):
		self.scope = ScopedSymbolTable(
						scope_name='global',
						scope_level=1
		)

	# ALGEBRA & FUNCTIONS
	def visit_BinaryOp(self, node):
		"""
		checks for semantic errors in left and right exprs
		"""
		self.visit(node.left)
		self.visit(node.right)

	def visit_UnaryOp(self, node):
		"""
		checks for semantic errors in expr
		"""
		self.visit(node.expr)

	def visit_Num(self, node):
		"""
		nothing to see here
		"""
		pass

	def visit_Var(self, node):
		"""
		checks that the variable has been defined. This check is identical to the 
		one in visit_VarDeclaration, but the description in the NameError differs
		"""
		var_name = node.value 
		var_symbol = self.scope.lookup(var_name)
		if var_symbol is None:
			raise NameError('name {} is not defined.'.format(repr(var_name)))

	def visit_Function(self, node):
		"""
		check the argument of the function
		"""
		self.visit(node.expr)


	# PARAMETERS
	def visit_Param(self, node):
		"""
		nothing to do just yet, will need to check that these are accepted
		parameters somehow, though
		"""
		pass

	# BLOCKS BEFORE MODEL
	def visit_VarDeclaration(self, node):
		"""
		retrieve the BuiltinSymbol based on 'vtype' attr of node, create a 
		VarSymbol based on that type and the name of the variable, and save
		to the Global Scope
		"""
		# getting built-in type
		vtype = node.vtype
		stype = self.scope.lookup(vtype)

		# variable name that's being declared & assigning its built-in type
		var_name = node.var_node.value
		var_symbol = VarSymbol(var_name, stype)

		# if this variable name has already been used, throw an error
		if self.scope.lookup(var_name) is not None:
			msg = 'name {} has already been declared.'
			raise NameError(msg.format(repr(var_name)))
		self.scope.insert(var_symbol)

	def visit_VarAssignment(self, node):
		"""
		checks variable has been declared & looks for semantic errors in its assigned 
		value (node.right). This check is identical to the one in visit_Var, but the 
		description in the NameError differs
		"""
		var_name = node.left.value
		var_symbol = self.scope.lookup(var_name)
		if var_symbol is None:
			msg = 'name {} was not declared before assignment.'
			raise NameError(msg.format(repr(var_name)))
		self.visit(node.right)


	# MODEL BLOCK
	def visit_ModelBlock(self, node):
		self.visit(node.parameters)
		self.visit(node.model)

	def visit_Model(self, node):
		node.describe()



	# MODFILE
	def visit_ModFile(self, node):
		"""
		checks the variable declarations node (a Compound), the variable
		assignments node (another Compound), and the model block (a 
		ModelBlock node)
		"""
		self.visit(node.declaration)
		self.visit(node.assignment)
		self.visit(node.model_block)

	def visit_Compound(self, node):
		for child in node.children:
			self.visit(child)


	def __repr__(self):
		return '<{class_name} Object>'.format(
				class_name=self.__class__.__name__
		)

	def __str__(self):
		header = '\n{class_name} with Underlying Table:'.format(
					class_name=self.__class__.__name__
		)
		return ''.join([header, self.scope.__str__()])