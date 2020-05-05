import pynare.parsing.functions as fnc

class Interpreter(fnc.Interpreter):

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