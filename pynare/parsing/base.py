""" 
module with keywords & key symbols, and validators used throughout 
pynare Parser 
"""

COMMA 		= 'COMMA' 
SEMI 		= 'SEMI'
COLON 		= 'COLON'
EOF 		= 'EOF'
ID 			= 'ID'
NUMBER 		= 'NUMBER'
STRING 		= 'STRING'

LPARE 		= '('
RPARE 		= ')'
LBRACKET 	= '['
RBRACKET 	= ']'
SQUOTE 		= '\''
DQUOTE 		= '\"'
POUND 		= '#'

PLUS 	= 'PLUS'
MINUS 	= 'MINUS'
MUL 	= 'MUL'
DIV 	= 'DIV'
POWER 	= 'POWER'
EQUALS 	= 'EQUALS'


def is_valid_varchar(c):
	# checks character is a valid character for a variable name: 
	# 	underscores, numbers, letters
	try:
		return c.isalnum() or c == '_'
	except AttributeError:
		return False

def is_valid_varname_firstchar(c):
	# variable names can begin with underscores or letters
	try:
		return c.isalpha() or c == '_'
	except AttributeError:
		return False

def is_valid_numchar(c):
	# numbers can begin with number or decimal
	try:
		return c.isdigit() or c == '.'
	except AttributeError:
		return False