""" module for recognizing, parsing, and interpreting different vtypes in pynare """

from pynare.parsing.tokens import Token

# DIFFERENT VARIABLE TYPES AS DEFINED IN DYNARE DOCUMENTATION
VAR 			= 'var'
VAREXO 			= 'varexo'
VAREXO_DET 		= 'varexo_det'
PARAMETERS 		= 'parameters'
PREDETERMINED 	= 'predetermined_variables'
TREND_VAR 		= 'trend_var'
LOG_TREND_VAR 	= 'log_trend_var'
CHANGE_TYPE 	= 'change_type'


RESERVED_KEYWORDS = {
	VAR: 			Token(VAR, VAR),
	VAREXO: 		Token(VAREXO, VAREXO),
	VAREXO_DET: 	Token(VAREXO_DET, VAREXO_DET),
	PARAMETERS: 	Token(PARAMETERS, PARAMETERS),
	PREDETERMINED: 	Token(PREDETERMINED, PREDETERMINED),
	TREND_VAR: 		Token(TREND_VAR, TREND_VAR),
	LOG_TREND_VAR: 	Token(LOG_TREND_VAR, LOG_TREND_VAR),
	CHANGE_TYPE: 	Token(CHANGE_TYPE, CHANGE_TYPE)
}

REQUIRED_DECLARATIONS = {VAR}

OPTIONAL_DECLARATIONS = {
	VAREXO,
	VAREXO_DET,
	PARAMETERS,
	PREDETERMINED,
	TREND_VAR,
	LOG_TREND_VAR,
	CHANGE_TYPE
}