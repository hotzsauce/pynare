""" module for tokenizing the different simulation commands """

from pynare.parsing.tokens import Token 


INITVAL 	= 'initval'
ENDVAL 		= 'endval'
HISTVAL 	= 'histval'
RESID 		= 'resid'
INITVALFILE = 'initval_file'
HISTVALFILE = 'histval_file'

SHOCKS 	= 'shocks'
PERIODS = 'periods'
VALUES 	= 'values'
CORR 	= 'corr'
STDERR 	= 'stderr'

RESERVED_KEYWORDS = {
	INITVAL: 		Token(INITVAL, INITVAL),
	ENDVAL: 		Token(ENDVAL, ENDVAL),
	HISTVAL:		Token(HISTVAL, HISTVAL),
	STEADY: 		Token(STEADY, STEADY),
	RESID: 			Token(RESID, RESID),
	INITVALFILE: 	Token(INITVALFILE, INITVALFILE),
	HISTVALFILE: 	Token(HISTVALFILE, HISTVALFILE),
	SHOCKS: 		Token(SHOCKS, SHOCKS),
	PERIODS:		Token(PERIODS, PERIODS),
	VALUES:			Token(VALUES, VALUES),
	CORR:			Token(CORR, CORR),
	STDERR:			Token(STDERR, STDERR)
}

CONDITION_TYPES = {
	INITVAL, 
	ENDVAL, 
	HISTVAL
}