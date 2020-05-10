""" module for tokenizing the different simulation commands """

from pynare.parsing.tokens import Token 


INITVAL 	= 'initval'
ENDVAL 		= 'endval'
HISTVAL 	= 'histval'
STEADY 		= 'steady'
RESID 		= 'resid'
INITVALFILE = 'initval_file'
HISTVALFILE = 'histval_file'

RESERVED_KEYWORDS = {
	INITVAL: 		Token(INITVAL, INITVAL),
	ENDVAL: 		Token(ENDVAL, ENDVAL),
	HISTVAL:		Token(HISTVAL, HISTVAL),
	STEADY: 		Token(STEADY, STEADY),
	RESID: 			Token(RESID, RESID),
	INITVALFILE: 	Token(INITVALFILE, INITVALFILE),
	HISTVALFILE: 	Token(HISTVALFILE, HISTVALFILE)
}

CONDITION_TYPES = {
	INITVAL, 
	ENDVAL, 
	HISTVAL
}