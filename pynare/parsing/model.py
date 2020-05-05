""" module for recognizing, parsing, and interpreting the model block of a modfile """

from pynare.parsing.tokens import Token

MODEL 	= 'model'
END 	= 'end'

RESERVED_KEYWORDS = {
	MODEL: Token(MODEL, MODEL),
	END: Token(END, END)
}