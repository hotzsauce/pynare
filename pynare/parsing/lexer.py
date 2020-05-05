""" 
module for the pynare Lexer - reads in source code and generates a stream of 
tokens reflecting the Dynare language. A pdf of the manual can be found here:
https://www.dynare.org/manual.pdf
"""

import pynare.parsing.base as base

import pynare.parsing.model as mdl
import pynare.parsing.variables as vbl 
import pynare.parsing.functions as fnc 

from pynare.parsing.tokens import Token


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

	def error(self, char):
		msg = 'Invalid character: {}'.format(repr(char))
		raise Exception(msg)

	def advance(self, n=1):
		""" advance the position and set the current character """
		self.pos += n
		try:
			self.current_char = self.text[self.pos]
		except IndexError:
			self.current_char = None

	def peek(self, n=1):
		""" returns the next character w/o incrementing past the current one """
		peek_pos = self.pos + n
		try:
			return self.text[peek_pos]
		except IndexError:
			return None

	def skip_whitespace(self):
		while self.current_char is not None and self.current_char.isspace():
			self.advance()

	def skip_comment(self):
		""" single line comments are signified by '%' """
		self.advance()
		while self.current_char is not '\n':
			self.advance()
		self.advance()
		self.skip_whitespace() # just in case

	def skip_multiline_comment(self):
		""" multiline comments begin '/*' and end '*/' """
		self.advance(n=2)
		while (self.current_char != '*') and (self.peek() != '/'):
			self.advance()
		self.advance(n=2)
		self.skip_whitespace() # just in case

	def _id(self):
		""" 
		recognizing characters or sequences of characters that may or may not be 
		recognized keywords. If they are not, we tag them as a Token of type ID,
		and assume they're a variable
		"""
		result = ''
		while base.is_valid_varchar(self.current_char):
			result += self.current_char
			self.advance()
		token = RESERVED_KEYWORDS.get(result, Token(base.ID, result))
		return token

	def _number(self):
		"""
		recognizing numbers. We treat everything as a float, although Dynare (and 
		Matlab underneath it) treats ints and floats as distinct types 
		"""
		result = ''
		while (
			self.current_char is not None 
			and base.is_valid_numchar(self.current_char)
		):
			result += self.current_char
			self.advance()
		num = float(result)
		return Token(base.NUMBER, num)

	def _string(self, quote):
		"""
		recognizing literal strings, which can be enclosed by single or double 
		quotes in Dynare/Matlab.
		"""
		result = ''
		self.advance()
		while self.current_char is not quote:
			result += self.current_char
			self.advance()
		self.advance()
		return Token(base.STRING, result)

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
				return self._number()

			if self.current_char == '\'':
				return self._string(quote='\'')

			if self.current_char == '\"':
				return self._string(quote='\"')

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