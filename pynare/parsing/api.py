"""
forward-facing functions for parsing the text that defines models
"""

from __future__ import annotations

from pynare.parsing.lexer import DynareLexer
from pynare.parsing.parser import DynareParser

from pynare.parsing.factory import ModelFactory


def read_model(
	model_text: str,
	language: str
) -> AST:
	"""
	The top-level function for translating a model from text -> AST. Called in 
	the Model classmethod 'from_path'


	Parameters
	----------
	model_text : str
		a list of the lines of the text that defines a model
	language : str
		the language the model is written is	

	Returns
	-------
	Abstract Syntax Tree of the model's definition
	"""
	return ModelFactory(model_text, language)
