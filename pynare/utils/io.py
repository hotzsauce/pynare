
from __future__ import annotations

from typing import Union

from pynare.examples import ExampleRegistry


def read_file(
	filepath: Union[str, Path] 
) -> str:
	"""
	Reads in a file line by line and returns those lines in a list

	Parameters
	----------
	filepath : str | Path
		a string or a pathlib.Path to the file

	Returns
	-------
	the file as a string
	"""
	file = open(filepath, 'r')
	file_str = file.read()
	file.close()
	return file_str


def read_example(
	ex_name: str
):
	"""
	Reads in an example model located in the pynare/examples directory

	Parameters
	----------
	ex_name : str
		the example file's name


	Returns
	-------
	list of lines in the file
	"""
	example_path = ExampleRegistry.get_example(ex_name)
	file = open(example_path, 'r')
	file_str = file.read()
	file.close()
	return file_str
