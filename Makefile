.PHONY: clean-pyc clean-build docs clean

help:
	@echo "clean - remove all build, test, coverage and Python artifacts"
	@echo "clean-pyc - remove Python file artifacts"
	@echo "clean-test - remove test and coverage artifacts"
	@echo "docs - generate Sphinx HTML documentation, including API docs"
	@echo "check - check all the dependencies"

clean: clean-pyc clean-test

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +
clean-test:
	cd test && bash run_test.sh -c
docs:
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	# open docs/_build/html/index.html
check:
	@echo "Not implemented :("
