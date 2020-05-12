# clean compiled python code (.pyc and .pyo files)
clean-pyc:
    # ToDo

# clean all the personal builds for the project including packaging related files
clean-build:
	rm --force --recursive build/
	rm --force --recursive dist/
	rm --force --recursive *.egg-info

# clean all
clean: clean-pyc clean-build

# create a personal build of the project
build: clean
	python setup.py build_ext --build-lib .
