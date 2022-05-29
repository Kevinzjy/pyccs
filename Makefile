.PHONY: help all activate clean develop lib install bin uninstall test

all: install

help:
	@echo "make install"
	@echo "         install pyccs"
	@echo "make clean"
	@echo "         clean previous build"
	@echo "make test"
	@echo "         test pyccs"
	@echo "make uninstall"
	@echo "         uninstall pyccs"

activate:
	. venv/bin/activate
	unset CONDA_PREFIX

clean:
	cargo clean

develop:	
	pip install maturin
	maturin develop

lib:
	pip install maturin
	maturin build --interpreter python

install:
	pip install ./target/wheels/pyccs-*-cp39-cp39-manylinux_2_28_x86_64.whl --force-reinstall

bin:
	maturin build --interpreter python -b bin

uninstall:
	pip uninstall pyccs

test:
	pip install pytest
	pytest
	# cargo test -- --nocapture

deploy_test:
	maturin publish -r https://test.pypi.org/legacy/