.PHONY: env-init clean jupyter

VIRTUALENV=venv-${APP}
PIP=${VIRTUALENV}/bin/pip
PYTHON=PYTHONPATH=$(pwd)/src ${VIRTUALENV}/bin/python3

env-init:
	python3 -m venv ${VIRTUALENV}
	${PIP} install --upgrade pip setuptools wheel
	${PIP} install -r requirements.txt

clean:
	@echo "-> cleaned build environment"
	@rm -rf $(VIRTUALENV) || true

jupyter:
	PYTHONPATH=`pwd`/src: ${VIRTUALENV}/bin/jupyter notebook