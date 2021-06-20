PACKAGE=mq_estimator
SAGE_BIN=sage

all: install

build:
	$(SAGE_BIN) -python -m build

uninstall:
	$(SAGE_BIN) -pip uninstall $(PACKAGE) -y

install: build
	$(SAGE_BIN) -pip install .

test: install
	$(SAGE_BIN) -t src/mpkc

doc: install
	cd docs/ && $(SAGE_BIN) -sh -c "make html"

clean-doc:
	cd docs/ && $(SAGE_BIN) -sh -c "make clean"

clean: clean-doc
	rm -rf build/
	rm -rf dist/
	rm -rf src/$(PACKAGE).egg-info/
