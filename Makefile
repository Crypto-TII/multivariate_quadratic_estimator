PACKAGE=mq_estimator
SAGE_BIN=sage

all: install

build:
	$(SAGE_BIN) -pip install build
	$(SAGE_BIN) -python -m build

uninstall:
	$(SAGE_BIN) -pip uninstall $(PACKAGE) -y

install: build
	$(SAGE_BIN) -pip install -r requirements.txt
	$(SAGE_BIN) -pip install .

test: install
	$(SAGE_BIN) -t -T 600 src/mpkc  # timeout is set to 600 seconds

doc: install
	cd docs/ && $(SAGE_BIN) -sh -c "make html"

clean-doc:
	cd docs/ && $(SAGE_BIN) -sh -c "make clean"

clean: clean-doc
	rm -rf build/
	rm -rf dist/
	rm -rf scripts/*.sage.py
	rm -rf src/$(PACKAGE).egg-info/
