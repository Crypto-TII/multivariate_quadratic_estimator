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
	$(SAGE_BIN) -t src/$(PACKAGE)

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf src/$(PACKAGE).egg-info/
