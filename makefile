# Build eBACOLI library and examples. Run tests on 'standard' examples

.PHONY: all debug opt src examples test opttest debugtest clean
all: src examples

debug:
	@echo "############################################################################"
	@echo "## Warning: may need to make clean to get debug symbols."
	@echo "############################################################################"
	cd src && make debug
	cd examples && make debug
	@echo "############################################################################"
	@echo "## Warning: may need to make clean to get debug symbols."
	@echo "############################################################################"

opt:
	@echo "############################################################################"
	@echo "## Warning: may need to make clean to get debug symbols."
	@echo "############################################################################"
	cd src && make opt
	cd examples && make opt
	@echo "############################################################################"
	@echo "## Warning: may need to make clean to get debug symbols."
	@echo "############################################################################"

src:
	cd src && make all

examples:
	cd examples && make all

test: src
	cd examples && time make test

opttest: opt
	cd examples && time make test

debugtest: debug
	cd examples && time make test

clean:
	cd src && make clean
	cd examples && make clean
