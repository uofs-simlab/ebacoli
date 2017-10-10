# Build examples.
#  Testing done solely through standard examples
#  - silence after executable name means it passed

.PHONY: all standard extended debug opt test opttest debugtest clean
all: standard extended

debug:
	cd standard && make debug
	cd extended && make debug

standard:
	cd standard && make

extended:
	cd extended && make

opt:
	cd standard && make opt
	cd extended && make opt

test:
	cd standard && make test

opttest: opt
	cd standard && make test

debugtest: debug
	cd standard && make test

clean:
	cd standard && make clean
	cd extended && make clean
