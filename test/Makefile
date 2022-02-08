.PHONY: all

all:
	mkdir -p build; cd build; cmake .. && make  

.PHONY: clean
clean:
	-rm -r build

.PHONY: test
test: 
	cd build; make test

.PHONY: doc
doc: 
	doxygen

.PHONY: dry
dry: clean all test

.PHONY: deepclean
deepclean:
	-rm -r build doc

