.PHONY: all

all:
	mkdir -p build; cd build; cmake .. && make  

.PHONY: clean
clean:
	-rm -r build

.PHONY: test
test: 
	cd build; make test

.PHONY: dry
dry: clean all test
