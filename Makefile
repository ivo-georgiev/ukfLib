.PHONY: all

all:
	mkdir -p build; cd build; cmake .. && make  

.PHONY: clean
clean:
	rm -r build
