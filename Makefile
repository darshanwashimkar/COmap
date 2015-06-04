CXX = g++

CFLAGS = -lboost_system -lboost_thread

all: bin/COmap

bin/COmap:
	mkdir bin
	$(CXX) --std=c++11 $(CFLAGS) -g src/COmap.cpp -o bin/COmap

clean:
	rm -rf bin
