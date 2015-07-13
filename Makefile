CXX = g++

CFLAGS = -lboost_system -lboost_thread

all: bin/COmap

bin/COmap:
	mkdir bin
	$(CXX) --std=c++11 -g -march=native src/*.cpp -o bin/COmap $(CFLAGS) 

clean:
	rm -rf bin
