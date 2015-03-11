CXX = g++

all: COmap

COmap:
	$(CXX) --std=c++11 COmap.cpp -o COmap

clean:
	rm -rf COmap
