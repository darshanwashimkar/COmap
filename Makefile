CXX = g++

CFLAGS = -lboost_system -lboost_thread

all: COmap

COmap:
	$(CXX) --std=c++11 $(CFLAGS) -g COmap.cpp -o COmap

clean:
	rm -rf COmap
