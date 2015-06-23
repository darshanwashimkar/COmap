#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

class Read{

	public:	
	std::vector<unsigned int> fragments;
	std::string name;
	std::string enzyme;
	std::string something;	
};
