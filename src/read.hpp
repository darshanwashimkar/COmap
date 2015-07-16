#define _READ_HPP

#include <iostream>
#include <vector>

class Read{

	public:	
	std::vector<double> fragments;
	std::string name;
	std::string enzyme;
	std::string something;	

	void printRead();	
};
