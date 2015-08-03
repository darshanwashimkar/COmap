#include <iostream>
#include <vector>
#include <unordered_map>
#include <boost/foreach.hpp>

#ifndef _UTIL_HPP
#include "util.hpp"
#endif

class KmerReadIndex{

	public:

	/* Data structure to store hashed information of Kmer */
	std::unordered_map<std::string, std::vector<unsigned int> > kmer_map;
	

	void readFileAndCreateIndex(std::ifstream &infile, std::vector<Read> &reads);
	void printKmerStatastcs();

};
