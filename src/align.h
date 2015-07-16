/*
 *      Author: Darshan Washimkar
 *      About:  Alignment related data structure and functions
 *
 */
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <math.h>
#include "./om_set1/alignment.h"
#include "./om_set1/msfl.h"

#ifndef _READ_HPP
#include "read.hpp"
#endif

class Aligner{
	unsigned int base_read;
	std::vector<unsigned int> tar_reads;
	std::vector<std::vector<unsigned int>> a_metrix;
    public:
	Aligner(unsigned int, std::unordered_map<unsigned int, uint8_t> &);

	/* alignSet() aligns targeted reads with base reads and save alignment in a_matrix (alignment matrix) */
	void alignSet(std::vector<Read> &);

	/* alignPair() aligns the pair of base read and target read */
	void alignPair(om_read &, om_read &);
	
	/* Helper function to create om_read object from the read object this may not required in later stage*/
	void createOMRead(om_read &,Read &);
};

