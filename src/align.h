/*
 *      Author: Darshan Washimkar
 *      About:  Alignment related data structure and functions
 *
 */
#include <iostream>
#include <fstream>
#include <string>

#include <vector>
#include <algorithm>
#include <math.h>

#include "RRI.hpp"
#include "./om_set1/alignment.h"
#include "./om_set1/msfl.h"
//#include "./om_set1/m_read.h"
//#include "./om_set1/scoring.h"

class Aligner{
	unsigned int base_read;
	std::vector<unsigned int> tar_reads;
	std::vector<std::vector<unsigned int>> a_metrix;
    public:
	void alignSet(RelatedReadsIndex &);
};

