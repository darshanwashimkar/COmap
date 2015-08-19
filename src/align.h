/*
 *      Author: Darshan Washimkar
 *      About:  Alignment related data structure and functions
 *
 */
#define _ALIGN_H
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <math.h>
#include <climits>
#include "./om_set1/alignment.h"
#include "./om_set1/msfl.h"

#ifndef _READ_HPP
#include "read.hpp"
#endif

class Aligner{
	unsigned int base_read;
	std::vector<unsigned int> tar_reads;

	/* Valuev doesn't produce alignment all the pairs hence need to store aligned reads */
	std::vector<unsigned int> aligned_reads;

	/* pair of "start of alignment in base read" and "difference of adjacent alignemnts" (-1, -1, 2)*/
	std::vector<std::pair<int, std::vector<int> > > a_diff; 

	/* Keep track of minimum and maximux index in a set of corresponding reads */
	int min_index, max_index;

	std::vector<std::vector< std::vector<int> > > a_metrix;

    public:
	Aligner(unsigned int, std::unordered_map<unsigned int, uint8_t> &);

	/* alignSet() aligns targeted reads with base reads and save alignment in a_matrix (alignment matrix) */
	void alignSet(std::vector<Read> &, std::vector<Read> &);

	/* alignPair() aligns the pair of base read and target read */
	std::vector< std::vector<int>> alignPair(om_read &, om_read &, unsigned int);
	
	/* Helper function to create om_read object from the read object this may not required in later stage*/
	void createOMRead(om_read &,Read &);

	/* This method corrects "Indel" errors from aligned set of reads based on output of valuev */
	void fixIndelErrors(std::vector<Read> &, std::vector<Read> &);

	/* Helper function to print content of aligner object */
	void print();
};

