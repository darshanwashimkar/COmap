#define _UTIL_HPP

#include <iostream>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <sstream>
#include <math.h>
#include "read.hpp"


void split(const std::string &s, char delim, std::vector<unsigned int> &elems, Read & read);
int readParameters(int, char **);
void printParameters();
void printReadStatastics(std::vector<Read> &);
void printReads(std::vector<Read> &);


extern bool p_error_count;
extern bool debug;
extern bool p_corrected_r;

extern int BIN_S;
extern int K;
extern int NO_OF_THREADS;
extern std::string OM_FILE;
extern int NUMBER_OF_BLOCKS;
extern int MIN_COMMON_K_IN_READS;
extern int MIN_CONSENSUS;
