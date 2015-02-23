/*
 *  	Author: Darshan Washimkar
 *	About:  COmap is a program that correct error in optical mapping data
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <stdlib.h> 
#include <unordered_map>
#include <boost/foreach.hpp>

using namespace std;
int bin_s = 300;
int k = 3;

/* Data structure to store hashed information */


/* Quantize the value */
void quantize(long int *val, int *bin_size){
    if (*val % *bin_size < *bin_size / 2.0)
        *val = *val - *val % *bin_size;
    else
        *val = *val - *val % *bin_size +  *bin_size;
    return;
}


/* split takes read string and convert it into long int vector after quantizing values */
std::vector<long int> &split(const std::string &s, char delim, std::vector<long int> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
	long int number = (atof(item.c_str())*1000);
	quantize(&number , &bin_s);
	if( number > 0){
		elems.push_back(number);
	}
    }
    return elems;
}

int main () {
	std::unordered_map<std::string, std::vector<long int> > kmer_map;
	std::vector<std::vector<long int> > reads;
	std::string line;
	std::ifstream infile("temp.txt");
	int line_number = 1;

	
	/* Read each line from file */
	while (std::getline(infile, line)){
		++line_number;
		
		/* Find required line from the file */
		if(line_number%3==0){
			std::vector<long int> elems;
			split(line, '\t' , elems);

			/* Save reads for later processing */
			reads.push_back(elems);
		
			/* Create K-mers from read */
			int head = 0;
			
			while(head+k <= elems.size()){
				string kmer = "";
				for(int i = 0; i < k; i++){
					kmer = kmer + " "+std::to_string(elems.at(head + i));					
				}
				
				/* Create vectore to store read number corresponding to given kmer */
				std::vector<long int> read_no(1, (long int)((line_number/3)- 1));

				/* Check if key already exit if not then insert */				
				std::pair<std::unordered_map<std::string, std::vector<long int> >::iterator, bool> iter;
				iter = kmer_map.insert (pair<std::string,std::vector<long int> >(kmer,read_no) ); 
 
				/* if key exist then add the read number to the hash map data structure */
				if (iter.second==false) {

				    iter.first->second.push_back((long int)((line_number/3)- 1));
				}
								
				head++;
			}
		}			
	}
	
	pair<std::string, std::vector<long int> > me; // what a map<int, int> is made of
	BOOST_FOREACH(me, kmer_map) {
	  cout << me.first<<"  :  ";
	  for(int i =0; i<me.second.size(); i++){
		  cout << me.second[i]<<"  ";
	  }
	  cout <<"\n";
	}


	infile.close();
	return 0;
}

























