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

using namespace std;
int bin_s = 300;
int k = 3;

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
			reads.push_back(elems);
		
			/* Create K-mers from read */
			int head = 0;
			
			while(head+k <= elems.size()){
				string kmer = "";
				for(int i = 0; i < k; i++){
					kmer = kmer + " "+std::to_string(elems.at(head + i));
					
				}
				cout<<kmer<<endl;
				head++;
			}

		}			
	}
	cout<<reads.at(1).size()<<endl;
	infile.close();
	return 0;
}

























