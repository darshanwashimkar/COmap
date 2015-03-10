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
string om_file = "/s/oak/b/nobackup/muggli/goat/whole_genome_mapping/goat_whole_genome.maps";

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

int main (int argc, char **argv) {
	int c;
	char *pEnd;

	/* Data structure to store hashed information of Kmer */
	std::unordered_map<std::string, std::vector<long int> > kmer_map;
	std::vector<std::vector<long int> > reads;
	std::string line;
	int line_number = 1;
	
	/* Parsing arguments */
	while ((c = getopt (argc, argv, "k:b:f:?")) != -1){
	    switch (c)
	      {
	      case 'k':
		k = strtol(optarg, &pEnd, 10);
		if (k<=0){
			std::cout<<"Please enter integer value greater than 0 for Kmer"<<std::endl;
			return(-1);
		}			
		break;

	      case 'b':
		bin_s = strtol(optarg, &pEnd, 10);
		if (bin_s<=0){
			std::cout<<"Please enter integer value greater than 0 for Bin Size"<<std::endl;
			return(-1);
		}
		break;
	      
	      case 'f':
		om_file = optarg;
		break;

	      case '?':
		std::cout<<"Usage: %%COmap [-k Kmer] [-b BinSize] [-f File-Name]"<<std::endl;
		return(0);
	      default:
		return(-1);
	      }	
	}

	/* Open File */
	std::ifstream infile(om_file);
	if(!infile){
		std::cout<<"File doesn't exist: "<<om_file<<std::endl;
		std::cout<<"Usage: %%COmap [-k Kmer] [-b Bin-Size] [-f File-Name]"<<std::endl;
		return(-1);
	}

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

	/* Data Structure to store related reads, size of which is equal to number of  number of reads */
	std::vector<std::unordered_map<long int, int> > rel_reads(reads.size());
	
	/* Build data structure of realated reads */
	std::pair<std::string, std::vector<long int> > pair_to_iter;
	std::pair<std::unordered_map<long int, int>::iterator, bool> iter;
	BOOST_FOREACH( pair_to_iter, kmer_map) {		
		for(int i =0; i< pair_to_iter.second.size()-1; i++){			
			for(int j=i+1; j< pair_to_iter.second.size();j++){	
				/* (R1-> R5, R5, R8...) | Check such condition */	
				if( pair_to_iter.second.at(i) ==  pair_to_iter.second.at(j)){ continue; }

				/* (R1-> R3, R5, R8...) | Insert R1 to> R3 */						
				iter = rel_reads[ pair_to_iter.second.at(i)].insert(std::pair<long int,int>(pair_to_iter.second.at(j),1));

				/* if key exist then increament count of time relation between R1 -> R3 exists */
				if (iter.second==false) {
					iter.first->second++; 				   
				}
				
				/* (R1-> R3, R5, R8...) | Insert R3 to> R1 */						
				iter = rel_reads[ pair_to_iter.second.at(j)].insert(std::pair<long int,int>(pair_to_iter.second.at(i),1));

				/* if key exist then increament count of time relation between R3 -> R1 exists */
				if (iter.second==false) {
					iter.first->second++; 				   
				}
			}
		}
	}	
	
	/* Printing related reads */
	cout<<rel_reads.size()<<endl;
	for(int i =0; i<rel_reads.size(); i++){
		cout<<i<<" -> ";
		for ( auto it = rel_reads.at(i).begin(); it != rel_reads.at(i).end(); ++it ){
			cout << it->first << "(" << it->second<<")";			
		}
		cout << endl;
	}


/* To Print statastic */	

/*        long int total_kmers = 0;
	pair<std::string, std::vector<long int> > me;
	BOOST_FOREACH(me, kmer_map) {
	  total_kmers = total_kmers + me.second.size();
	}

        std::cout<<"Total Kmers: "<<total_kmers<<std::endl;
	std::cout<<"No of Kmer with K = "<<k<<" are : "<<kmer_map.size()<<std::endl;
	std::cout<<"Average number of reads associated: "<<(double)((double)total_kmers/(double)kmer_map.size())<<std::endl;
*/	


//	int count = 1;
//	pair<std::string, std::vector<long int> > me; // what a map<int, int> is made of
//	BOOST_FOREACH(me, kmer_map) {
//	  cout<<count<<" "<<me.second.size()<<endl;
//	  cout << me.first<<"  :  ";
//	  for(int i =0; i<me.second.size(); i++){
//		  cout << me.second[i]<<"  ";		
//	  }
//		count++;
//	}
	

	infile.close();
	return 0;
}

























