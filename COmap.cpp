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
#include <boost/thread.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <stdint.h>

class KmerReadIndex;

using namespace std;
int bin_s = 300;
int k = 3;
int no_of_threads = 1;
string om_file = "/s/oak/b/nobackup/muggli/goat/whole_genome_mapping/goat_whole_genome.maps";
int NUMBER_OF_BLOCKS = 100;
uint8_t MIN_RREADS = 2;

/* Quantize the value */
void quantize(unsigned int *val, int *bin_size){
    if (*val % *bin_size < *bin_size / 2.0)
        *val = *val - *val % *bin_size;
    else
        *val = *val - *val % *bin_size +  *bin_size;
    return;
}


/* split takes read string and convert it into unsigned int vector after quantizing values */
std::vector<unsigned int> &split(const std::string &s, char delim, std::vector<unsigned int> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
	unsigned int number = (atof(item.c_str())*1000);
	quantize(&number , &bin_s);
	if( number > 0){
		elems.push_back(number);
	}
    }
    return elems;
}

/* Read execution parameters from command line */
int readParameters(int argc, char **argv){

	int c;
	char *pEnd;

	/* Parsing arguments */
	while ((c = getopt (argc, argv, "k:b:f:t:?")) != -1){
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

	      case 't':
		no_of_threads = strtol(optarg, &pEnd, 10);
		break;

	      case '?':
		std::cout<<"Usage: %%COmap [-k Kmer] [-b BinSize] [-f File Name] [-t No of Threads]"<<std::endl;
		return(0);
	      default:
		return(-1);
	      }	
	}
	
	return(1);
}


/* Print parameters for the run */
void printParameters(){
	std::cout<<"==================================================="<<std::endl;
	std::cout<<"Running with following parameters,";
	std::cout<<"\nK : "<<k<<std::endl;
	std::cout<<"Bin size : "<<bin_s<<std::endl;
	std::cout<<"File : "<<om_file<<std::endl;
	std::cout<<"No of threads : "<<no_of_threads<<std::endl;
	std::cout<<"==================================================="<<std::endl;	
}


class KmerReadIndex{

	public:

	/* Data structure to store hashed information of Kmer */
	std::unordered_map<std::string, std::vector<unsigned int> > kmer_map;
	
	void readFileAndCreateIndex(std::ifstream &infile, std::vector<std::vector<unsigned int> > *reads){

		std::string line;
		int line_number = 1;
		
		/* Read each line from file */
		while (std::getline(infile, line)){
			++line_number;
		
			/* Find required line from the file */
			if(line_number%3==0){
				std::vector<unsigned int> elems;
				split(line, '\t' , elems);

				/* remove all short reads - having read_length <= 15 */
				if(elems.size()<=15){
					continue;
				}
				
				/* Save reads for later processing */
				reads->push_back(elems);
		
				/* Create K-mers from read */
				int head = 0;
			
				while(head+k <= elems.size()){
					string kmer = "";
					for(unsigned int i = 0; i < k; i++){
						kmer = kmer + " "+std::to_string(elems.at(head + i));					
					}
				
					/* Create vectore to store read number corresponding to given kmer */
					std::vector<unsigned int> read_no(1, (unsigned int)((line_number/3)- 1));

					/* Check if key already exit if not then insert */				
					std::pair<std::unordered_map<std::string, std::vector<unsigned int> >::iterator, bool> iter;
					iter = kmer_map.insert (pair<std::string,std::vector<unsigned int> >(kmer,read_no) ); 
	 
					/* if key exist then add the read number to the hash map data structure */
					if (iter.second==false) {

					    iter.first->second.push_back((unsigned int)((line_number/3)- 1));
					}
								
					head++;
				}
			}			
		}			
	}

	void printKmerStatastcs(){

		/* To Print statastics */	
		unsigned int total_kmers = 0;
		pair<std::string, std::vector<unsigned int> > me;
		BOOST_FOREACH(me, kmer_map) {
		  total_kmers = total_kmers + me.second.size();
		}

		std::cout<<"Total Kmers: "<<total_kmers<<std::endl;
		std::cout<<"No of Distinct Kmer with K = "<<k<<" are : "<<kmer_map.size()<<std::endl;
		std::cout<<"Average number of reads associated: "<<(double)((double)total_kmers/(double)kmer_map.size())<<std::endl;
	}
};


class RelatedReadsIndex {
	
	public:	
	/* Data Structure to store related reads, size of which is equal to number of  number of reads */
	std::vector<std::unordered_map<unsigned int, uint8_t> > rel_reads;	
	boost::ptr_vector<boost::mutex> mutex_pv;

	

	RelatedReadsIndex(unsigned int no_of_reads){
		rel_reads.resize(no_of_reads);
		for(unsigned int i = 0; i< no_of_reads; i++){
			mutex_pv.push_back(new boost::mutex);
		}
	}

	void buildRelatedReadsIndex(KmerReadIndex &KRI, std::pair<unsigned int,unsigned int> read_block, std::pair<unsigned int,unsigned int> kmer_block){
		/* Build data structure of realated reads */
		std::pair<std::string, std::vector<unsigned int> > pair_to_iter;
		std::pair<std::unordered_map<unsigned int, uint8_t>::iterator, bool> iter;
	
		/* Variable to make thread process part of kmer index*/
		unsigned int count = 0;
		BOOST_FOREACH( pair_to_iter, KRI.kmer_map) {		

			/* Check if count is out of specified range */
			if(count < kmer_block.first){				
				count++;
				continue;
			}	
			
			if(count > kmer_block.second){
				return;
			}

			count++;

			/* If single read associated with kmer */
			if(pair_to_iter.second.size() == 1){
				continue;
			}

			for(unsigned int i =0; i< pair_to_iter.second.size()-1; i++){			
				for(int j=i+1; j< pair_to_iter.second.size();j++){	
					/* (R1-> R5, R5, R8...) | Check such condition */	
					if( pair_to_iter.second.at(i) ==  pair_to_iter.second.at(j)){ continue; }
					
					else if(pair_to_iter.second.at(i) <  pair_to_iter.second.at(j)){
						/* Check indexed read is in slot location */						
						if((pair_to_iter.second.at(i) < read_block.first) || (pair_to_iter.second.at(i) >= read_block.second)){
							continue;
						}

						/* Obtain lock on related reads index (R1->R2,R4) at R1 */
						mutex_pv[pair_to_iter.second.at(i)].lock();

						/* (R1-> R3, R5, R8...) | Insert R1 to> R3 */						
						iter = rel_reads[ pair_to_iter.second.at(i)].insert(std::pair<unsigned int,uint8_t>(pair_to_iter.second.at(j),1));

						/* if key exist then increament count of time relation between R1 -> R3 exists */
						if (iter.second==false) {
							iter.first->second++; 				   
						}

						/* unlock on related reads index (R1->R2,R4) at R1 */
						mutex_pv[pair_to_iter.second.at(i)].unlock();					
					}
		
					else{
						/* Check indexed read is in slot location */						
						if((pair_to_iter.second.at(j) < read_block.first) || (pair_to_iter.second.at(j) >= read_block.second)){
							continue;
						}

						/* Obtain lock on related reads index (R1->R2,R4) at R2 to insert R2->R1 relation */
						mutex_pv[pair_to_iter.second.at(j)].lock();

						/* (R1-> R3, R5, R8...) | Insert R3 to> R1 */						
						iter = rel_reads[pair_to_iter.second.at(j)].insert(std::pair<unsigned int,uint8_t>(pair_to_iter.second.at(i),1));

						/* if key exist then increament count of time relation between R3 -> R1 exists */
						if (iter.second==false) {
							iter.first->second++; 				   
						}

						/* unlock on related reads index (R1->R2,R4) at R2 to after inserting R2->R1 relation */
						mutex_pv[pair_to_iter.second.at(j)].unlock();
					}
				}
			}
	
		}
	}
	
	void trimRelatedReadIndex(std::pair<unsigned int,unsigned int> read_block){		
		
		for(unsigned int i = read_block.first; i < read_block.second; i++){			
			for (auto it = rel_reads.at(i).begin(); it != rel_reads.at(i).end();) {
			   if(it->second <= MIN_RREADS) {
			      it = rel_reads.at(i).erase(it);
			   }
			   else
			      it++;
			}

		}

		
	}	
	
	std::vector< std::pair<unsigned int,unsigned int> > createBlocks(unsigned int total_size, unsigned int no_of_blocks){
		std::vector< std::pair<unsigned int,unsigned int> > blocks;
		unsigned int head;
		unsigned int tail;
		for(unsigned int i = 0; i < no_of_blocks -1 ; i++){
			head = (total_size/no_of_blocks) * i;
			tail = (total_size/no_of_blocks) * (i+1);
			blocks.push_back(std::make_pair(head,tail));
		}
		head = (total_size/no_of_blocks) * (no_of_blocks -1);
		tail = total_size;
		blocks.push_back(std::make_pair(head,tail));
		return(blocks);
	}

	/* Print realated reads in this form 0 -> 8(1)2(1)9(2)3(4)5(8) */
	/*void printRelatedReads(){
		std::cout<<"==================================================="<<std::endl;
		for(int i =0; i<rel_reads.size(); i++){
			cout<<i<<" ->";
			for ( int j = 0; j < rel_reads.at(i).size(); j++ ){
				cout<<" "<< rel_reads.at(i).at(j);
			}
			cout<<endl;
		}
		std::cout<<"==================================================="<<std::endl;
	}*/

	/* Print realated reads in this form 0 -> 8(1)2(1)9(2)3(4)5(8) */
	void printRelatedReads(){
		std::cout<<"==================================================="<<std::endl;
		for(unsigned int i =0; i<rel_reads.size(); i++){
			cout<<i<<" -> ";
			for ( auto it = rel_reads.at(i).begin(); it != rel_reads.at(i).end(); ++it ){
				cout << it->first << "(" << it->second<<")";			
			}
			cout << endl;
		}
		std::cout<<"==================================================="<<std::endl;
	}

	/* Print number of comman kmers between reads */
	void printNumberCommanKmerBetweenReads(){
		for(unsigned int i =0; i<rel_reads.size(); i++){
			for(auto kv : rel_reads[i]) {
				cout<<i<<" "<<(unsigned)kv.second<<endl;
			}
		}
	}
	
	/* Print common kmer */
	void printCommonKmerBetweenReads(unsigned int start, unsigned int end){		
		ofstream outfile(std::to_string(start));		

		std::pair<std::unordered_map<unsigned int, unsigned int>::iterator, bool> iter;

		for(int i = start; i<= end; i++){	


			std::unordered_map<unsigned int, unsigned int> temp_map;				
			for ( int j = 0; j < rel_reads.at(i).size(); j++ ){
				iter = temp_map.insert(std::pair<unsigned int,unsigned int>(rel_reads.at(i).at(j),1));

				if (iter.second==false) {
					iter.first->second++; 				   
				}				
			}
	
			for ( auto it = temp_map.begin(); it != temp_map.end(); ++it ){
				outfile << i << " "<< it->second<<std::endl;			
			}		

		}	

		outfile.close();
	}
		
};



void printReadStatastics(std::vector<std::vector<unsigned int> > &reads){
	std::cout<<"Number of reads: "<<reads.size()<<std::endl;
	double total_read_length = 0;
	for(int i =0; i < reads.size(); i++){
		total_read_length += reads[i].size();
	}
	std::cout<<"Average size of read: "<<(total_read_length/reads.size())<<std::endl;
}

int main (int argc, char **argv) {

	if(readParameters(argc,argv)<1){
		return(0);
	}

	/* Open File */
	std::ifstream infile(om_file);
	if(!infile){
		std::cout<<"File doesn't exist: "<<om_file<<std::endl;
		std::cout<<"Usage: %%COmap [-k Kmer] [-b BinSize] [-f File Name] [-t No of Threads]"<<std::endl;
		return(-1);
	}

	/* Print parameters of run */
	//printParameters();

	/* Data structure to store reads */
	std::vector<std::vector<unsigned int> > reads;
	
	/* Create Kmer - Read index */
	KmerReadIndex KRI;
	KRI.readFileAndCreateIndex(infile,&reads);

//	cout<<"\nDone with mapping Kmers----\n";

	/* Print Total number of distinct Kmers and average number of reads associated with each Kmers */
//	KRI.printKmerStatastcs();
	
	/* Print statastics of the reads of optical mapping data */	
	printReadStatastics(reads);

	/* Create realated read index */
	RelatedReadsIndex RRI(reads.size());

	/* Building related read index */	
	/* Create blocks of related read index */
	std::vector< std::pair<unsigned int,unsigned int> > read_blocks;
	read_blocks = RRI.createBlocks(RRI.rel_reads.size(), NUMBER_OF_BLOCKS);
		
	for(int i =0; i < read_blocks.size(); i++){
		/* Create threads to parallarize building related read index */
		std::vector<boost::thread *> thread_pool;

		std::vector< std::pair<unsigned int, unsigned int> > kmer_blocks;
		kmer_blocks = RRI.createBlocks(KRI.kmer_map.size(), no_of_threads);
		for(int j =0; j < kmer_blocks.size(); j++){
			thread_pool.push_back(new boost::thread(boost::bind(&RelatedReadsIndex::buildRelatedReadsIndex,  &RRI, KRI, read_blocks.at(i), kmer_blocks.at(j))));
		}
		
		/* Join threads */
		for (k = 0; k < no_of_threads; k++){
			thread_pool.at(k)->join();			
		}		
	
		/* To reduce the size of related read data structure */					
		RRI.trimRelatedReadIndex(read_blocks.at(i));
	}


	/* Code to Printing related reads */
	RRI.printNumberCommanKmerBetweenReads();













//-----------------------------------------------------------------------------------------------------

	/* Create threads to parallarize building related read index */
//	std::vector<boost::thread *> thread_pool;
	

	/* Find start and end for traversal in kmer index for each thread, and create thread */
/*	unsigned int start;
	unsigned int end;
	int i;
	for (i = 0; i < no_of_threads-1; ++i){ 
	    start = (i*KRI.kmer_map.size())/no_of_threads;
	    end = (((i+1)*KRI.kmer_map.size())/no_of_threads)-1;
	    thread_pool.push_back(new boost::thread(boost::bind(&RelatedReadsIndex::buildRelatedReadsIndex, &RRI, KRI.kmer_map, start, end)));
		
	}
	start = (i*KRI.kmer_map.size())/no_of_threads;
	end = KRI.kmer_map.size();
	thread_pool.push_back(new boost::thread(boost::bind(&RelatedReadsIndex::buildRelatedReadsIndex, &RRI, KRI.kmer_map, start, end)));
	
	for (i = 0; i < no_of_threads; ++i){
		thread_pool.at(i)->join();
		delete thread_pool.at(i);
	}

	
	/* Code to Printing related reads */
//	RRI.printRelatedReads();
	
	/* Code to print number of common kmers between reads, R1 <-> R2. Multithreading printing values in different files*/
/*	thread_pool.erase(thread_pool.begin(),thread_pool.end());

	for (i = 0; i < no_of_threads-1; ++i){
		start = (i*RRI.rel_reads.size())/no_of_threads;
		end = (((i+1)*RRI.rel_reads.size())/no_of_threads)-1;
		thread_pool.push_back(new boost::thread(boost::bind(&RelatedReadsIndex::printCommonKmerBetweenReads, &RRI, start, end)));
	}
	start = (i*RRI.rel_reads.size())/no_of_threads;
	end = RRI.rel_reads.size() - 1;
	thread_pool.push_back(new boost::thread(boost::bind(&RelatedReadsIndex::printCommonKmerBetweenReads, &RRI, start, end)));

	for (i = 0; i < no_of_threads; ++i){
		thread_pool.at(i)->join();
		delete thread_pool.at(i);
	}


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
	
*/
	infile.close();
	return (0);
}

























