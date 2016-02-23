/*
 *  	Author: Darshan Washimkar
 *	About:  COmap is a program that correct error in optical mapping data
 */
#include <iostream>
#include <vector>
#include <unordered_map>
#include "RRI.hpp"


using namespace std;


int main (int argc, char **argv) {

	if(readParameters(argc,argv)<1){
		return(0);
	}

	/* Open File */
	std::ifstream infile(OM_FILE);
	if(!infile){
		std::cout<<"File doesn't exist: "<<OM_FILE<<std::endl;
		std::cout<<"Usage: %%COmap [-k Kmer] [-b BinSize] [-f File Name] [-t No of Threads]"<<std::endl;
		return(-1);
	}

	/* Print parameters of run */
	//printParameters();

	/* Data structure to store reads */
	std::vector<Read> reads;
	
	/* Create Kmer - Read index */
	KmerReadIndex KRI;
	KRI.readFileAndCreateIndex(infile, reads);


	/* Print Total number of distinct Kmers and average number of reads associated with each Kmers */
//	KRI.printKmerStatastcs();
	
	/* Print statastics of the reads of optical mapping data */	
//	printReadStatastics(reads);
//	cout<<"---------------------"<<endl;
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
		kmer_blocks = RRI.createBlocks(KRI.kmer_map.size(), NO_OF_THREADS);
		
		if(NO_OF_THREADS == 1){
			RRI.buildRelatedReadsIndex(KRI, read_blocks.at(i), kmer_blocks.at(0));						
		}
		else{
			for(int j =0; j < NO_OF_THREADS; j++){
				thread_pool.push_back(new boost::thread(boost::bind(&RelatedReadsIndex::buildRelatedReadsIndex,  &RRI, KRI, read_blocks.at(i), kmer_blocks.at(j))));
			}
			
			/* Join threads */
			for (int k = 0; k < NO_OF_THREADS; k++){
				thread_pool.at(k)->join();			
			}		
		}
		/* To reduce the size of related read data structure */					
		RRI.trimRelatedReadIndex(read_blocks.at(i));
	}


	/* Code to Printing related reads */
//	RRI.printNumberCommanKmerBetweenReads();
	if(debug){
		RRI.printRelatedReads();
	}

	/* Store corrected reads */
	/* Create Copy for reads */
	std::vector<Read> corrected_reads(reads);

	RRI.correctReads(reads, corrected_reads);

	if(p_corrected_r){
		/* Print reads */
		printReads(reads);
	}

	infile.close();
	return (0);
}





















