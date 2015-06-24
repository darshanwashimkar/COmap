/*
 *  	Author: Darshan Washimkar
 *	About:  COmap is a program that correct error in optical mapping data
 */

#include "COmap.hpp"
#include "read.hpp"
#include "KRI.hpp"
#include "RRI.hpp"

using namespace std;
int BIN_S = 300;
int K = 3;
int NO_OF_THREADS = 1;
string OM_FILE = "/s/oak/b/nobackup/muggli/goat/whole_genome_mapping/goat_whole_genome.maps";
int NUMBER_OF_BLOCKS = 20;
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
std::vector<unsigned int> &split(const std::string &s, char delim, std::vector<unsigned int> &elems, std::string &enzyme, std::string &something) {
    std::stringstream ss(s);
    std::string item;
    int i = 0;
    while (std::getline(ss, item, delim)) {	
	/* To get enzyme name */
	if(i==1){
		enzyme = item;
		i++;
		continue;
	}
	/* To get something string of a read */
	if(i==2){
		something = item;
		i++;
		continue;
	}

	unsigned int number = (atof(item.c_str())*1000);
	quantize(&number , &BIN_S);
	if( number > 0){
		elems.push_back(number);
	}

	i++;
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
		K = strtol(optarg, &pEnd, 10);
		if (K<=0){
			std::cout<<"Please enter integer value greater than 0 for Kmer"<<std::endl;
			return(-1);
		}			
		break;

	      case 'b':
		BIN_S = strtol(optarg, &pEnd, 10);
		if (BIN_S<=0){
			std::cout<<"Please enter integer value greater than 0 for Bin Size"<<std::endl;
			return(-1);
		}
		break;
	      
	      case 'f':
		OM_FILE = optarg;
		break;

	      case 't':
		NO_OF_THREADS = strtol(optarg, &pEnd, 10);
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
	std::cout<<"\nK : "<<K<<std::endl;
	std::cout<<"Bin size : "<<BIN_S<<std::endl;
	std::cout<<"File : "<<OM_FILE<<std::endl;
	std::cout<<"No of threads : "<<NO_OF_THREADS<<std::endl;
	std::cout<<"==================================================="<<std::endl;	
}


void printReadStatastics(std::vector<Read> &reads){
	std::cout<<"Number of reads: "<<reads.size()<<std::endl;
	double total_read_length = 0;
	for(int i =0; i < reads.size(); i++){
		total_read_length += reads[i].fragments.size();
		cout<<reads[i].fragments.size()<<"  ";
		cout<<reads[i].name<<" ";
		cout<<reads[i].enzyme<<" "<<endl;
	}
	std::cout<<"Average size of read: "<<(total_read_length/reads.size())<<std::endl;
}




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
	RRI.printNumberCommanKmerBetweenReads();


	infile.close();
	return (0);
}





















