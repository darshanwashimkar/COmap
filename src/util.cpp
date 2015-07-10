#include "util.hpp"

int BIN_S = 300;
int K = 3;
int NO_OF_THREADS = 1;
std::string OM_FILE = "/s/oak/b/nobackup/muggli/goat/whole_genome_mapping/goat_whole_genome.maps";
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
void split(const std::string &s, char delim, std::vector<unsigned int> &elems, Read &read) {
    std::stringstream ss(s);
    std::string item;
    int i = 0;
    while (std::getline(ss, item, delim)) {	
	/* To get enzyme name */
	if(i==1){
		read.enzyme = item;
		i++;
		continue;
	}
	/* To get something string of a read */
	if(i==2){
		read.something = item;
		i++;
		continue;
	}
	
	double frag = atof(item.c_str());
	unsigned int number = ((frag)*1000);
	quantize(&number , &BIN_S);
	if( number > 0){
		elems.push_back(number);
		read.fragments.push_back(frag);
	}

	i++;
    }
    return;
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






