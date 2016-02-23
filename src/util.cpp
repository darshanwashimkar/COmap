#define _UTIL_CPP
#include "util.hpp"

bool p_error_count = false;
bool debug = false;
bool p_corrected_r = false;

int BIN_S = 300;
int K = 3;
int NO_OF_THREADS = 1;
std::string OM_FILE = "/s/oak/b/nobackup/muggli/goat/whole_genome_mapping/goat_whole_genome.maps";
int NUMBER_OF_BLOCKS = 20;
int MIN_COMMON_K_IN_READS = 3;
int MIN_CONSENSUS = 4;
double S_VARIENCE = 0.33;

using namespace std;


/* Quantize the value  // OLD
void quantize(unsigned int *val, int *bin_size){
    if (*val % *bin_size < *bin_size / 2.0)
        *val = *val - *val % *bin_size;
    else
        *val = *val - *val % *bin_size +  *bin_size;
    return;
}

*/
/*
unsigned int nextRange(unsigned int cur_range, unsigned int bin_size,double s_varience){
    if(cur_range < 2000)
        return(cur_range + bin_size);
    else
        return(round(((double)(cur_range/1000.0) + 3 * sqrt(((double)cur_range/1000.0) * s_varience))*1000));
}


void quantize(unsigned int *val, unsigned int bin_size, double s_varience){
	static vector<double> quantizeList;
	if(*val == 0){
		return;
	}

	if(quantizeList.size() == 0){
		quantizeList.push_back(0.0);
		quantizeList.push_back(nextRange(0.0, bin_size, s_varience));
	}

	for(int i = 0; i < quantizeList.size(); i++){ 
		if(*val < quantizeList[i]){
			*val = quantizeList[i-1];
			return;
		}
	}

	while(1){
		quantizeList.push_back(nextRange(quantizeList[quantizeList.size() - 1], bin_size, s_varience));
		if(*val < quantizeList[quantizeList.size() - 1]){
			*val = quantizeList[quantizeList.size() - 2];
			return;
		}
	}
}
*/

/* standard deviation with Bionano data */
double signmaS(double len){
    double singma = (0.2*0.2) + len * 0.1 * (-0.1) + len*len*0.04*0.04;
    return(sqrt(singma));
}

unsigned int nextRange(unsigned int pre_val, int folds){
    double S = signmaS(float(pre_val/1000.0));
    double d = (2.0 * folds * S) * 1000;
    return(int(pre_val + d));
}


void quantize(unsigned int *val, int folds, unsigned int start){
    static vector<double> quantizeList;
    
    if(*val <= start)
        return;
        
    if(quantizeList.size() == 0){
        quantizeList.push_back(start);
        quantizeList.push_back(nextRange(start, folds));
    }
    
    for(int i = 0; i < quantizeList.size(); i++){
    	if(*val < quantizeList[i]){        
            *val = quantizeList[i-1];
            return;
    	}
    }
        
    while(1){
        quantizeList.push_back(nextRange(quantizeList[quantizeList.size() - 1], folds));
        if(*val < quantizeList[quantizeList.size() - 1]){
            *val = quantizeList[quantizeList.size() - 2];
            return;
        }
    }
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
	quantize(&number , 4, 500);
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
	while ((c = getopt (argc, argv, ":k:b:f:t:m:d:xyz:?")) != -1){
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

	      case 'm':
			MIN_COMMON_K_IN_READS = strtol(optarg, &pEnd, 10);
			break;		  

		  case 'x':
			p_error_count = true;
			break;			

		  case 'y':
			debug = true;
			break;			

		  case 'z':
			p_corrected_r = true;
			break;			

		  case 'd':
			MIN_CONSENSUS = strtol(optarg, &pEnd, 10);
			break;

	      case '?':
			std::cout<<"Usage: %%COmap [-k Kmer] [-b BinSize] [-f File Name] [-t No of Threads] [-m Min Common k between reads] [-d Number of similar alignment to creat CONSENSUS] [-x enable printing error count] [-y enable printing debug] [-z enable printing corrected reads]"<<std::endl;
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
	std::cout<<"NUMBER_OF_BLOCKS : "<<NUMBER_OF_BLOCKS<<std::endl;
	std::cout<<"MIN_COMMON_K_IN_READS : "<<MIN_COMMON_K_IN_READS<<std::endl;
	std::cout<<"MIN_CONSENSUS : "<<MIN_CONSENSUS<<std::endl;
	std::cout<<"p_error_count : "<<p_error_count<<std::endl;
	std::cout<<"debug : "<<debug<<std::endl;
	std::cout<<"p_corrected_r : "<<p_corrected_r<<std::endl;
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

void printReads(std::vector<Read> & reads){
	cout<<std::setprecision(3)<<std::fixed;
	for(int i = 0; i < reads.size(); i++){
		reads.at(i).printRead();
		cout<<endl<<endl;
	}
}




