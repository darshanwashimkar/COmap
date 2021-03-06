#include "KRI.hpp"

extern int K;
extern int MIN_FRAG_IN_READ;

using namespace std;

void KmerReadIndex::readFileAndCreateIndex(std::ifstream &infile, std::vector<Read> &reads){

	std::string line;
	int line_number = 1;
	std::string read_name = "";

	/* Read each line from file */
	while (std::getline(infile, line)){
		++line_number;
	
		if((line_number+1)%3==0){
			read_name = line;
			continue;
		}

		/* Find required line from the file */
		if(line_number%3==0){
			Read temp_read;
			std::vector<unsigned int> qfrags;
			split(line, '\t' , qfrags, temp_read);

			/* remove all short reads - having read_length < 10 */
			if(qfrags.size() < MIN_FRAG_IN_READ){
				continue;
			}				
			temp_read.name = read_name;

			/* Save reads for later processing */
			reads.push_back(temp_read);
	
			/* Create K-mers from read in forward direction */
			int head = 0;
		
			while(head+K <= qfrags.size()){
				string kmer = "";
				for(unsigned int i = 0; i < K; i++){
					kmer = kmer + " "+std::to_string(qfrags.at(head + i));					
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

			/* Create K-mers from read in reverse direction */
			head = qfrags.size() - 1;
			
			while(head-(K-1) >= 0){
				string kmer = "";
				for(unsigned int i = 0; i < K; i++){
					kmer = kmer + " "+std::to_string(qfrags.at(head - i));
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
							
				head--;
			}
		}			
	}			
}


void KmerReadIndex::printKmerStatastcs(){

	/* To Print statastics */	
	unsigned int total_kmers = 0;
	pair<std::string, std::vector<unsigned int> > me;
	BOOST_FOREACH(me, kmer_map) {
	  total_kmers = total_kmers + me.second.size();
	  std::cout<<me.first<<std::endl;
	}

	std::cout<<"Total Kmers: "<<total_kmers<<std::endl;
	std::cout<<"No of Distinct Kmer with K = "<<K<<" are : "<<kmer_map.size()<<std::endl;
	std::cout<<"Average number of reads associated: "<<(double)((double)total_kmers/(double)kmer_map.size())<<std::endl;
}
