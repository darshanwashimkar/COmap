#include <boost/foreach.hpp>
#include <boost/thread.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "KRI.hpp"

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

	void buildRelatedReadsIndex(KmerReadIndex &KRI, std::pair<unsigned int,unsigned int> read_block, std::pair<unsigned int,unsigned int> kmer_block);
	void trimRelatedReadIndex(std::pair<unsigned int,unsigned int> read_block);
	std::vector< std::pair<unsigned int,unsigned int> > createBlocks(unsigned int total_size, unsigned int no_of_blocks);
	void printRelatedReads();
	void printNumberCommanKmerBetweenReads();
	void printCommonKmerBetweenReads(unsigned int start, unsigned int end);
		
};
