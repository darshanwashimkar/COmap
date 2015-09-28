#include "RRI.hpp"

extern int MIN_COMMON_K_IN_READS;
unsigned int count_updated_r = 0;


using namespace std;

void RelatedReadsIndex::buildRelatedReadsIndex(KmerReadIndex &KRI, std::pair<unsigned int,unsigned int> read_block, std::pair<unsigned int,unsigned int> kmer_block){
	
	/* If nothing to process */
	if(read_block.first == read_block.second){
		return;
	}

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

		unsigned int prev_read;
		for(unsigned int i =0; i< pair_to_iter.second.size()-1; i++){

			/* If same read have multiple k-mers in comman */
			if((i != 0) && prev_read == pair_to_iter.second.at(i)){
				continue;
			}
			prev_read = pair_to_iter.second.at(i);

			for(unsigned int j=i+1; j< pair_to_iter.second.size();j++){	


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
						std::cout<<pair_to_iter.second.size()<<" - ";
						if(pair_to_iter.second.size() == 2)
							std::cout<<pair_to_iter.second.at(0)<<" "<<pair_to_iter.second.at(1);
						std::cout<<std::endl;
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
	
	//temp Code
	//cout<<"Done with itteration"<<endl;
}






void RelatedReadsIndex::trimRelatedReadIndex(std::pair<unsigned int,unsigned int> read_block){		
	
	for(unsigned int i = read_block.first; i < read_block.second; i++){		

		if(rel_reads.at(i).empty()){
		   continue;
		}

		for (auto it = rel_reads.at(i).begin(); it != rel_reads.at(i).end();) {
		   if(it->second < MIN_COMMON_K_IN_READS) {
		      it = rel_reads.at(i).erase(it);
		   }
		   else
		      it++;
		}
	}
}


std::vector< std::pair<unsigned int,unsigned int> > RelatedReadsIndex::createBlocks(unsigned int total_size, unsigned int no_of_blocks){
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
void RelatedReadsIndex::printRelatedReads(){
	std::cout<<"==================================================="<<std::endl;
	for(unsigned int i =0; i<rel_reads.size(); i++){
		cout<<i<<" -> ";
		for ( auto it = rel_reads.at(i).begin(); it != rel_reads.at(i).end(); ++it ){
			cout << it->first << "(" << (unsigned)it->second<<")";			
		}
		cout <<" size:"<<rel_reads.at(i).size()<< endl;
	}
	std::cout<<"==================================================="<<std::endl;
}


/* Print number of comman kmers between reads */
void RelatedReadsIndex::printNumberCommanKmerBetweenReads(){
	for(unsigned int i =0; i<rel_reads.size(); i++){
		for(auto kv : rel_reads[i]) {
			cout<<i<<" "<<(unsigned)kv.second<<endl;
		}
	}
}

void RelatedReadsIndex::correctReads(std::vector<Read> & reads, std::vector<Read> & corrected_reads){
	for(unsigned int i =0; i<rel_reads.size(); i++){
		if(rel_reads[i].size()){
			Aligner aligner(i, rel_reads[i]);
			aligner.alignSet(reads, corrected_reads);			
		}
	}
//	std::cout<<"Count of updated reads:"<<count_updated_r<<std::endl;
}


/* Print common kmer */
void RelatedReadsIndex::printCommonKmerBetweenReads(unsigned int start, unsigned int end){		
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



