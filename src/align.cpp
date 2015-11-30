/*
 *      Author: Darshan Washimkar
 *      About:  Alignment related data structure and functions
 *
 */

#include "align.h"

extern int BIN_S;
extern int K;
extern int NO_OF_THREADS;
extern std::string OM_FILE;
extern int NUMBER_OF_BLOCKS;
extern int MIN_COMMON_K_IN_READS;
extern int MIN_CONSENSUS;
extern unsigned int count_updated_r;

using namespace std;

Aligner::Aligner(unsigned int br, std::unordered_map<unsigned int, uint8_t> &rel_read_map){
	this->base_read = br;

	for(auto kv : rel_read_map) {
		this->tar_reads.push_back(kv.first);
	}
	//min_index = 0;
	//max_index = INT_MAX;
}

void Aligner::createOMRead(om_read &om_r, Read & r){
	om_r.read_name = r.name;
	om_r.Enz_name = r.enzyme;
	om_r.Enz_acr = r.something;
	om_r.map_read = r.fragments;
}

void Aligner::alignSet(std::vector<Read> & reads, std::vector<Read> & corrected_reads){
	om_read br;
	createOMRead(br, reads.at(this->base_read));
	
	for(int i = 0 ; i < this->tar_reads.size(); i++){
		om_read tr;
		createOMRead(tr, reads.at(this->tar_reads.at(i)));
		alignPair(br,tr, this->tar_reads.at(i));
	}

	printMultiAlignInfo(reads);

	/* Check if we have minimux number of reads to form consensus */
	if(multi_align_info.size() >=  MIN_CONSENSUS){	
		fixIndelErrors(reads, corrected_reads);
	}
}

void Aligner::alignPair(om_read &br, om_read &tr, unsigned int tar_r_no){
	
	std::vector<std::vector<int>> alignment(br.map_read.size());

	AdjAlignmentDifference alignment_data;
	alignment_data.diff.resize(br.map_read.size(), 0);

	om_read rev_tr = tr.reverse();
	scoring_params sp(0.2,1.2,.9, 7,17.43,0.58, 0.01, 0.75, 1, 3);
	rm_alignment for_alignment(br, tr, sp);
	rm_alignment rev_alignment(br, rev_tr, sp);


	for_alignment.optimized_overlap_alignment();
	rev_alignment.optimized_overlap_alignment();

	for_alignment.overlap_t_score();
	rev_alignment.overlap_t_score();

	double for_score = for_alignment.Smax;
	double rev_score = rev_alignment.Smax;

	double for_t_score = for_alignment.Tmax;
	double rev_t_score = rev_alignment.Tmax;

	double score_thresh = 25;
	double t_score_thresh = 0;
	double t_mult = 0;
	
//	std::cout<<"fs: \t"<<for_score<<"  \trs: \t" <<rev_score<<std::endl;
//	std::cout<<"fs_t: \t"<<for_t_score<<"  \t\trs_t: \t"<<rev_t_score<<std::endl;

	int b_ptr = 0;

	if(for_score > rev_score && for_t_score > t_score_thresh && for_score > score_thresh){
		if(for_alignment.ref_restr_al_sites.size() > 0){
			b_ptr = for_alignment.ref_restr_al_sites[for_alignment.ref_restr_al_sites.size()-1];		     	
			alignment_data.start = for_alignment.tar_restr_al_sites[for_alignment.tar_restr_al_sites.size()-1];

/*			if(b_ptr > this->min_index)
				min_index = b_ptr;
*/
//			std::cout<<"\n=="<<for_alignment.tar_restr_al_sites[for_alignment.tar_restr_al_sites.size()-1]<<"\n"<<std::endl;
		}

		for(int k=for_alignment.ref_restr_al_sites.size()-1; k>0; k--){
			int ref_diff = for_alignment.ref_restr_al_sites[k-1] - for_alignment.ref_restr_al_sites[k];
			int tar_diff = for_alignment.tar_restr_al_sites[k-1] - for_alignment.tar_restr_al_sites[k];
//			std::cout<<"  - "<< ref_diff <<"  "<<tar_diff<<std::endl;
			while(ref_diff>1){
				alignment_data.diff.at(b_ptr) = -1;
				b_ptr++;
				ref_diff--;
			}

			alignment_data.diff.at(b_ptr) = tar_diff;
			b_ptr++;

//-----------------
/*			if(k!=for_alignment.ref_restr_al_sites.size()-1)	
			std::cout<<" ";
			std::cout<<for_alignment.ref_restr_al_sites[k];
			std::cout<<" ";
			std::cout<<for_alignment.tar_restr_al_sites[k];
*/
		}

/*		if(b_ptr < max_index)
			max_index = b_ptr;
*/
		//cout<<endl<<endl;
		for_alignment.output_alignment(cout);	
	}
	else if(for_score <= rev_score && rev_t_score > t_score_thresh && rev_score > score_thresh ){

		return; /* Need to update this */
//-----------------
		for(int k=rev_alignment.ref_restr_al_sites.size()-1; k>=0; k--){
			if(k!=rev_alignment.ref_restr_al_sites.size()-1)
			std::cout<<" ";
			std::cout<<rev_alignment.ref_restr_al_sites[k];
			std::cout<<" ";
			std::cout<<rev_alignment.tar_restr_al_sites[k];
		}
		cout<<endl<<endl;
		//for_alignment.output_alignment(cout);
	}
	else{
		return;
	}

	alignment_data.a_read = tar_r_no;
	multi_align_info.push_back(alignment_data);
	
	return;
}


void Aligner::printMultiAlignInfo(std::vector<Read> &reads){	
	cout<<endl<<"------------------------"<<endl;
	cout<<"Base Read: "<<base_read<<" "<<reads.at(base_read).name<<endl;	
	cout<<"Target Reads: ";
	for(int i = 0; i < tar_reads.size(); i++)
		cout<<tar_reads.at(i)<<" ";
	cout<<endl;
	cout<<"Printing alignment info"<<endl;

	for(int i = 0; i < multi_align_info.size(); i++){
		cout<<"Base read aligned to: "<<multi_align_info.at(i).a_read<<" "<< reads.at(multi_align_info.at(i).a_read).name<<endl;
		cout<<"Starting at: "<<multi_align_info.at(i).start<<" -> ";
		for(int j = 0; j < multi_align_info.at(i).diff.size(); j++){
			cout<<multi_align_info.at(i).diff.at(j)<<" ";
		}
		cout<<endl;
	}
	
}


void Aligner::fixIndelErrors(std::vector<Read> & reads, std::vector<Read> & corrected_reads){
	std::vector<int> cur_read_ptr;
	std::vector<double> corrected_base_frag;//(reads.at(base_read).fragments.begin(), reads.at(base_read).fragments.begin() + min_index);
	std::vector<double> temp_deleted_frag;

//	std::cout<<std::endl<<"Min_index: "<<this->min_index<<" Max_index "<<this->max_index<<std::endl;

	std::vector<std::pair<int, int> > consensus; // varible to keep track of maximux number of fragment overlap, pair<consensus, count>
	bool consensus_happening = true;

	/* Initialize 'start' fragment according for each target read as per generated alignment to alignment */
/*	for(int i = 0; i < min_index; i++){
		for(int j = 0; j < multi_align_info.size(); j++){
			if(multi_align_info.at(j).diff.at(i) != -1)
				multi_align_info.at(j).start += multi_align_info.at(j).diff.at(i);
		}
	}
*/
//	printMultiAlignInfo(reads);
	int deletion_corrected = 0;
	int insertion_error = 0;
	int no_of_minus_one = 0;


	/* not processing first and last fragment */
	for(int b_ptr = 0; b_ptr < reads.at(base_read).fragments.size() - 1; b_ptr++){
		std::vector< std::vector<int> > con_o(10); // Can get segfault here // [-1, 0, 1, 2, 3, 4, 5] count occurrences
		for(int j = 0; j < multi_align_info.size(); j++){
			// because 0 means overlap
			if(multi_align_info.at(j).diff.at(b_ptr) != 0){
				con_o.at(multi_align_info.at(j).diff.at(b_ptr) + 1).push_back(j);
			}
		}

		int max_count_index = -1;
		int max_count = 0;
		for(int k = 0; k < con_o.size(); k++){
			if(con_o.at(k).size() > max_count){
				max_count_index = k;
				max_count = con_o.at(k).size();
			}
		}
		

		/* Check if alignment satisfies the requirement of minimum consensus between reads */
		/* Don't correct any error in such type of alignemnt */
		/* '-2' value in 'consensus' means there are no consensus as per user entered min_consensus value */		
		if(max_count < MIN_CONSENSUS){

			if(temp_deleted_frag.size() != 0){
				std::copy(temp_deleted_frag.begin(), temp_deleted_frag.end(), std::back_inserter(corrected_base_frag));
				temp_deleted_frag.clear();
			}

			consensus.push_back(std::make_pair(-2, max_count));
			corrected_base_frag.push_back(reads.at(base_read).fragments.at(b_ptr));	
			consensus_happening = false;		
			no_of_minus_one = 0;
		}

		// if insertion error
		else if((max_count_index - 1) == -1){ 
			consensus.push_back(std::make_pair((max_count_index-1), max_count));		
			temp_deleted_frag.push_back(reads.at(base_read).fragments.at(b_ptr));
			consensus_happening = true;			
			no_of_minus_one++;
		}

		else {			
			temp_deleted_frag.clear();

			/* find average for only [(max_count_index-1) > 1] (one fragment aligning to more than one fragment) */
			if(consensus_happening && (max_count_index-1) > 1){

				/* Finding average of fragments forming consusus */
				for(int m = 0; m < (max_count_index-1); m++){

					double add = 0.0;
					AdjAlignmentDifference *t_info; // pointer to the alignment_info to make code easy to understand

					for(int n = 0; n < con_o.at(max_count_index).size(); n++){
						t_info = &multi_align_info.at(con_o.at(max_count_index).at(n));
						add += reads.at(t_info->a_read).fragments.at(t_info->start + m); 
						// m is added because multiple fragments can align to single fragment from base read. 
						// So such case we are no updating the start hence we need to add m.					
					}
					
					/* Adding average to as corrected read */
					corrected_base_frag.push_back(floor((add/con_o.at(max_count_index).size()) * 1000) / 1000);
					
				}
			}
			else{
				consensus_happening = true;
				corrected_base_frag.push_back(reads.at(base_read).fragments.at(b_ptr));
			}

			if(max_count_index > 1){
				if((max_count_index - no_of_minus_one) > 2){
					deletion_corrected += (max_count_index - no_of_minus_one - 2);					
				}
				else if((max_count_index - no_of_minus_one) < 2){
					insertion_error += -(max_count_index - no_of_minus_one - 2);
				}
				no_of_minus_one = 0;
			}
			consensus.push_back(std::make_pair((max_count_index-1), max_count));
		}

		/* Update start of every alignment info */
		for(int j = 0; j < multi_align_info.size(); j++){
			if(multi_align_info.at(j).diff.at(b_ptr) != -1)
				multi_align_info.at(j).start += multi_align_info.at(j).diff.at(b_ptr);
		}
	}
	
	/* Copy if anything is remaining in temp_deleted_frag */
	if(temp_deleted_frag.size() > 0){
		std::copy(temp_deleted_frag.begin(), temp_deleted_frag.end(), std::back_inserter(corrected_base_frag));		
	}
	/* Copy last fragment */
	corrected_base_frag.push_back(reads.at(base_read).fragments.at(reads.at(base_read).fragments.size() - 1));


	/* Copy if any remaining fragment-lengths */
/*	for(int k = this->max_index; k < reads.at(base_read).fragments.size(); k++){		
		corrected_base_frag.push_back( reads.at(base_read).fragments.at(k) );
	}
*/

	/* Update the corrected Read */
	reads.at(base_read).fragments.swap(corrected_base_frag);

	count_updated_r++; // delete this
	static int number_of_reads_with_more_than_five_alignment = 0;
	if(multi_align_info.size() >= 4){
		number_of_reads_with_more_than_five_alignment++;
	}
//	std::cout<<"-> "<<number_of_reads_with_more_than_five_alignment<<endl;
	if(deletion_corrected > 0 || insertion_error > 0){
		std::cout<<base_read<<" == "<<deletion_corrected<<"  ";
		std::cout<<base_read<<" == "<<insertion_error<<std::endl;
	}

        /* print corrected reads */
	std::cout<<std::endl;
	for(int z = 0; z < reads.at(base_read).fragments.size(); z++){
		
		cout<<reads.at(base_read).fragments.at(z)<<"\t";
	}
	cout<<endl;



	cout<<"-*-*-*-*-*-*"<<endl;
	// Printing consensus
	//cout<< "start from : "<<min_index<<endl;
	for(int z = 0; z < consensus.size(); z++ ){
		cout<<consensus.at(z).first<<" ";
	}
	cout<<endl;

	for(int z = 0; z < consensus.size(); z++ ){
		cout<<consensus.at(z).second<<" ";
	}
	cout<<" - "<<consensus.size();

	cout<<endl<<"-*-*-*-*-*-*"<<endl;


}

