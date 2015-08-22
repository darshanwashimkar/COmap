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
extern int MIN_RREADS;
extern int MIN_CONSENSUS;

using namespace std;

Aligner::Aligner(unsigned int br, std::unordered_map<unsigned int, uint8_t> &rel_read_map){
	this->base_read = br;

	for(auto kv : rel_read_map) {
		this->tar_reads.push_back(kv.first);
	}
	min_index = 0;
	max_index = INT_MAX;
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
	//print();
	printMultiAlignInfo();
	fixIndelErrors(reads, corrected_reads);
}

void Aligner::alignPair(om_read &br, om_read &tr, unsigned int tar_r_no){
	
	std::vector<std::vector<int>> alignment(br.map_read.size());
	std::pair<int, std::vector<int>> ad;

	AdjAlignmentDifference alignment_data;
	alignment_data.diff.resize(br.map_read.size(), 0);

	ad.second.resize(br.map_read.size(), 0);

	om_read rev_tr = tr.reverse();
	scoring_params sp(.2,1.2,.9,7,17.43,0.58, 0.0015, 0.8, 1, 3);
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
	double t_score_thresh = 1;
	double t_mult = 0;
	
//	std::cout<<for_score<<"  " <<rev_score<<std::endl;
//	std::cout<<for_t_score<<"  "<<rev_t_score<<std::endl;

	int b_ptr = 0;

	if(for_score > rev_score && for_t_score > t_score_thresh && for_score > score_thresh){
		
		if(for_alignment.ref_restr_al_sites.size() > 0){
			b_ptr = for_alignment.ref_restr_al_sites[for_alignment.ref_restr_al_sites.size()-1];
		     	ad.first = for_alignment.tar_restr_al_sites[for_alignment.tar_restr_al_sites.size()-1];
			alignment_data.start = for_alignment.tar_restr_al_sites[for_alignment.tar_restr_al_sites.size()-1];

			if(b_ptr > this->min_index)
				min_index = b_ptr;

//			std::cout<<"\n=="<<for_alignment.tar_restr_al_sites[for_alignment.tar_restr_al_sites.size()-1]<<"\n"<<std::endl;
		}

		for(int k=for_alignment.ref_restr_al_sites.size()-1; k>0; k--){
			int ref_diff = for_alignment.ref_restr_al_sites[k-1] - for_alignment.ref_restr_al_sites[k];
			int tar_diff = for_alignment.tar_restr_al_sites[k-1] - for_alignment.tar_restr_al_sites[k];
//			std::cout<<"  - "<< ref_diff <<"  "<<tar_diff<<std::endl;
			while(ref_diff>1){
				ad.second.at(b_ptr) = -1;
				alignment_data.diff.at(b_ptr) = -1;
				b_ptr++;
				ref_diff--;
			}
			ad.second.at(b_ptr) = tar_diff;
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

		if(b_ptr < max_index)
			max_index = b_ptr;

		//cout<<endl<<endl;
		for_alignment.output_alignment(cout);	
	}
	else if(for_score <= rev_score && rev_t_score > t_score_thresh && rev_score > score_thresh ){

//-----------------
		for(int k=rev_alignment.ref_restr_al_sites.size()-1; k>=0; k--){
			if(k!=rev_alignment.ref_restr_al_sites.size()-1)
			std::cout<<" ";
			std::cout<<rev_alignment.ref_restr_al_sites[k];
			std::cout<<" ";
			std::cout<<rev_alignment.tar_restr_al_sites[k];
		}
		cout<<endl<<endl;
		for_alignment.output_alignment(cout);
	}
	else{
		return;
	}

	alignment_data.a_read = tar_r_no;
	multi_align_info.push_back(alignment_data);

	a_diff.push_back(ad);
	aligned_reads.push_back(tar_r_no);
	
	return;
}


void Aligner::print(){
	std::cout<<"\nAligher Object with base read: "<<base_read<<std::endl<<"Target Reads: ";
	for(int i = 0; i < tar_reads.size(); i++){
		std::cout<<tar_reads.at(i);
		if(i != (tar_reads.size() - 1)){
			std::cout<<", ";
		}
	}
	std::cout<<std::endl<<"Aligned reads and adj difference alighment"<<std::endl;

	for(int i = 0; i < aligned_reads.size(); i++){
		std::cout<<aligned_reads.at(i)<<" : "<<a_diff.at(i).first<<" - ";
		for(int j = 0; j <a_diff.at(i).second.size(); j++){
			std::cout<<a_diff.at(i).second.at(j)<<" ";
		}
		std::cout<<std::endl;
	}		
}

void Aligner::printMultiAlignInfo(){	
	cout<<endl<<"------------------------"<<endl;
	cout<<"Base Read: "<<base_read<<endl;	
	cout<<"Target Reads: ";
	for(int i = 0; i < tar_reads.size(); i++)
		cout<<tar_reads.at(i)<<" ";
	cout<<endl;
	cout<<"Printing alignment info"<<endl;

	for(int i = 0; i < multi_align_info.size(); i++){
		cout<<"Base read aligned to: "<<multi_align_info.at(i).a_read<<endl;
		cout<<"Starting at: "<<multi_align_info.at(i).start<<" -> ";
		for(int j = 0; j < multi_align_info.at(i).diff.size(); j++){
			cout<<multi_align_info.at(i).diff.at(j)<<" ";
		}
		cout<<endl;
	}
	
}


void Aligner::fixIndelErrors(std::vector<Read> & reads, std::vector<Read> & corrected_reads){
	std::vector<int> cur_read_ptr;
	std::vector<double> corrected_base_fragments(reads.at(base_read).fragments.begin(), reads.at(base_read).fragments.begin() + min_index);

	for(int i = 0; i < a_diff.size(); i++){
		cur_read_ptr.push_back(a_diff.at(i).first);
	}
	std::cout<<std::endl<<"Min_index: "<<this->min_index<<" Max_index "<<this->max_index<<std::endl;

	std::vector<std::pair<int, int> > max_con_o_t; // varible to keep track of occurrences of fragment overlap
	
	
	for(int b_ptr = min_index; b_ptr < max_index; b_ptr++){
		std::vector< std::vector<int> > con_o(7); // Can get segfault here // [-1, 0, 1, 2, 3, 4, 5] count occurrences
		for(int j = 0; j < a_diff.size(); j++){
//			if(a_diff.at(j).second.at(b_ptr) != -1)
				con_o.at(a_diff.at(j).second.at(b_ptr) + 1).push_back(aligned_reads.at(j));
		}
		
		int max_con_index = -1;
		int max_con = 0;
		for(int k = 0; k < con_o.size(); k++){
			if(con_o.at(k).size() > max_con){
				max_con_index = k;
				max_con = con_o.at(k).size();
			}			
		}
		
		/* Check if alignment satisfies the requirement of minimum consensus between reads */
		/* Don't correct any error in such type of alignemnt */
		/* '-2' value in 'max_con_o_t' means there are no consensus as per user entered min_consensus value */
		if(max_con < MIN_CONSENSUS){
			max_con_o_t.push_back(std::make_pair(-2, max_con));
			continue;
		}

		/*  */

		if((max_con_index - 1) == -1){ // if insertion error
			max_con_o_t.push_back(std::make_pair((max_con_index-1), max_con));
		}

//		cout<<(max_con_index - 1)<<" -> "<<max_con<<endl;
		
	}	
	//cout<<"Size: "<<a_diff.size()<<"and Size: "<<aligned_reads.size();
}
