/*
 *      Author: Darshan Washimkar
 *      About:  Alignment related data structure and functions
 *
 */

#include "align.h"

using namespace std;

Aligner::Aligner(unsigned int br, std::unordered_map<unsigned int, uint8_t> &rel_read_map){
	this->base_read = br;

	for(auto kv : rel_read_map) {
		this->tar_reads.push_back(kv.first);
	}
	
}

void Aligner::createOMRead(om_read &om_r, Read & r){
	om_r.read_name = r.name;
	om_r.Enz_name = r.enzyme;
	om_r.Enz_acr = r.something;
	om_r.map_read = r.fragments;
}

void Aligner::alignSet(std::vector<Read> & reads){
	om_read br;
	createOMRead(br, reads.at(this->base_read));
	
	for(int i = 0 ; i < this->tar_reads.size(); i++){
		om_read tr;
		createOMRead(tr, reads.at(this->tar_reads.at(i)));
		alignPair(br,tr, this->tar_reads.at(i));
	}
}

std::vector<std::vector<int>> Aligner::alignPair(om_read &br, om_read &tr, unsigned int tar_r_no){
	
	std::vector<std::vector<int>> alignment(br.map_read.size());
	std::pair<int, std::vector<int>> ad;
	ad.second.resize(br.map_read.size(), 0);

	om_read rev_tr = tr.reverse();
	scoring_params sp(.2,1.2,.9,7,17.43,0.58, 0.0015, 0.8, 2, 1);
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
	double t_score_thresh = 8;
	double t_mult = 0;
	
	std::cout<<"sssssssss"<<std::endl;
	std::cout<<for_score<<"  " <<rev_score<<std::endl;
	std::cout<<for_t_score<<"  "<<rev_t_score<<std::endl;

	int b_ptr = 0;

	if(for_score > rev_score && for_t_score > t_score_thresh && for_score > score_thresh){
		
		if(for_alignment.ref_restr_al_sites.size() > 0){
			ad.first = for_alignment.ref_restr_al_sites[for_alignment.ref_restr_al_sites.size()-1];
			b_ptr = ad.first;
			std::cout<<"\n=="<<ad.first<<"\n"<<std::endl;
		}

		for(int k=for_alignment.ref_restr_al_sites.size()-1; k>0; k--){
			int ref_diff = for_alignment.ref_restr_al_sites[k-1] - for_alignment.ref_restr_al_sites[k];
			int tar_diff = for_alignment.tar_restr_al_sites[k-1] - for_alignment.tar_restr_al_sites[k];
			std::cout<<"  - "<< ref_diff <<"  "<<tar_diff<<std::endl;
			while(ref_diff>1){
				ad.second.at(b_ptr) = -1;
				b_ptr++;
				ref_diff--;
			}
			ad.second.at(b_ptr) = tar_diff;			
			b_ptr++;

//-----------------
			if(k!=for_alignment.ref_restr_al_sites.size()-1)	
			std::cout<<" ";
			std::cout<<for_alignment.ref_restr_al_sites[k];
			std::cout<<" ";
			std::cout<<for_alignment.tar_restr_al_sites[k];
		}
		cout<<endl<<endl;
		for_alignment.output_alignment(cout);	
	}
	else if(for_score <= rev_score && rev_t_score > t_score_thresh && rev_score > score_thresh ){
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
		return(alignment);
	}

	a_diff.push_back(ad);
	aligned_reads.push_back(tar_r_no);
	print();
	return(alignment);
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
