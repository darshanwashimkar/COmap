/*
 *      Author: Darshan Washimkar
 *      About:  Alignment related data structure and functions
 *
 */

#include "align.h"

using namespace std;

void Aligner::alignSet(RelatedReadsIndex &RRI){

	double readarry1[] = {4.123, 5.123, 2.315, 4.152, 9.432, 11.212, 8.271, 7.456, 4.564, 9.179, 34.473, 2.682, 13.115, 6.171, 4.170, 22.259, 3.319, 3.281, 7.171, 3.170, 12.259, 9.319, 12.181};

	double readarry2[] = {11.212, 8.271, 7.456, 4.564, 9.179, 34.473, 2.682, 13.115, 6.171, 4.170, 22.259, 3.319, 3.281, 7.171, 3.170, 12.259, 9.319, 12.181, 15.213, 5.122, 11.123, 13.212, 4.123, 5.112};
	
	om_read read1, read2;
	read1.read_name = "5038535_0_1317";
	read1.Enz_name = "SpeI";
	read1.Enz_acr = "1317";
	read1.map_read.assign(readarry1, readarry1+(sizeof(readarry1)/sizeof(readarry1[0])));

	read2.read_name = "5038579_0_14";
	read2.Enz_name = "SpeI";
	read2.Enz_acr = "1317";
	read2.map_read.assign(readarry2, readarry2+(sizeof(readarry2)/sizeof(readarry2[0])));

	scoring_params sp(.2,1.2,.9,7,17.43,0.58, 0.0015, 0.8, 1, 3);
	
	om_read tar_map = read1; 
	om_read for_map = read2; 
	om_read rev_map = for_map.reverse();

	rm_alignment for_alignment(tar_map, for_map, sp);
	rm_alignment rev_alignment(tar_map, rev_map, sp);


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


	if(for_score > rev_score && for_t_score > t_score_thresh && for_score > score_thresh){
		for(int k=for_alignment.ref_restr_al_sites.size()-1; k>=0; k--){
			if(k!=for_alignment.ref_restr_al_sites.size()-1)
			std::cout<<" ";
			std::cout<<for_alignment.ref_restr_al_sites[k];
			std::cout<<" ";
			std::cout<<for_alignment.tar_restr_al_sites[k];
		}
		cout<<endl<<endl;
		for_alignment.output_alignment(cout);	
	}
	if(for_score <= rev_score && rev_t_score > t_score_thresh && rev_score > score_thresh ){
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
}

