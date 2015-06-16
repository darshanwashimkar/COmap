#include <iostream>
#include <fstream>
#include <string>

#include <vector>
#include <algorithm>
#include <math.h>

using namespace std;

#include "./../om_set1/msfl.cpp"
#include "./../om_set1/m_read.cpp"
#include "./../om_set1/scoring.cpp"
#include "./../om_set1/alignment.cpp"

int main(int argc, char *argv[])    
{
    double readarry1[] = {8.271, 2.456, 1.564, 9.179, 34.473, 2.682, 13.115, 6.171, 4.170, 22.259, 3.319, 2.181, 7.171, 3.170, 12.259, 9.319, 12.181};

    double readarry2[] = {11.212, 8.271, 2.456, 1.564, 9.179, 34.473, 2.682, 13.115, 6.171, 4.170, 22.259, 3.319, 2.181, 7.171, 3.170, 12.259, 9.319, 12.181};
    om_read read1, read2;
    read1.read_name = "5038535_0_1317";
    read1.Enz_name = "SpeI";
    read1.Enz_acr = "1317";
    read1.map_read.assign(readarry1, readarry1+(sizeof(readarry1)/sizeof(readarry1[0])));

    read2.read_name = "5038579_0_14";
    read2.Enz_name = "SpeI";
    read2.Enz_acr = "1317";
    read2.map_read.assign(readarry2, readarry2+(sizeof(readarry2)/sizeof(readarry2[0])));
    
    ifstream ifs;
   
    
    om_read_collection maps(ifs);
   



    scoring_params sp(.2,1.2,.9,7,17.43,0.58, 0.0015, 0.8, 1, 3);

	  om_read tar_map = read1; //maps.collection[i];
	  om_read for_map = read2; //maps.collection[j];
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

	  double for_ovlp_size = for_alignment.ovlp_size();
	  double rev_ovlp_size = rev_alignment.ovlp_size();

	  cout<<"fs: "<<for_score<<" ft: "<<for_t_score;
	  //cout<<" p_v: "<<for_p_value;
	  cout<<endl;
	  cout<<"rs: "<<rev_score<<" rt: "<<rev_t_score;
	  //cout<<" p_v: "<<rev_p_value;
	  cout<<endl;

	  //rev_alignment.output_alignment(cout);

	  double score_thresh = 25;
	  double t_score_thresh = 8;
	  double t_mult = 0;

	  if(for_score > rev_score){
	    //for_alignment.output_alignment(cout);
	  }
	  else{
	    //rev_alignment.output_alignment(cout);
	  }

	  if(for_score > rev_score &&
	     for_t_score > t_score_thresh &&
	     for_score > score_thresh){
	    //&& for_t_score > t_mult*for_ovlp_size){

	    int ovlp_start1 = 
	      for_alignment.ref_restr_al_sites
	      [for_alignment.ref_restr_al_sites.size()-1];
	    int ovlp_end1 = for_alignment.ref_restr_al_sites[0];

	    int ovlp_start2 = 
	      for_alignment.tar_restr_al_sites
	      [for_alignment.tar_restr_al_sites.size()-1];
	    int ovlp_end2 = for_alignment.tar_restr_al_sites[0];

	    cout<<tar_map.read_name.c_str();	    
	    cout<<" "<<for_map.read_name.c_str();
	    cout<<" "<<tar_map.map_read.size();
	    cout<<" "<<for_map.map_read.size();
	    cout<<" 1 1 "<<for_score;
	    cout<<" "<<for_t_score<<endl;//" ";
	    for(int k=for_alignment.ref_restr_al_sites.size()-1; k>=0; k--){
	      if(k!=for_alignment.ref_restr_al_sites.size()-1)
		cout<<" ";
	      cout<<for_alignment.ref_restr_al_sites[k];
	      cout<<" ";
	      cout<<for_alignment.tar_restr_al_sites[k];
	    }

	    cout<<endl<<endl;

	    //if(for_score>score_thresh)
	    //if(for_t_score > t_score_thresh)
	    for_alignment.output_alignment(cout);
	  }
	  if(for_score <= rev_score && 
	     rev_t_score > t_score_thresh &&
	     rev_score > score_thresh ){
	    //&& rev_t_score > t_mult*rev_ovlp_size){

	    int rev_map_size = rev_map.map_read.size();
	    int ovlp_start1 = 
	      rev_alignment.ref_restr_al_sites
	      [rev_alignment.ref_restr_al_sites.size()-1];
	    int ovlp_end1 = rev_alignment.ref_restr_al_sites[0];
	    
	    int ovlp_start2 = 
	      rev_alignment.tar_restr_al_sites
	      [rev_alignment.tar_restr_al_sites.size()-1];
	    int ovlp_end2 = rev_alignment.tar_restr_al_sites[0];

	    assert(ovlp_start2 >=0 && ovlp_end2 >= 0);


	    cout<<tar_map.read_name.c_str();	    
	    cout<<" "<<rev_map.read_name.c_str();
	    cout<<" "<<tar_map.map_read.size();
	    cout<<" "<<rev_map.map_read.size();
	    cout<<" 1 0 "<<rev_score;
	    cout<<" "<<rev_t_score<<endl;//" ";
	    for(int k=rev_alignment.ref_restr_al_sites.size()-1; k>=0; k--){
	      if(k!=rev_alignment.ref_restr_al_sites.size()-1)
		cout<<" ";
	      cout<<rev_alignment.ref_restr_al_sites[k];
	      cout<<" ";
	      cout<<rev_alignment.tar_restr_al_sites[k];
	    }

	    cout<<endl<<endl;
	    	    
	    rev_alignment.output_alignment(cout);
	  }

  

  return 0;  
}
  