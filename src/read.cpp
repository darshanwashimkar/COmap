#include "read.hpp"

using namespace std;

void Read::printRead(){
	cout<<this->name<<std::endl<<"\t"<<this->enzyme<<"\t"<<this->something;
	for(int i = 0; i < fragments.size(); i++){
		cout<<"\t"<<fragments.at(i);
	}
}
