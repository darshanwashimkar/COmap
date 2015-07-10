#include "read.hpp"

using namespace std;

void Read::printRead(){
	cout<<this->name<<"\t"<<this->enzyme<<"\t"<<this->something<<endl;
	for(int i = 0; i < fragments.size(); i++){
		cout<<"\t"<<fragments.at(i);
	}
}
