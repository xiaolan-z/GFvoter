#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<algorithm>
#include<vector>
#include <utility>
#include<map>
using namespace std;
ofstream out;
map<pair<string,string>, vector<string> > fusion_info_map;
void load_records_and_merge(char*file1,char*file2)
{
    ifstream in(file1);
    ifstream in_(file2);
    string s;
    istringstream istr;
    while(getline(in,s))
    {
	istr.str(s);
	string g1,g2,temp;
	istr>>temp>>g1>>temp>>temp>>temp
	    >>temp>>temp>>temp>>temp>>g2;
	istr.clear();
	pair<string,string> P = make_pair(g1,g2);

	if(fusion_info_map.find(P) ==fusion_info_map.end())
	{
		vector<string> v(1,s);
		fusion_info_map[P] = v;
	}
	else fusion_info_map[P].push_back(s);
    }
    while(getline(in_,s))
    {
	istr.str(s);
	string g1,g2,temp;
	istr>>temp>>g1>>temp>>temp>>temp
	    >>temp>>temp>>temp>>temp>>g2;
	istr.clear();
	pair<string,string> P = make_pair(g1,g2);

	if(fusion_info_map.find(P) ==fusion_info_map.end())
	{
		vector<string> v(1,s);
		fusion_info_map[P] = v;
	}
	else fusion_info_map[P].push_back(s);
    }
    in.close();
    in_.close();
    cout<<"Fusion Size: "<<fusion_info_map.size()<<endl;

    for(map<pair<string,string>, vector<string> >::iterator i = fusion_info_map.begin(); i!=fusion_info_map.end(); i++)
    {
	for(size_t j=0;j<i->second.size();j++)
	    out<<i->second[j]<<endl;
    }

    return;
}
int main(int argc,char* argv[])
{
    if(argc == 1)
    {
	cout<<"Usage: "<<endl;
	cout<<" mergeAndsort records1 records2 merged.records";
	cout<<endl;
	return 0;
    }
    out.open(argv[3]);
    load_records_and_merge(argv[1],argv[2]);   
    out.close();
    return 0;
}
