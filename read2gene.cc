#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<algorithm>
#include<vector>
#include <utility>
#include<map>
using namespace std;



struct read_info
{
    string read_id;
    string strand;
    int soft_left;
    int soft_right;
    int mapped_length;
};
struct gene_info
{
    int left;
    int right;
};


map<string, map<pair<int,int>, string> > Chr_Range_gid_map;
//map<pair<int,int>, vector<string> >  Read_Range_map;
map<string, map<pair<int,int>, vector<read_info> > > Chr_Read_map;
map<string,vector<pair<string,string> > > CandRead_rinfo_gid_map;

map<string, gene_info> gene_info_map;
map<string, read_info> read_info_map;



map<string, map<vector<int>, bool> > Chr_vecExon_map;
ofstream outcand;
map<string,bool> read_with_prim_and_supp_map;

void load_gene(char* file)
{
    cout<<"Loading gene..."<<endl;
    ifstream in(file);
    string s;
    istringstream istr;
    while(getline(in,s))
    { 
 	replace(s.begin(), s.end(), ',', ' ');
	istr.str(s);
	string gene_id, chr, strand, gl_, gr_;
  	int gl,gr;
	istr>>gene_id>>chr>>gl>>gr>>strand;
	//cout<<gene_id<<" "<<chr<<" "<<gl<<" "<<gr<<" "<<strand<<endl;
	istr.clear();
	//chr=chr+strand;	
	gene_id += strand;
	pair<int,int> range = make_pair(gl,gr);
	if(Chr_Range_gid_map.find(chr) == Chr_Range_gid_map.end())
	{
	    map<pair<int,int>, string> range_gid_map;
	    range_gid_map[range] = gene_id;
	    Chr_Range_gid_map[chr]=range_gid_map;
	}
	else
	{
	    Chr_Range_gid_map[chr][range] = gene_id;
	}
    }
    cout<<Chr_Range_gid_map.size()<<endl;
}

void load_read(char*file)
{
    cout<<"Loading reads..."<<endl;
    ifstream in(file);
    string s;
    istringstream istr;
    //if primary and supplementary
    map<string, pair<bool,bool> > read_prim_and_supp_map;
    int i = 0;
    while(getline(in,s))
    {
	if(i % 500000 ==0) cout<<" loading "<<i<<" lines"<<endl;
	i++;
  	istr.str(s);
	string read_id, chr,strand,temp;
        int start_pos,end_pos, soft_l, soft_r, mapped_len;
        string flag, supp_flag;
        istr>>read_id>>chr>>strand>>start_pos>>end_pos>>soft_l>>soft_r>>temp>>flag>>supp_flag>>mapped_len;
        istr.clear();
	if(flag == "secondary") continue;
	
	if(read_prim_and_supp_map.find(read_id) == read_prim_and_supp_map.end())
	{
		pair<bool,bool> flag = make_pair(0,0);
		if(supp_flag == "not_supplementary") 
			flag.first = 1;
		else flag.second = 1;
		read_prim_and_supp_map[read_id] = flag;
	}
	else 
	{
		if(supp_flag == "not_supplementary")
			 read_prim_and_supp_map[read_id].first = 1;
		else read_prim_and_supp_map[read_id].second = 1;
	}
       
    }
    int Size = 0;
    for(map<string, pair<bool,bool> >::iterator i = read_prim_and_supp_map.begin(); i!= read_prim_and_supp_map.end(); i++)
    {
	if(i->second.first & i->second.second)
		Size++;
    }
    cout<<"Potential candidates size: "<<Size<<endl;
    in.close();
    in.open(file);
    cout<<"Getting reads info..."<<endl;
    i = 0;
    while(getline(in,s))
    {
	if(i % 500000 ==0) cout<<" loading "<<i<<" lines"<<endl;
        i++;
	istr.str(s);
	string read_id, chr,strand,temp;
	int start_pos,end_pos, soft_l, soft_r, mapped_len;
	string flag, supp_flag;
	istr>>read_id>>chr>>strand>>start_pos>>end_pos>>soft_l>>soft_r>>temp>>flag>>supp_flag>>mapped_len;
	istr.clear();
	
	if(flag == "secondary") continue;
	
	map<string, pair<bool,bool> >::iterator it = read_prim_and_supp_map.find(read_id);

	if(it->second.first & it->second.second)//primary and supplementary
	{
	    if(supp_flag == "not_supplementary") read_id += "N";
	    else read_id += "Y";

	    
	    //chr+=strand;

	    pair<int,int> read_range = make_pair(start_pos, end_pos);
	    read_info ri={read_id,strand,soft_l,soft_r,mapped_len};
	    if(Chr_Read_map.find(chr) == Chr_Read_map.end())
            {
              vector<read_info> reads(1,ri);
              map<pair<int,int>, vector<read_info> > range_readsid_map;
	      range_readsid_map[read_range] = reads;
              Chr_Read_map[chr] = range_readsid_map;

	      
            }
            else
            {
              if(Chr_Read_map[chr].find(read_range) == Chr_Read_map[chr].end())
              {
                vector<read_info> reads(1,ri);
                Chr_Read_map[chr][read_range] = reads;
              }
              else Chr_Read_map[chr][read_range].push_back(ri);
            }
/*
	    if(Read_Range_map.find(read_range) == Read_Range_map.end())
	    {
		vector<string> reads(1,read_id);
		Read_Range_map[read_range] = reads;
	    }
	    else Read_Range_map[read_range].push_back(read_id);
*/
	}	


    }
    //cout<<Read_Range_map.size()<<endl;
}


void process_one_gene()
{
}

typedef map<string, map<pair<int,int>, string> >::iterator iter1;//chr gene
typedef map<pair<int,int>, vector<read_info> >::iterator iter2;
typedef map<string, map<pair<int,int>, vector<read_info> > >::iterator iter3;//chr read
void get_candfusion()
{
 
    cout<<"Getting candidate fusion..."<<endl;
    for(iter3 m = Chr_Read_map.begin(); m != Chr_Read_map.end(); m++) //for each chromosome: m->first
    {

      map<pair<int,int>, vector<read_info> > & Read_Range_map = m->second;

      iter1 i = Chr_Range_gid_map.find(m->first);
      if(i == Chr_Range_gid_map.end()) continue;
      map<pair<int,int>, string> & Range_map = i->second;//the gene range of chromosome: m->first

      map<pair<int,int>, string>::iterator K = Range_map.begin();

      for(iter2 j = Read_Range_map.begin(); j != Read_Range_map.end();j++)//for each read;
      {

	    int rl = j->first.first, rr = j->first.second;

	    map<pair<int,int>, string>::iterator k = K;	    
	    //map<pair<int,int>, string>::iterator k = Range_map.begin();	    
	    for(;k!=Range_map.end();k++)
	    {
		vector<read_info> reads;
		//vector<read_info> reads_=j->second;
		/*
		cout<<" read-range: "<<rl<<" "<<rr<<" - ";
		for(size_t a = 0;a<reads_.size();a++) cout<<reads_[a]<<" ";
		cout<<endl;

		cout<<" gene-range: "<<k->first.first<<" "<<k->first.second<<" - "<<k->second<<" - "<<m->first<<endl;
		*/
		//if(k->first.first <= rl ) K = k;
		if(k->first.first <= rl && rr <= k->first.second)
		{
			reads = j->second;
			//next read may appear in the previous gene.
			int L = k->first.first;
			int M = 0;
			for(map<pair<int,int>, string>::iterator b = k; b != Range_map.begin();b--)
			{
				//cout<<"  "<<b->second<<": "<<endl;
				//cout<<"  "<<L<<" "<<b->first.second<<endl;
				M++;
				if(b != Range_map.begin() && (b->first.second >= L || M <= 5) )
				{
				    //cout<<"  "<<L<<" "<<b->first.second<<endl;
				    if(L > b->first.first) 
					L = b->first.first;
				}
				else {
				    K = b;
				    break;
				}
			} 
		}
		else if( k->first.first > rl) 
	  	{
			break;
		}
		
		if(!reads.empty())
		{
			/*
			cout<<i->first<<" "<<k->second<<": ";
		  	for(size_t a = 0;a<reads.size();a++)
                	{
                            cout<<reads[a]<<" ";
                	}
                	cout<<endl;
			*/
			for(size_t a = 0;a<reads.size();a++)
			{
			    string gene_id_ = k->second;
			    string gene_id = gene_id_.substr(0,gene_id_.length() - 1);
		 	    string gene_str = gene_id_.substr(gene_id_.length() - 1,1);//gene strand

			    string read_info=m->first;
			    stringstream ss;
			    ss<<rl<<" "<<rr<<" ";
			    //ss<<reads[a].strand<<" ";
			    ss<<gene_str<<" ";//strand
			    ss<<reads[a].soft_left<<" "<<reads[a].soft_right<<" "<<reads[a].mapped_length;
			    read_info = read_info + " " + ss.str();

			    string rid = reads[a].read_id;

			    

			    pair<string,string> rinfo_gid = make_pair(read_info,gene_id);

			
			    if(CandRead_rinfo_gid_map.find(rid) == CandRead_rinfo_gid_map.end())
			    {	
				vector<pair<string,string> > vec_rinfo_gid(1,rinfo_gid);
			      	CandRead_rinfo_gid_map[rid] = vec_rinfo_gid;
			    }
			    else CandRead_rinfo_gid_map[rid].push_back(rinfo_gid);
			}

		}
		
	    }//for each gene range
	    
      }//for each read 
    }//for each chromosome    
    cout<<"CandRead_rinfo_gid_map.size: "<<CandRead_rinfo_gid_map.size()<<endl;
    
    map<string,vector<pair<string,string> > >::iterator i = CandRead_rinfo_gid_map.begin();
    for(;i!=CandRead_rinfo_gid_map.end();i++)
    {
	string read_id = i->first;
	vector<pair<string,string> > rinfo_genes = i->second;
	if(read_id[read_id.length() - 1] == 'N')//primary&NotSupplementary
	{
		string read_id_ = read_id.substr(0,read_id.length() - 1) + "Y";
		map<string,vector<pair<string, string> > >::iterator j = CandRead_rinfo_gid_map.find(read_id_);

		if(j != CandRead_rinfo_gid_map.end())
		{
		    vector<pair<string,string> > rinfo_genes_ = j->second;

		    if(rinfo_genes.size() == 1 && rinfo_genes_.size() == 1)
		    {//2025.1.6

		      for(int k = 0;k<rinfo_genes.size();k++)
		      {
			for(int m = 0;m<rinfo_genes_.size();m++)
			{
			    istringstream istr,istr_;
			    istr.str(rinfo_genes[k].first);
			    istr_.str(rinfo_genes_[m].first);

			    int sl_k,sr_k,sl_m,sr_m;//softclip left; softclip right; k and m 
			    string temp;
			    istr>>temp>>temp>>temp>>temp>>sl_k>>sr_k;
			    istr_>>temp>>temp>>temp>>temp>>sl_m>>sr_m;
			    if(sl_k < sl_m && sr_k > sr_m)
			    {
				//outcand<<sl_k<<" "<<sr_k<<" - "<<sl_m<<" "<<sr_m<<endl;
			    	outcand<<read_id.substr(0,read_id.length() - 1)<<" "
			    	<<rinfo_genes[k].second<<" "<<rinfo_genes[k].first<<" "
			    	<<rinfo_genes_[m].second<<" "<<rinfo_genes_[m].first<<endl;
			    }
			    else if(sl_k > sl_m && sr_k < sr_m)
			    {
				outcand<<read_id.substr(0,read_id.length() - 1)<<" "
			 	<<rinfo_genes_[m].second<<" "<<rinfo_genes_[m].first<<" "
				<<rinfo_genes[k].second<<" "<<rinfo_genes[k].first<<endl;
			    }
			}
		      }
		    }//2025.1.6
		}
	}
	else continue;
    }

}
int main(int argc,char* argv[])
{
    if(argc == 1)
    {
	cout<<"Usage: "<<endl;
	cout<<" read2gene sum_gene_info_sort_nonrepeat.txt1  alignment.info1 candidates ";
	cout<<endl;
	return 0;
    }
    load_gene(argv[1]);
    load_read(argv[2]);
    outcand.open(argv[3]);
    get_candfusion();
    outcand.close();

    return 0;
}
