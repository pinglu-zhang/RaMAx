#include "wfmash_router.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
using namespace WfmashRouterDetail;
using Pair = std::tuple<std::string,uint64_t,std::string,uint64_t>;
using Pairs = std::set<Pair>;
void check(bool v, const char* message) { if (!v) throw std::runtime_error(message); }
ParsedPafRecord record(std::string t, uint64_t ts, uint64_t qs, std::string cg, bool reverse=false) {
    uint64_t q=0,n=0,c=0,m=0,x=0; std::istringstream in(cg); char op;
    while (in>>x>>op) { c+=x; if(op!='I') n+=x; if(op!='D')q+=x; if(op=='='||op=='M')m+=x; }
    std::ostringstream s; s<<"query\t10000\t"<<qs<<'\t'<<qs+q<<'\t'<<(reverse?'-':'+')<<'\t'<<t<<"\t10000\t"<<ts<<'\t'<<ts+n<<'\t'<<m<<'\t'<<c<<"\t60\tcg:Z:"<<cg;
    return parsePafLine(s.str(),true);
}
Pairs pairs(const std::vector<ParsedPafRecord>& records) {
    Pairs p; std::map<std::string,std::set<uint64_t>> qt,tt;
    for (const auto& r:records) {
        uint64_t t=r.target_start,q=r.strand==Strand::FORWARD?r.query_start:r.query_end;
        for(uint64_t x=r.query_start;x<r.query_end;++x) check(qt[r.query_name].insert(x).second,"query collision");
        for(uint64_t x=r.target_start;x<r.target_end;++x) check(tt[r.target_name].insert(x).second,"target collision");
        size_t paired=0;
        for(auto u:r.cigar) {char op;uint32_t n;intToCigar(u,op,n);for(uint32_t j=0;j<n;++j){
            if(op!='D'&&r.strand!=Strand::FORWARD)--q;
            if(op!='D'&&op!='I'){p.emplace(r.target_name,t,r.query_name,q);++paired;}
            if(op!='I')++t;
            if(op!='D'&&r.strand==Strand::FORWARD)++q;
        }}
        check(t==r.target_end,"target consumption");check(q==(r.strand==Strand::FORWARD?r.query_end:r.query_start),"query consumption");check(paired>0,"pure gap fragment");
    } return p;
}
void emit(const ParsedPafRecord& r) {
    std::cout<<r.query_name<<'\t'<<r.query_length<<'\t'<<r.query_start<<'\t'<<r.query_end<<'\t'<<(r.strand==Strand::FORWARD?'+':'-')<<'\t'<<r.target_name<<'\t'<<r.target_length<<'\t'<<r.target_start<<'\t'<<r.target_end<<'\t'<<r.matches<<'\t'<<r.block_length<<'\t'<<r.mapq<<"\tcg:Z:"<<r.cigar_text<<'\n';
}
int main(int argc,char**argv) {try {
    if(argc==2){std::ifstream f(argv[1]);check(bool(f),"input missing");std::string line;std::vector<ParsedPafRecord> r;
        while(std::getline(f,line))if(!line.empty())r.push_back(parsePafLine(line,true));
        normalizePafForGraph(r);for(const auto& x:r)emit(x);return 0;}
    std::mt19937 rng(20261003);
    for(int reverse=0;reverse<2;++reverse)for(int trial=0;trial<500;++trial){
        std::string cg="3I";for(int j=0;j<25;++j){char ops[]={'=','X','M','I','D'};cg+=std::to_string(1+rng()%7)+ops[rng()%5];}cg+="4=3D";
        auto raw=record("B",20,10,cg,reverse);
        std::vector<ParsedPafRecord> input;std::set<uint64_t> occupied;uint64_t target=0;
        for(uint64_t q=10;q<raw.query_end;){if(rng()%3==0){auto n=std::min<uint64_t>(1+rng()%5,raw.query_end-q);input.push_back(record("A",target,q,std::to_string(n)+"="));target+=n+1;for(uint64_t k=q;k<q+n;++k)occupied.insert(k);q+=n;}else ++q;}
        Pairs expected=pairs(input);auto rawpairs=pairs({raw});
        for(const auto& p:rawpairs)if(!occupied.count(std::get<3>(p)))expected.insert(p);
        input.push_back(raw);auto copy=input;normalizePafForGraph(input);normalizePafForGraph(copy);
        check(pairs(input)==expected,"basewise oracle mismatch");check(pairs(copy)==pairs(input),"not deterministic");
        check(input.size()==copy.size(),"record count nondeterministic");
        for(size_t i=0;i<input.size();++i)check(input[i].cigar_text==copy[i].cigar_text&&input[i].query_start==copy[i].query_start&&input[i].target_start==copy[i].target_start,"record nondeterministic");
    }
    // Target-only middle overlap; preserve both the accepted blocker and suffix.
    auto a=record("T",40,100,"10=");auto b=record("T",0,0,"100=");
    std::vector<ParsedPafRecord> v{a,b};auto st=normalizePafForGraph(v);check(st.recovered_fragments==2,"target split");check(pairs(v).size()==100,"target split pairs");
    // A later baseline record is reserved before a rejected record is rescued.
    auto c=record("U",0,20,"10=");v={a,b,c};normalizePafForGraph(v);auto got=pairs(v);for(auto p:pairs({c}))check(got.count(p),"baseline later record lost");
    // Duplicates, containment, complete exhaustion and chromosome endpoints.
    a=record("T",0,0,"100=");b=record("T",10,10,"20=");v={a,a,b};normalizePafForGraph(v);check(v.size()==1&&pairs(v).size()==100,"duplicate/containment");
    a=record("T",9990,9990,"10=",true);v={a,a};normalizePafForGraph(v);check(v.size()==1&&pairs(v).size()==10,"boundary reverse");
    a=record("A",0,0,"5="); b=record("B",0,0,"10M");b.matches=6;
    v={a,b};normalizePafForGraph(v);check(v.size()==2&&v[1].matches==3,"M aggregate match estimate");
    a=record("A",0,0,"3=");b=record("B",0,0,"3=2I2D");v={a,b};normalizePafForGraph(v);check(v.size()==1,"pure gap suffix");
    std::cout<<"PASS: 1000 randomized strand/CIGAR oracle cases and target overlap, priority, duplicate, containment, boundary cases\n";
} catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;} }
