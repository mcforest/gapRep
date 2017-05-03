#include <iostream>
#include <string>
#include <math.h>
#include "strinalyze.cpp"
#include <iomanip>
#include <glog/logging.h>
#include "substring.hpp"
#include "checked_vector.hpp"
#include <sdsl/rmq_support.hpp>
#include <boost/utility/string_ref.hpp>



using namespace std;
using namespace sdsl;



//vergleicht die Laenge der Arme ueber Zeichenvergleiche
//nutzt word packing mit uint64_t, "verbesserte" Variante
//i und j sind Positionen im Text, length ist bereits die geforderte Mindestlaenge
int naivComp(const char* text, size_t i, size_t j, size_t textLength){
    bool gr = 1;
    uint64_t* l;
    uint64_t* r;
//     char* lbuff[8];
//     char* rbuff[8];
    size_t k;
    size_t length = 0;
    
    for (k=0; gr==1; k++){
        if (i + 8*k <= j && j + 8*k <= textLength){
//             memcpy( lbuff, &text[i+8*k], 8);
//             memcpy( lbuff, &text[j+8*k], 8);
//             l = (uint64_t*) (lbuff);
//             r = (uint64_t*) (rbuff);
            l = (uint64_t*) (text+(i+8*k));
            r = (uint64_t*) (text+(j+8*k));
            if (l[0]!=r[0]){
                gr = 0;
                length = k*8;
                while ( text[i+length] == text[j+length] && i + length <= j && j + length <= textLength ){
                    length++;
                }
            }
        }
        else {
            gr = 0;
            length = (k-1)*8;
            while ( text[i+length+1] == text[j+length+1] && i + length+1 <= j && j + length+1 <= textLength){
                length++;
            }
        }
        
    }
    return length;
}

//vergleicht mit LCP-Array die Laenge der Arme
//i und j sind Positionen im S/LCP-Array
int lcpMin(const StringStats &stats, size_t i, size_t j){ //bei Uebergabe muss i<j gelten
    int length = abs(stats.sa[j] - stats.sa[i]);
    for(size_t k = i+1; k<=j; k++){
        length = min(length, stats.lcp[k]);
    }
    return length;
}

//vergleich mit LCP-Array und RMQ die Laenge der Arme
//i und j sind Positionen im S/LCP-Array
int lcpRmqMin(const StringStats &stats, rmq_succinct_sada<> &rmq, size_t i, size_t j){
    
    return stats.lcp[rmq(i+1,j)];
} 

//LCE-Datenstruktur
struct lceDataStructure {
#ifdef NDEBUG
		typedef checked_vector<int> vektor_type;
#else
		typedef std::vector<int> vektor_type;
		//#define vektor_type sdsl::int_vector<>
#endif

	const std::string text;        //Text
	const vektor_type sa;          //Suffix-Array
	const vektor_type isa;         //inverted Suffix-Array
	const vektor_type lcp;         //LCP-Array
    rmq_succinct_sada<> rmq;        //RMQ-Datenstruktur auf Suffix-Array
    const std::string mtext;        //Mirror-Image
    const vektor_type msa;
	const vektor_type misa;
	const vektor_type mlcp;
    rmq_succinct_sada<> mrmq;
    
    lceDataStructure(const std::string&& ttext) 
		: text(ttext)
		, sa(create_sa<vektor_type>(text, FLAGS_stripDollar))
		, isa(inverse<vektor_type>(sa))
		, lcp(create_lcp<vektor_type>(text, sa, isa))
        , rmq(&lcp)
        , mtext(string ( text.rbegin(), text.rend() ))
        , msa(create_sa<vektor_type>(mtext, FLAGS_stripDollar))
		, misa(inverse<vektor_type>(msa))
		, mlcp(create_lcp<vektor_type>(mtext, msa, misa))
        , mrmq(&mlcp)
	{
	}    
};

//Berechnet die Startpostionen der Cluster des Suffix-Arrays fuer die Berechnung
//von alpha-gapped repeats mit einer Armlaenge von 1
template<typename vektor_type>
int* calcClusterStarts(const std::string text, const vektor_type sa){
	int size = sa.size();
	int *clusterStarts = new int[size];
	clusterStarts[0]=0;
	int pos = 1;
	for(int i=1; i < size; i++){
		if(text[sa[i-1]] != text[sa[i]]){
			clusterStarts[pos]=i;
			pos++;
		}
	}


//	delete [] clusterStarts;
	return clusterStarts;
}


//Cluster des Suffix-Arrays fuer die Berechnung von alpha-gapped repeats
//mit einer Armlaenge von 1
struct saClusters {
#ifdef NDEBUG
		typedef checked_vector<int> vektor_type;
#else
		typedef std::vector<int> vektor_type;
		//#define vektor_type sdsl::int_vector<>
#endif

	const vektor_type sa;       //Suffix-Array des Texts
	int* clusterStarts;      	//Array mit den Startpositionen der einzelnen Cluster
    
    saClusters(const lceDataStructure lce) 
		: sa(lce.sa)
		, clusterStarts(calcClusterStarts(lce.text, sa))
	{
	}    
};
