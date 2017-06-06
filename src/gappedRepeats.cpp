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

const size_t constGamma = 5;

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

//vergleicht mit LCP-Array und RMQ die Laenge der Arme
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
	const size_t length;
	const vektor_type sa;          //Suffix-Array
	const vektor_type isa;         //inverted LCP-Array
	const vektor_type lcp;         //LCP-Array
    rmq_succinct_sada<> rmq;        //RMQ-Datenstruktur auf Suffix-Array
    const std::string mtext;        //Mirror-Image
    const vektor_type msa;
	const vektor_type misa;
	const vektor_type mlcp;
    rmq_succinct_sada<> mrmq;
    
    lceDataStructure(const std::string& ttext) 
		: text(ttext)
		, length(text.size())
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


//alpha-gapped repeat als Tripel dargestellt
struct alphaGappedRepeat {

	const size_t lArm;		//Startposition linker Arm
	const size_t rArm;		//Startposition rechter Arm
	const size_t length;		//Armlaenge
    
    alphaGappedRepeat(size_t l, size_t r, size_t len) 
		: lArm(l)
		, rArm(r)
		, length(len)
	{
	}    
};


//Berechnet die Startpostionen der Cluster des Suffix-Arrays fuer die Berechnung
//von alpha-gapped repeats mit einer Armlaenge von 1
template<typename vektor_type>
vektor_type calcClusterStarts(const std::string text, const vektor_type sa){
	int size = sa.size();
	vektor_type clusterStarts = new vektor_type(size);
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

//Funktion zur Berechnung alpha-gapped repeats mit Armen der Laenge 1
template<typename vektor_type>
int calc1Arm (lceDataStructure lce, size_t alpha){
	vektor_type clusterStarts = calcClusterStarts(lce.text, lce.sa);
	size_t m = clusterStarts.size();
	size_t n = lce.length;
	vektor_type sa = lce.sa;

	//innerhalb der Cluster des SA nach Startpositionen sortieren
	for (size_t i = 0; i<m-1; i++){
		sort(sa.begin()+clusterStarts[i], sa.begin()+clusterStarts[i+1]);
	}

	//Suche der alpha-gapped repeats mit Armlaenge 1
	for (size_t i = 0; i<n-1; i++){
		for (size_t j = i+1; j<=i+alpha && j<n; j++){
				if (sa[i]+alpha >= sa[j] 
					&& lce.text[sa[i]]==lce.text[sa[j]]){
				//gapRep gefunden
			}
			else{
				j = i + alpha + 1;
			}
		}
	}
	return 1;
}

//Funktion zur Berechnung alpha-gapped repeats mit kurzen Armen
template<typename vektor_type>
int calcShortArm (lceDataStructure lce, size_t alpha){
	size_t n = lce.length;
	size_t sbBegin;				//Anfangsposition von Superblock
	size_t sbEnd;				//Endposition
	vektor_type leftArms;
	size_t raBegin;				//Anfangsposition des rechten Arms
	size_t raEnd;				//Endposition
	//TODO erstelle leere Liste von gapReps oder uebergebe oben einen Pointer auf eine Liste
	
	//TODO richtiger Logarithmus?
	//fuer jeden Superblock
	for (size_t m = constGamma * alpha; m <= n/log2(n) - constGamma -1; m++){
		//TODO Verschiebung richtig?
		//Superblock bestimmen
		if( (m - constGamma * alpha) * log2(n) < 1 ){
			sbBegin = 0;
			sbEnd = (m + constGamma + 1) * log2(n) - 1;
		}
		else{
			sbBegin = (m - constGamma * alpha) * log2(n);
			sbEnd = (m + constGamma + 1) * log2(n) - 1;
		}
		
		//TODO BasisfaktorBitvektor bestimmen (Lemma 22)
		for ( size_t k = 0; k <= log2( constGamma * log2(n) ); k++){
			//Armlaenge auf 2^k fixiert
			
			//gehe alle moeglichen rechten Arme durch
			//TODO rechte Arme so richtig festgelegt?
			for ( size_t j = 1; (j+1) * 2^k - 1 <= sbEnd ; j++ ){
				raBegin = j * 2^k;
				raEnd = (j+1) * 2^k - 1;
				
				//TODO leftArms = Lister der moeglichen linken Arme
				for ( size_t i = 0; i < leftArms.size(); i++ ){
					//TODO linken und rechten Arm erweitern
					//TODO if gueltig: in Liste der gapped repeats einfuegen


					size_t p = lcSuffix(lce, j*2^k -1, leftArms[i] -1); //TODO richtig um 1 verringert?

					size_t s = lcPrefix(lce, (j+1)+2^k , leftArms[i] + 2^k );
	
					//TODO bestimme rl und rr mit Lemma 17
	
					//case a
	
					//case b
	
					//case c
	
					//case d
	
					//case e
	
					//case f
	
					//case g
				}
			}	
		}
	}

	return 1;
}

//TODO
//bekommt rechten Arm und suche moegliche linke Arme
//erstellt dafuer Basisfaktor-Bitvektor
//TODO BasisfaktorBitvektor muss uebergeben werden
template<typename vektor_type>
vektor_type findLeftArms (lceDataStructure lce, size_t sbBegin, size_t sbEnd, size_t raBegin, size_t raEnd){
	//TODO Arme finden
	return 0;
}

//Lemma 17:
//Bekommt LCE-Datenstruktur, Startposition i eines Faktors und Periode p
//gibt Laenge des laengsten Faktors aus
size_t findLongestPeriod (lceDataStructure lce, size_t i, size_t p){
	size_t length = lce.lcp[lce.rmq(i+1,i+p)];
	return length;
}

//TODO um "schnellere" LCE-Anfragen erweitern
// gibt longest common prefix von 2 Woertern aus, die an Position i und j beginnen
int lcPrefix(lceDataStructure lce, size_t i, size_t j){    
    return lce.lcp[ lce.rmq( lce.isa[i]+1, lce.isa[j] ) ];
} 

//TODO um "schnellere" LCE-Anfragen erweitern
// gibt longest common suffix von 2 Woertern aus, die an Position i und j beginnen
int lcSuffix(lceDataStructure lce, size_t i, size_t j){  
	size_t n = lce.length;  
    return lce.mlcp[ lce.mrmq( lce.misa[n-i-1]+1, lce.misa[n-j-1] ) ];
} 

//TODO wird nicht verwendet sondern in calcShortArm direkt ausgefuehrt
//Funktion erweitert gapped repeat maximal nach links und rechts
//bekommt LCE-Datenstruktur, Start- und Endposition des Superblocks, 
//		  k (Armlaenge 2^k), j, "Startpositionen" beider y-Arme
/*
alphaGappedRepeat* extendArms (lceDataStructure lce, size_t sbBegin, size_t sbEnd, 
								size_t k, size_t j, size_t raBegin, size_t laBegin) {
	alphaGappedRepeat* rep = 0; //new alphaGappedRepeat(0,0,0);

	size_t p = lcSuffix(lce, j*2^k -1, laBegin -1); //TODO richtig um 1 verringert?

	size_t s = lcPrefix(lce, (j+1)+2^k , laBegin + 2^k );
	
	//TODO bestimme rl und rr mit Lemma 17
	
	//case a
	
	//case b
	
	//case c
	
	//case d
	
	//case e
	
	//case f
	
	//case g
		
	return rep;
}
*/

//Funktion zur Berechnung alpha-gapped repeats mit kurzen Armen
template<typename vektor_type>
int calcLongArm (lceDataStructure lce, size_t alpha){
	size_t n = lce.length;
	vektor_type wB = 0; //TODO bestimme Block-Representation
	//TODO lceB = lce-Datenstruktur der Blockrep
	for ( size_t k = 0; k <= (n/log2(n)); k++){
		vektor_type kBlocks = 0; //TODO kBloecke Bestimmen
		for ( size_t i = 0; i < kBlocks.length; i++){
			for ( size_t start = 0; start < log2(n); start++ ){
				//TODO linken Arm y fixieren
				//TODO rechte Arme mit Binaersuche finden
				//TODO alle rechten Arme pruefen und erweitern
			}
		}
	}
	
	return 0;
}
