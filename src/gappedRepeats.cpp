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
#include "../build/src/grenzwerte.hpp"




using namespace std;
using namespace sdsl;


const size_t constGamma = 2;

//vergleicht die Laenge der Arme ueber Zeichenvergleiche
//nutzt word packing mit uint64_t, "verbesserte" Variante
//i und j sind Positionen im Text, length ist bereits die geforderte Mindestlaenge
//gibt tatsaechlichen laengsten Faktor ohne Ueberschneidung zurueck
int naivComp(const char* text, size_t i, size_t j, size_t textLength){
    bool gr = 1;		//Abbruchbedingung der Schleife (wenn linkes Teilwort das rechte schneidet oder ungleich ist)
    uint64_t* l;
    uint64_t* r;
    size_t k;
    size_t length = 0;
    
	
    for (k=0; gr==1; k++){
        if (i + 8*(k+1) < j && j + 8*(k+1) < textLength){

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
            length = k*8;
            while ( text[i+length] == text[j+length] && i + length < j && j + length < textLength){
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
    int length = abs(stats.sa[j] - stats.sa[i]);
    return min(length, stats.lcp[rmq(i+1,j)]);
}

string invertieren (string text){
	string invert;
	for (int i = text.size()-1 ; i >= 0 ; i--){
		invert = invert + text[i];
	}
	return invert;
}

//LCE-Datenstruktur
struct lceDataStructure {
#ifdef NDEBUG
		typedef checked_vector<int> vektor_type;
#else
		typedef std::vector<int> vektor_type;
		//#define vektor_type sdsl::int_vector<>
#endif

	const std::string &text;        //Text
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
		, length(ttext.size())
		, sa(create_sa<vektor_type>(text, FLAGS_stripDollar))
		, isa(inverse<vektor_type>(sa))
		, lcp(create_lcp<vektor_type>(text, sa, isa))
        , rmq(&lcp)
		, mtext( invertieren(ttext) )
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

bool isValid(alphaGappedRepeat *gr, size_t alpha){
	return ( gr->lArm + alpha * gr->length >= gr->rArm && gr->lArm + gr->length <= gr->rArm && gr->length > 0 );
}

//Gibt eine alphaGappedRepeat in der Konsole aus
int printGappedRepeat(alphaGappedRepeat *gr){
	cout << "(" << gr->lArm << "," << gr->rArm << "," << gr->length << ")" << endl;
	return 0;
}

//Gibt ein Array aller alphaGappedRepeats in der Konsole aus
template<typename alphaGappedRepeat>
int printGappedRepeat(vector<alphaGappedRepeat*>& vec){
	for(size_t i = 0; i<vec.size(); i++){
		printGappedRepeat(vec[i]);
	}
	return 0;
}

//vergleicht die Laenge der Arme ueber Zeichenvergleiche
//nutzt word packing mit uint64_t, "verbesserte" Variante
//i und j sind Positionen im Text, length ist bereits die geforderte Mindestlaenge
//gibt lcp-Wert aus, ignoriert ueberschneidung
int naivLCP(const char* text, size_t i, size_t j, size_t textLength){
    bool gr = 1;		//Abbruchbedingung der Schleife (wenn linkes Teilwort das rechte schneidet oder ungleich ist)
    uint64_t* l;
    uint64_t* r;
    size_t k;
    size_t length = 0;
    
	
    for (k=0; gr==1; k++){
        if (j + 8*(k+1) < textLength){

            l = (uint64_t*) (text+(i+8*k));
            r = (uint64_t*) (text+(j+8*k));
            if (l[0]!=r[0]){
                gr = 0;
                length = k*8;
                while ( text[i+length] == text[j+length] && j + length <= textLength ){
					
                    length++;
                }
            }
        }
        else {
            gr = 0;
            length = k*8;
            while ( text[i+length] == text[j+length] && j + length < textLength){
                length++;
            }
        }
        
    }
    return length;
}


// gibt longest common prefix von 2 Woertern aus, die an Position i und j beginnen
int lcPrefix(lceDataStructure*& lce, size_t i, size_t j){
	size_t left = min(lce->isa[i],lce->isa[j]);
	size_t right = max(lce->isa[i],lce->isa[j]);
	if (left == right){
		return 0;
	}
	size_t realX = j-i;
	size_t realY = abs(right - left);
	
	if( realX < x && (realX * y) < (realY * x) ){         //wenig zeichenvergleiche noetig	
		//durchlaeuft Text
		return naivLCP(lce->text.c_str(),i,j,lce->length); 
	}
	else if( realY < y){                         //kurzer Abstand im Suffix-Array
		//durchlaeuft LCP-Array
		int length = lce->lcp[left+1];
    	for(size_t k = left+1; k<=right; k++){
        	length = min(length, lce->lcp[k]);
    	}
		return length;
	}
	else{
		//benutzt RMQ
    	return lce->lcp[ lce->rmq( left+1, right ) ];
	}

} 



// gibt longest common suffix von 2 Woertern aus, die an Position i und j enden
int lcSuffix(lceDataStructure*& lce, size_t i, size_t j){ 
	size_t n = lce->length;
	size_t left = min( lce->misa[n-i-1], lce->misa[n-j-1] );
	size_t right = max( lce->misa[n-i-1], lce->misa[n-j-1] );
	if (left == right){
		return 0;
	}
	size_t realX = j-i;
	size_t realY = abs(right - left);
	
	if( realX < x && (realX * y) < (realY * x) ){         //wenig zeichenvergleiche noetig	
		//durchlaeuft Text
		return naivLCP(lce->mtext.c_str(),n-j-1,n-i-1,lce->length); 
	}
	else if( realY < y){                         //kurzer Abstand im Suffix-Array
		//durchlaeuft LCP-Array
		int length = lce->mlcp[left+1];
    	for(size_t k = left+1; k<=right; k++){
        	length = min(length, lce->mlcp[k]);
    	}
		return length;
	}
	else{
		//benutzt RMQ
    	return lce->mlcp[ lce->mrmq( left+1, right ) ];
	}

} 


//Berechnet die Startpostionen der Cluster des Suffix-Arrays fuer die Berechnung
//von alpha-gapped repeats mit einer Armlaenge von 1
//template<typename vektor_type>
vector<int> calcClusterStarts(const std::string text, const vector<int> sa){
	int size = sa.size();
	vector<int> clusterStarts;
	clusterStarts.reserve(size);
	clusterStarts.push_back(0);
	size_t pos = 1;
	for(int i=1; i < size; i++){
		if(text[sa[i-1]] != text[sa[i]]){
			clusterStarts.push_back(i);
			pos++;
		}
	}

	return clusterStarts;
}

//Funktion zur Berechnung alpha-gapped repeats mit Armen der Laenge 1
//template<typename vektor_type>
int calc1Arm(lceDataStructure*& lce, size_t alpha, vector<alphaGappedRepeat*> *grList){
	size_t n = lce->text.size();
	vector<int> clusterStarts = calcClusterStarts(lce->text, lce->sa);
	clusterStarts.push_back(n+1);
	size_t m = clusterStarts.size();
	vector<int> sa = lce->sa;
	alphaGappedRepeat* gappedRep;

	
	//innerhalb der Cluster des SA nach Startpositionen sortieren
	for (size_t i = 0; i<m-1; i++){
		sort(sa.begin()+clusterStarts[i], sa.begin()+clusterStarts[i+1]);

	}


	//Suche der alpha-gapped repeats mit Armlaenge 1
	for (size_t i = 0; i<n; i++){
		//alpha*2, da unser suffix array groeßer ist als angenommen 
		//(da der text einmal vorw und einmal rueckw gespeichert ist
		for (size_t j = i+1; j<=i+alpha && j<n+1; j++){
			
				if (lce->text[sa[i]]==lce->text[sa[j]]){
					if ( ( i==0 && lce->text[sa[i]+1]!=lce->text[sa[j]+1] )
						|| ( j==n-1 && lce->text[sa[i]-1]!=lce->text[sa[j]-1] )
						|| i+1 == j
						|| ( lce->text[sa[i]-1]!=lce->text[sa[j]-1] && lce->text[sa[i]+1]!=lce->text[sa[j]+1] ) ) {
						gappedRep = new alphaGappedRepeat(sa[i], sa[j] , 1);
						if (gappedRep->lArm <= lce->length && gappedRep->rArm <= lce->length){
							if (isValid(gappedRep, alpha)){
								//grList->push_back(gappedRep);	//auskommentiert fuer Messungen
							}
						}
						delete gappedRep;
						gappedRep =0;
					}
				}
				
				//sobald cluster verlassen wird nicht weiter pruefen
				else{
					j = i + alpha + 1;
				}
				
		}
	}
	return 0;
}

//Lemma 17 aus Paper von Gawrychowski:
//Bekommt LCE-Datenstruktur, Startposition i eines Faktors und Periode p
//gibt Laenge des laengsten Faktors aus
size_t findLongestPeriod (lceDataStructure*& lce, size_t i, size_t p){
	size_t length = lcPrefix(lce, i, i+p);
	if ( p <= length ){
		return length+p;
	}
	else{
		return 0;
	}
}


//KMP Algorithmus findet alle Vorkommen von y im Superblock
//Startposition y_rho ist letzter Wert im ausgegebenen Array
vector<int> kmpMatching (lceDataStructure*& lce, size_t sbStart, size_t sbEnd, size_t raStart, size_t raLen){
	
	const char* text = lce->text.c_str();
	string rArm = (string)(text+raStart);
    vector<int> lcs(raLen + 1, -1);
	vector<int> lArms;

	//Pattern vorverarbeiten
	for(int i = 1; i <= (int)raLen; i++)
	{
		int pos = lcs[i - 1];
		while(pos != -1 && rArm[pos] != rArm[i - 1]) pos = lcs[pos];
		lcs[i] = pos + 1;
	}

	int textPos = sbStart;
	int armPos = 0;

	//Pattern suchen
	while(textPos < (int)sbEnd+1 && textPos < (int)(raStart+raLen) )
	{
		while(armPos != -1 && (armPos == (int)raLen || rArm[armPos] != text[textPos])) armPos = lcs[armPos];
		armPos++;
		textPos++;
		if(armPos == (int)raLen) lArms.push_back(textPos - raLen);
	}
	return lArms;
}




//Funktion zur Berechnung alpha-gapped repeats mit kurzen Armen
//schnellere Berechnung fuer Perioden bisher nicht enthalten
//template<typename vektor_type>
int calcShortArm (lceDataStructure*& lce, size_t alpha, vector<alphaGappedRepeat*> *grList, bool v2){
	alphaGappedRepeat* gappedRep;
	size_t n = lce->length;
	size_t sbBegin;				//Anfangsposition von Superblock
	size_t sbEnd;				//Endposition

	vector<int> leftArms;		//Liste der moeglichen linken Arme
	size_t raBegin;				//Anfangsposition des rechten Arms
	size_t laBegin;
	
	int a;						//Suffix des Arms (Laenge)
	int s;						//Prefix des Arms (Laenge)

	size_t p;					//Periode
	size_t rLLength;			//r_lambda laenge (startet an laBegin)
	size_t rRLength;			//r_rho laenge (startet an raBegin)
	size_t helpLength;
	bool oneblock = false;		//falls Text nur einen Superblock enthält
	
	vector<int> *occs;			// occs[h]=ph wenn leftArms[h] erstes Vorkommen des Runs ist und ph die Periode
								// occs[h]=0 wenn leftArms[h] spaeteres Vorkommen eines Runs ist
								// occs[h]=-1 wenn leftArms[h] single occurance ist
								// occs[h]=-2 wenn leftArms[h] im selben run wie y_rho ist
								// 				dann gilt p = leftArms[n-1] - leftArms[n-2] mit n=leftArms.size()
	size_t ph;					//vorlaeufige Periode p
	size_t lcHelp;
	size_t lastPos;				//speichert letzte Position des aktuellen runs
	size_t rRhoBegin;

	

	if ( constGamma*alpha <= n/log2(n) - constGamma -1 || v2 ){
		
		oneblock = true;
	}

	//fuer jeden Superblock
	for (size_t m = constGamma * alpha; m <= n/log2(n) - constGamma -1 || oneblock; m++){


		//Superblock bestimmen
		if (oneblock || v2){ //falls Text nur einen Superblock enthält (zu kurz oder v2)
			sbBegin = 0;
			sbEnd = n-1;
			oneblock = false; 	//damit Schleife ueber m nur einmal durchlaufen wird, wenn text zu kurz ist
			m = n/log2(n) +1; 
		}
		else if( (m - constGamma * alpha) * log2(n) < 1 ){
			sbBegin = 0;
			sbEnd = (m + constGamma + 1) * log2(n) - 1;
		}
		else{
			sbBegin = (m - constGamma * alpha) * log2(n);
			sbEnd = (m + constGamma + 1) * log2(n) - 1;
		}
		
		for ( size_t k = 0; k <= log2( constGamma * log2(n) ); k++){
			//Armlaenge auf 2^k fixiert

			
			//gehe alle moeglichen rechten Arme durch
			for ( size_t j = 1; (j+1) * (1<<k) - 1 <= sbEnd ; j++ ){
				raBegin = sbBegin + j * (1<<k);
				
				leftArms = kmpMatching(lce, sbBegin, sbEnd, raBegin, (1<<k));
				rRhoBegin = raBegin;
				lastPos = 0;
				p=0;

				//Laeufe und Periode bestimmen
				occs = new vector<int>(leftArms.size(),0);
				for ( size_t h = 0; h < leftArms.size() ; h++ ){

					if ( h == leftArms.size()-1 ){ //Grenzfall: letzter Wert, leftArms[h] == raBegin
						if ( leftArms[h] >= (int)lastPos ){ //liegt nicht im vorherigen run, also auf neuen run pruefen
							if ( p==0 || leftArms[h]+p >= lce->length || lcPrefix(lce, leftArms[h], leftArms[h]+p) < (int)p ){
								(*occs)[h] = -1;
								
							}
							else {
								(*occs)[h] = p;
							}
						}
					}
					else{
						if ( leftArms[h] >= (int)lastPos ){  //falls nicht im letzten run enthalten

							ph = leftArms[h+1]-leftArms[h];
							lcHelp = lcPrefix(lce, leftArms[h], leftArms[h]+ph);
							//im selben run wie y_rho
							if ( lcHelp >= ph && leftArms[h]+ph+lcHelp>=raBegin){
								(*occs)[h] = -2;
								rRhoBegin = min(rRhoBegin,(size_t)leftArms[h]);
								p = ph;
								lastPos = leftArms[h+1]+lcHelp;
							}
							//falls kein run, also single occ
							else if (lcHelp < ph || lcPrefix(lce, leftArms[h], raBegin) < (int)ph){
								(*occs)[h] = -1;

							}
							//falls im run, aber in anderem als y_rho
							else {
								(*occs)[h] = ph;
								p = ph;
								lastPos = leftArms[h+1]+lcHelp;
							}

						}
					}
		
				}

				
				for ( size_t i = 0; i < leftArms.size() && leftArms[i] < (int)raBegin; i++){
					// falls kein run
					if ( (*occs)[i]==-1 || (*occs)[ (*occs).size()-1 ]==-1 ){

						laBegin = leftArms[i];
						a = lcSuffix(lce , (laBegin -1), j*(1<<k) -1 );
						s = lcPrefix(lce, (laBegin + (1<<k)), (j+1)*(1<<k) );
						gappedRep = new alphaGappedRepeat(laBegin-a, raBegin-a, a+s+(1<<k));

						
						//auf negative Luecke pruefen und s. "To avoid duplicates (S.15)"
						//wenn gueltig: einfuegen
						if ( isValid(gappedRep,alpha) && (size_t)(1<<(k+1)) <= gappedRep->length && gappedRep->length < (size_t)(1<<(k+2)) 
							&& j*(1<<k)+sbBegin < gappedRep->rArm+(1<<k) && j*(1<<k)+sbBegin >= gappedRep->rArm ){
							//grList->push_back(gappedRep);	//auskommentiert fuer Messungen
						}
						delete gappedRep;
						gappedRep =0;
						
					}

					// periodischer Fall mit Fallunterscheidung
					else if ( (*occs)[i] > 0 && (*occs)[i]==(int)p ){

						laBegin = leftArms[i];
						rLLength = findLongestPeriod(lce, laBegin, p);
						rRLength = findLongestPeriod(lce, rRhoBegin, p);
						
						a = lcSuffix(lce, laBegin, rRhoBegin)-1;
						s = lcPrefix(lce, laBegin+rLLength, rRhoBegin+rRLength);

						//case a
						if ( rRLength > rLLength ){
							gappedRep = new alphaGappedRepeat(laBegin, rRhoBegin+rRLength-rLLength , rLLength+s);
							//pruefen und hinzufuegen
							if ( isValid(gappedRep,alpha) && (size_t)(1<<(k+1)) <= gappedRep->length && gappedRep->length < (size_t)(1<<(k+2)) 
								&& j*(1<<k)+sbBegin < gappedRep->rArm+(1<<k) && j*(1<<k)+sbBegin >= gappedRep->rArm ){
								//grList->push_back(gappedRep);	//auskommentiert fuer Messungen
							}
							delete gappedRep;
							gappedRep =0;
						}

						//case b
						helpLength = min(rLLength,rRLength)/p;
						if (rRLength > helpLength){
							for (size_t h = 0; (int)helpLength-(int)(p*h) > 0 &&
									rRhoBegin+rRLength-(helpLength-p*h)-laBegin <= alpha*(helpLength-p*h) ; h++){
								gappedRep = new alphaGappedRepeat( laBegin, rRhoBegin+rRLength-(helpLength-p*h) , 
																	helpLength-p*h );
								//pruefen und hinzufuegen
								if ( isValid(gappedRep,alpha) && (size_t)(1<<(k+1)) <= gappedRep->length && gappedRep->length < (size_t)(1<<(k+2))
									&&	j*(1<<k)+sbBegin < gappedRep->rArm+(1<<k) && j*(1<<k)+sbBegin >= gappedRep->rArm ){
									//grList->push_back(gappedRep);	//auskommentiert fuer Messungen
								}
								delete gappedRep;
								gappedRep =0;
							}
						}
						
						//case c
						for (size_t h = 1; rLLength+p*h<=rRLength && rRhoBegin+p*h-laBegin <= alpha*rLLength; h++){
							gappedRep = new alphaGappedRepeat(laBegin, rRhoBegin+p*h, rLLength);
							//pruefen und hinzufuegen
							if ( isValid(gappedRep,alpha) && (size_t)(1<<(k+1))<=gappedRep->length && gappedRep->length<(size_t)(1<<(k+2))
								&& j*(1<<k)+sbBegin < gappedRep->rArm+(1<<k) 
								&& j*(1<<k)+sbBegin >= gappedRep->rArm ){
								//grList->push_back(gappedRep);	//auskommentiert fuer Messungen
							}
							delete gappedRep;
							gappedRep =0;
						}

						//case d
						helpLength = min(rLLength,rRLength)/p;
						gappedRep = new alphaGappedRepeat(laBegin+rLLength-helpLength, rRhoBegin, helpLength);
						//pruefen und hinzufuegen
						if ( isValid(gappedRep,alpha) && (size_t)(1<<(k+1)) <= gappedRep->length && gappedRep->length < (size_t)(1<<(k+2)) 
							&& j*(1<<k)+sbBegin < gappedRep->rArm+(1<<k) && j*(1<<k)+sbBegin >= gappedRep->rArm ){
							//grList->push_back(gappedRep);	//auskommentiert fuer Messungen
						}
						delete gappedRep;
						gappedRep =0;
						
						//case e
						for (size_t h = 1; rRLength+p*h<=rLLength; h++){
							gappedRep = new alphaGappedRepeat(laBegin+p*h, rRhoBegin, rRLength);
							if ( isValid(gappedRep,alpha) && (size_t)(1<<(k+1))<=gappedRep->length && gappedRep->length<(size_t)(1<<(k+2))
								&& j*(1<<k)+sbBegin < gappedRep->rArm+(1<<k) 
								&& j*(1<<k)+sbBegin >= gappedRep->rArm ){
								//grList->push_back(gappedRep);	//auskommentiert fuer Messungen
							}
							delete gappedRep;
							gappedRep =0;
						}


						//case f
						if (rLLength > rRLength){
							gappedRep = new alphaGappedRepeat(laBegin+rLLength+-rRLength, rRhoBegin, rRLength+s);
							//pruefen und hinzufuegen
							if ( isValid(gappedRep,alpha) && (size_t)(1<<(k+1)) <= gappedRep->length && gappedRep->length < (size_t)(1<<(k+2))
								&& j*(1<<k)+sbBegin < gappedRep->rArm+(1<<k) && j*(1<<k)+sbBegin >= gappedRep->rArm ){
								//grList->push_back(gappedRep);	//auskommentiert fuer MessungenO
							}
							delete gappedRep;
							gappedRep =0;
						}

						//case g
						if ( a>0 ){ //sonst wird bedingung nicht erfüllt
							helpLength = min(rLLength,rRLength)/p;
							gappedRep = new alphaGappedRepeat(laBegin-a, rRhoBegin-a, helpLength);
						
							//pruefen und hinzufuegen
							if ( isValid(gappedRep,alpha) && (size_t)(1<<(k+1)) <= gappedRep->length && gappedRep->length < (size_t)(1<<(k+2))
								&& j*(1<<k)+sbBegin < gappedRep->rArm+(1<<k) && j*(1<<k)+sbBegin >= gappedRep->rArm ){
								//grList->push_back(gappedRep);	//auskommentiert fuer Messungen
							}
							delete gappedRep;
							gappedRep =0;
						}
					}

					//beide Arme im gleichen run
					else if ( (*occs)[i] == -2 ) {

						//rRhoBegin schon gegeben
						//p gegeben
						rRLength = findLongestPeriod(lce, rRhoBegin, p); //Laenge von rRho

						//Anzahl Occurances = rRLength/p
						//Squares
						for ( size_t h = 0; h <= rRLength -2; h++){ //Schleife ueber Position
							for ( size_t g = 1; (h+g) +g <= rRLength ; g++){ //Schleife ueber Laenge
								gappedRep = new alphaGappedRepeat(rRhoBegin+h, rRhoBegin+h+g, g);
								if ( isValid(gappedRep,alpha) && (size_t)(1<<(k+1)) <= gappedRep->length 
									&& gappedRep->length < (size_t)(1<<(k+2)) && j*(1<<k)+sbBegin < gappedRep->rArm+(1<<k) 
									&& j*(1<<k)+sbBegin >= gappedRep->rArm ){
									//grList->push_back(gappedRep);	//auskommentiert fuer Messungen
								}
								delete gappedRep;
								gappedRep =0;
							}
						}
						
						//lArm am Anfang, rArm am Ende vom run, unterschiedl. Laengen
						//for (size_t h = 1; h < ((float)rRLength/p)/2 ; h++){
						for (size_t h = rRLength/(p*(alpha+1)); h < ((float)rRLength/p)/2 ; h++){
							gappedRep = new alphaGappedRepeat(rRhoBegin, rRhoBegin+rRLength-h*p, h*p );
							if ( isValid(gappedRep,alpha) && (size_t)(1<<(k+1)) <= gappedRep->length 
								&& gappedRep->length < (size_t)(1<<(k+2)) && j*(1<<k)+sbBegin < gappedRep->rArm+(1<<k)
								&& j*(1<<k)+sbBegin >= gappedRep->rArm ){
								//grList->push_back(gappedRep);	//auskommentiert fuer Messungen
							}
							delete gappedRep;
							gappedRep =0;
							
						}
					
					}
					
				}
				delete occs;
				occs = 0;

			}	
		}
	}

	return 1;
}


//VARIANTE VON CALCSHORTARMS, DIE NUR SQUARES BERECHNET,
//WIRD NUR ZUM VERGLEICH VON MESSWERTEN VERWENDET
int calcSQUARES (lceDataStructure*& lce, vector<alphaGappedRepeat*> *grList, bool v2){
	size_t alpha = 1;
	alphaGappedRepeat* gappedRep;
	size_t n = lce->length;
	size_t sbBegin;				//Anfangsposition von Superblock
	size_t sbEnd;				//Endposition

	vector<int> leftArms;		//Liste der moeglichen linken Arme
	size_t raBegin;				//Anfangsposition des rechten Arms

	size_t p;					//Periode
	size_t rRLength;			//r_rho laenge (startet an raBegin)		//fuer Fall c und e
	bool oneblock = false;		//falls Text nur einen Superblock enthält
	
	vector<int> *occs;			// occs[h]=ph wenn leftArms[h] erstes Vorkommen des Runs ist und ph die Periode
								// occs[h]=0 wenn leftArms[h] spaeteres Vorkommen eines Runs ist
								// occs[h]=-1 wenn leftArms[h] single occurance ist
								// occs[h]=-2 wenn leftArms[h] im selben run wie y_rho ist
								// 				dann gilt p = leftArms[n-1] - leftArms[n-2] mit n=leftArms.size()
	size_t ph;					//vorlaeufige Periode p
	size_t lcHelp;
	size_t lastPos;				//speichert letzte Position des aktuellen runs
	size_t rRhoBegin;

	

	if ( constGamma*alpha <= n/log2(n) - constGamma -1 || v2 ){
		
		oneblock = true;
	}

	//fuer jeden Superblock
	for (size_t m = constGamma * alpha; m <= n/log2(n) - constGamma -1 || oneblock; m++){


		//Superblock bestimmen
		if (oneblock || v2){ //falls Text nur einen Superblock enthält (zu kurz oder v2)
			sbBegin = 0;
			sbEnd = n-1;
			oneblock = false; 	//damit Schleife ueber m nur einmal durchlaufen wird, wenn text zu kurz ist
			m = n/log2(n) +1; 
		}
		else if( (m - constGamma * alpha) * log2(n) < 1 ){
			sbBegin = 0;
			sbEnd = (m + constGamma + 1) * log2(n) - 1;
		}
		else{
			sbBegin = (m - constGamma * alpha) * log2(n);
			sbEnd = (m + constGamma + 1) * log2(n) - 1;
		}
		
		for ( size_t k = 0; k <= log2( constGamma * log2(n) ); k++){
			//Armlaenge auf 2^k fixiert

			
			//gehe alle moeglichen rechten Arme durch
			for ( size_t j = 1; (j+1) * (1<<k) - 1 <= sbEnd ; j++ ){
				raBegin = sbBegin + j * (1<<k);
				
				leftArms = kmpMatching(lce, sbBegin, sbEnd, raBegin, (1<<k));
				rRhoBegin = raBegin;
				lastPos = 0;
				p=0;


				occs = new vector<int>(leftArms.size(),0);
				for ( size_t h = 0; h < leftArms.size() ; h++ ){

					if ( h == leftArms.size()-1 ){ //Grenzfall: letzter Wert, leftArms[h] == raBegin
						if ( leftArms[h] >= (int)lastPos ){ //liegt nicht im vorherigen run, also auf neuen run pruefen
							if ( p==0 || leftArms[h]+p >= lce->length || lcPrefix(lce, leftArms[h], leftArms[h]+p) < (int)p ){
								(*occs)[h] = -1;
								
							}
							else {
								(*occs)[h] = p;
							}
						}
					}
					else{
						if ( leftArms[h] >= (int)lastPos ){  //falls nicht im letzten run enthalten

							ph = leftArms[h+1]-leftArms[h];
							lcHelp = lcPrefix(lce, leftArms[h], leftArms[h]+ph);
							//im selben run wie y_rho
							if ( lcHelp >= ph && leftArms[h]+ph+lcHelp>=raBegin){
								(*occs)[h] = -2;
								rRhoBegin = min(rRhoBegin,(size_t)leftArms[h]);
								p = ph;
								lastPos = leftArms[h+1]+lcHelp;
							}
							//falls kein run, also single occ
							else if (lcHelp < ph || lcPrefix(lce, leftArms[h], raBegin) < (int)ph){
								(*occs)[h] = -1;

							}
							//falls im run, aber in anderem als y_rho
							else {
								(*occs)[h] = ph;
								p = ph;
								lastPos = leftArms[h+1]+lcHelp;
							}

						}
					}
		
				}

				
				for ( size_t i = 0; i < leftArms.size() && leftArms[i] < (int)raBegin; i++){
					// falls kein run


					// periodischer Fall mit Fallunterscheidung


					//beide Arme im gleichen run
					if ( (*occs)[i] == -2 ) {

						//rRhoBegin schon gegeben
						//p gegeben
						rRLength = findLongestPeriod(lce, rRhoBegin, p); //Laenge von rRho

						//Anzahl Occurances = rRLength/p
						//Squares
						for ( size_t h = 0; h <= rRLength -2; h++){ //Schleife ueber Position
							for ( size_t g = 1; (h+g) +g <= rRLength ; g++){ //Schleife ueber Laenge
								gappedRep = new alphaGappedRepeat(rRhoBegin+h, rRhoBegin+h+g, g);
								if ( isValid(gappedRep,alpha) && (size_t)(1<<(k+1)) <= gappedRep->length 
									&& gappedRep->length < (size_t)(1<<(k+2)) && j*(1<<k)+sbBegin < gappedRep->rArm+(1<<k) 
									&& j*(1<<k)+sbBegin >= gappedRep->rArm ){
									//grList->push_back(gappedRep);	//auskommentiert fuer Messungen
								}
								delete gappedRep;
								gappedRep =0;
							}
						}
						
						//lArm am Anfang, rArm am Ende vom run, unterschiedl. Laengen

					
					}
					
				}
				delete occs;
				occs = 0;

			}	
		}
	}

	return 1;
}




//Berechnet Block-Repraesentation auf w
//template<typename vektor_type>
vector<int> calcBlockRep (lceDataStructure*& lce){
	size_t size = lce->length;
	double lgn = log2((double)size);

	vector<int> blockRep;
	blockRep.reserve(size/log2(size)+1);

	vector<int> clusters(size);
	clusters[0]=0;
	
	//Cluster des Suffix-Arrays bestimmen
	for(size_t i = 1; i < size; i++){
		if( lce->sa[i] > lgn ){
			clusters[i] = clusters[i-1];
		}
		else{
			clusters[i] = clusters[i-1]+1;
		}
	}

	//Bestimme Blockrepraesentation
	for(size_t i = 0; i <= size/log2(size)+1; i++){
		blockRep.push_back( clusters[lce->isa[i*log2(size)]] );
	}
	return blockRep;
}


//Berechnet k-Bloecke auf w fuer fixes k
//template<typename vektor_type>
vector<int> calcKBlock (lceDataStructure*& lce, size_t k){
	size_t size = lce->length;

	vector<int> blockRep;
	blockRep.reserve(size/ ( (1<<k) * log2(size)) +1 );

	
	//Bestimme k-Bloecke
	for(size_t i = 0; i <= size/ ( (1<<k) * log2(size)) +1; i++){
		blockRep.push_back( i*(1<<k)*log2(size) );
	}
	return blockRep;
}




//binaere Suche, durchlaeuft das Suffix-Array von der Blockrep w'
//sucht nach vorkommen von y in w'
//template<typename vektor_type>
int binarySearch(lceDataStructure*& lce, lceDataStructure*& lceBlock, size_t yStart, size_t yLength){
	vector<int> arms;
	size_t left = 0;
	size_t right = lceBlock->length; 
	size_t mid;
	//size_t neighbour;
	double lgn = log2((double)lce->length);
	while (left <= right){
		mid = left + (( right - left) / 2);
		if ( lcPrefix(lce,yStart,lce->text[mid*lgn]) >= (int)yLength){
			//Wert in Liste der Arme einfuegen
			return mid*lgn;
			
		}
		else {
			if ( lce->text[yStart] < lce->text[mid*lgn] ){
				right = mid - 1;
			}
			else {
				left = mid + 1;
			}
			
		}
	}
	return -1;
}



//Funktion zur Berechnung alpha-gapped repeats mit langen Armen
//Funktioniert nicht
template<typename lceDataStructure>
int calcLongArm (lceDataStructure*& lce, size_t alpha, vector<alphaGappedRepeat*> *grList){
	size_t n = lce->length;
	vector<int> blockRep = calcBlockRep(lce);
	stringstream ss;
    copy( blockRep.begin(), blockRep.end(), ostream_iterator<int>(ss, " ")); //blockRep in string umwandeln
    string blockString = ss.str();
	lceDataStructure* lceBlock = new lceDataStructure(blockString); //lce-ds auf string der block rep erstellen

	size_t yStart;	//Startpostion rechter Arm
	size_t yLength;
	int y2; //Startposition von y'
	//size_t suffix;
	//size_t prefix;
	alphaGappedRepeat* gappedRep;
	vector<int> leftArms;
	vector<int> kBlocks;
	for ( size_t k = 0; k <= (n/log2(n)); k++){
		kBlocks = calcKBlock(lce, k);
		for ( size_t i = 0; i < kBlocks.size(); i++){
			for ( size_t start = 0; start < log2(n); start++ ){
				// rechten Arm y fixieren
				yStart = kBlocks[i] + start;
				yLength = ((1<<(k-1))) * log2(n);
				
				y2 = binarySearch(lce, lceBlock, yStart, yLength); //TODO erzeugt Speicherzugriffsfehler
				if (y2 != -1){

					leftArms = kmpMatching (lceBlock, 0, lceBlock->length, y2, yLength);
					for ( size_t j = 0; j < leftArms.size(); j++ ){
						/*
						Arme verlaengern, pruefen und ausgeben wie in calcShortArm
						*/
					}
				}
			}
		}
	}
	
	return 0;
}


