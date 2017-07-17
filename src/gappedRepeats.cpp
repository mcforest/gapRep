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


const size_t constGamma = 2;

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
	//Mirror-Image von ttext beginnt ab text[size+1]

    //const std::string mtext;        //Mirror-Image
    //const vektor_type msa;
	//const vektor_type misa;
	//const vektor_type mlcp;
    //rmq_succinct_sada<> mrmq;
    
    lceDataStructure(const std::string& ttext) 
		: text(ttext + '\0' + string ( ttext.rbegin(), ttext.rend() ))
		, length(ttext.size())
		, sa(create_sa<vektor_type>(text, FLAGS_stripDollar))
		, isa(inverse<vektor_type>(sa))
		, lcp(create_lcp<vektor_type>(text, sa, isa))
        , rmq(&lcp)
        //, mtext(string ( text.rbegin(), text.rend() ))
        //, msa(create_sa<vektor_type>(mtext, FLAGS_stripDollar))
		//, misa(inverse<vektor_type>(msa))
		//, mlcp(create_lcp<vektor_type>(mtext, msa, misa))
        //, mrmq(&mlcp)
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


//Berechnet die Startpostionen der Cluster des Suffix-Arrays fuer die Berechnung
//von alpha-gapped repeats mit einer Armlaenge von 1
//template<typename vektor_type>
vector<int> calcClusterStarts(const std::string text, const vector<int> sa){
	int size = sa.size();
	vector<int> clusterStarts;
	clusterStarts.reserve(size);
	clusterStarts.push_back(0);
	int pos = 1;
	for(int i=1; i < size; i++){
		if(text[sa[i-1]] != text[sa[i]]){
			clusterStarts.push_back(i);
			pos++;
		}
	}


//	delete [] clusterStarts;
	return clusterStarts;
}

//Funktion zur Berechnung alpha-gapped repeats mit Armen der Laenge 1
//template<typename vektor_type>
int calc1Arm(lceDataStructure*& lce, float alpha, vector<alphaGappedRepeat*> *grList){
	vector<int> clusterStarts = calcClusterStarts(lce->text, lce->sa);
	size_t m = clusterStarts.size();
	//size_t n = lce->length;
	size_t n = lce->text.size();
	vector<int> sa = lce->sa;
	alphaGappedRepeat* gappedRep;

	//innerhalb der Cluster des SA nach Startpositionen sortieren
	for (size_t i = 0; i<m-1; i++){
		sort(sa.begin()+clusterStarts[i], sa.begin()+clusterStarts[i+1]);
	}

	//Suche der alpha-gapped repeats mit Armlaenge 1
	for (size_t i = 0; i<n-1; i++){
		//alpha*2, da unser suffix array groeßer ist als angenommen 
		//(da der text einmal vorw und einmal rueckw gespeichert ist
		for (size_t j = i+1; j<=i+alpha*2 && j<n; j++){				
				//cout << "l: " << sa[i] << " r: " << sa[j] << endl;
				if (lce->text[sa[i]]==lce->text[sa[j]]){
					//cout << "i: " << i << " j: " << j << endl;
					//cout << "hi" << endl;
					//gapRep gefunden
					if ( (i==0 || lce->text[sa[i]-1]!=lce->text[sa[j]-1]) &&
						 (j==n-1 || lce->text[sa[i]+1]!=lce->text[sa[j]+1])){
						gappedRep = new alphaGappedRepeat(sa[i], sa[j] , 1);
						if (gappedRep->lArm <= lce->length && gappedRep->rArm <= lce->length){
							if (gappedRep->lArm+(alpha*gappedRep->length) >= gappedRep->rArm ){
								grList->push_back(gappedRep);
							}
						}
					}
				}
				
				//sobald cluster verlassen wird nicht weiter pruefen
				else{
					j = i + alpha*2 + 1;
				}
				
		}
	}
	return 0;
}


//Lemma 17:
//Bekommt LCE-Datenstruktur, Startposition i eines Faktors und Periode p
//gibt Laenge des laengsten Faktors aus
size_t findLongestPeriod (lceDataStructure lce, size_t i, size_t p){
	size_t length = lce.lcp[lce.rmq(i+1,i+p)];
	return length;
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



vector<int> kmpMatching (lceDataStructure*& lce, size_t sbStart, size_t sbEnd, size_t raStart, size_t raLen){
	
	const char* text = lce->text.c_str();//TODO funktioniert das?
	string rArm = (string)(text+raStart);
    vector<int> lcs(raLen + 1, -1);
	vector<int> lArms;


	for(int i = 1; i <= raLen; i++)
	{
		int pos = lcs[i - 1];
		while(pos != -1 && rArm[pos] != rArm[i - 1]) pos = lcs[pos];
		lcs[i] = pos + 1;
	}

	int textPos = sbStart;
	int armPos = 0;
	while(textPos < sbEnd)
	{
		while(armPos != -1 && (armPos == raLen || rArm[armPos] != text[textPos])) armPos = lcs[armPos];
		armPos++;
		textPos++;
		if(armPos == raLen) lArms.push_back(textPos - raLen);
	}
	
	return lArms;
}



/*
vector<int> kmpMatching (lceDataStructure lce, size_t sbStart, size_t sbEnd, size_t raStart, size_t raLen){
	vector<int> x;
	return x;
}
*/


//TODO um "schnellere" LCE-Anfragen erweitern
// gibt longest common prefix von 2 Woertern aus, die an Position i und j beginnen
int lcPrefix(lceDataStructure*& lce, size_t i, size_t j){    
    return lce->lcp[ lce->rmq( lce->isa[i]+1, lce->isa[j] ) ];
} 


//TODO um "schnellere" LCE-Anfragen erweitern
// gibt longest common suffix von 2 Woertern aus, die an Position i und j beginnen
int lcSuffix(lceDataStructure*& lce, size_t i, size_t j){  
	size_t n = lce->length;  
    return lce->lcp[ lce->rmq( lce->isa[i+n+1]+1, lce->isa[j+n+1] ) ];
} 

//Funktion zur Berechnung alpha-gapped repeats mit kurzen Armen
//schnellere Berechnung fuer Perioden bisher nicht enthalten
//template<typename vektor_type>
int calcShortArm (lceDataStructure*& lce, float alpha, vector<alphaGappedRepeat*> *grList){
	size_t n = lce->length;
	size_t sbBegin;				//Anfangsposition von Superblock
	size_t sbEnd;				//Endposition
	vector<int> leftArms;
	size_t raBegin;				//Anfangsposition des rechten Arms
	size_t raEnd;				//Endposition
	alphaGappedRepeat* gappedRep;
	int a;						//Suffix des Arms
	int s;						//Prefix des Arms

	//TODO naechste zeile entfernen
	//vector<int> testvar = kmpMatching(lce, sbBegin, sbEnd, raBegin, 0);
	
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
			for ( size_t j = 1; (j+1) * pow(2,k) - 1 <= sbEnd ; j++ ){
				raBegin = j * pow(2,k);
				raEnd = (j+1) * pow(2,k) - 1;
				
				//TODO leftArms = Lister der moeglichen linken Arme
				leftArms = kmpMatching(lce, sbBegin, sbEnd, raBegin, pow(2,k));
				//leftArms = kmpMatching();
				for ( size_t i = 0; i < leftArms.size(); i++ ){

					a = lcSuffix(lce, j*pow(2,k) -1, (leftArms[i] -1) ); //TODO richtig um 1 verringert?

					s = lcPrefix(lce, (j+1)+pow(2,k) , (leftArms[i] + pow(2,k)) );

					
					//TODO muessen Positionen noch um 1 verringert werden? Laenge richtig?
					gappedRep = new alphaGappedRepeat(leftArms[i]-a, raBegin-a, a+s+pow(2,k));

					//auf negative Luecke pruefen und s. "To avoid duplicates (S.15)"
					//wenn gueltig: einfuegen
					if( gappedRep->lArm < gappedRep->rArm && j*pow(2,k)+sbBegin <= gappedRep->rArm+pow(2,k) ){
						grList->push_back(gappedRep);
					}
					//TODO Fallunterscheidung fuer periodischen Fall

				}
			}	
		}
	}

	return 1;
}







//Funktion erweitert gapped repeat maximal nach links und rechts
//bekommt LCE-Datenstruktur, Start- und Endposition des Superblocks, 
//		  k (Armlaenge 2^k), j, "Startpositionen" beider y-Arme
/*
alphaGappedRepeat* extendArms (lceDataStructure lce, size_t sbBegin, size_t sbEnd, 
								size_t k, size_t j, size_t raBegin, size_t laBegin) {
	alphaGappedRepeat* rep = 0; //new alphaGappedRepeat(0,0,0);

	size_t p; //TODO Periode bestimmen / auslesen
	//bestimme rl und rr mit Lemma 17 (nur fuer periodischen Fall)
	size_t rl = findLongestPeriod(lce, laBegin, p);
	size_t rr = findLongestPeriod(lce, raBegin, p);
	
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
		//xxxxxxxxxx//blockRep.push_back(lce.isa[i*log2(size)]);
	}
	return blockRep;
}


//Berechnet k-Bloecke auf w fuer fixes k
//template<typename vektor_type>
vector<int> calcKBlock (lceDataStructure*& lce, size_t k){
	size_t size = lce->length;
	double lgn = log2((double)size);

	//vector<int> blockRep;
	//blockRep.reserve(size/ ( (pow(2,k)) * log2(size)) +1 );

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

	/* folgendes ueberfluessig, es werden nur startpositionen benoetigt
	//Bestimme Blockrepraesentation
	for(size_t i = 0; i <= size/ ( (pow(2,k)) * log2(size)) +1; i++){
		blockRep.push_back( clusters[lce->isa[i*(2^k)*log2(size)]] );
		//xxxxxxxxxxxxxxx//blockRep.push_back(lce.isa[i*log2(size)]);
	}
	return blockRep;
	*/
	return clusters;
}




//binaere Suche, durchlaeuft das Suffix-Array von der Blockrep w'
//sucht nach vorkommen von y in w'
//template<typename vektor_type>
int binarySearch(lceDataStructure*& lce, lceDataStructure*& lceBlock, size_t yStart, size_t yLength){
	vector<int> arms;
	size_t left = 0;
	size_t right = lceBlock->length; //TODO richtig runden?
	size_t mid;
	//size_t neighbour;
	double lgn = log2((double)lce->length); //TODO runden?
	while (left <= right){
		mid = left + (( right - left) / 2);
		if ( lcPrefix(lce,yStart,lce->text[mid*lgn]) >= yLength){
			//Wert in Liste der Arme einfuegen
			return mid*lgn;
			//arms.push_back(mid*lgn);
			/*
			//Nachbarn untersuchen
			neighbour = mid - 1;
			while (lcPrefix(lce,yStart,lce->text[neighbour*lgn]) >= yLength && neighbour >= 0){
				arms.push_back(neighbour*lgn);
				neighbour--;
			}
			neighbour = mid + 1;
			while (lcPrefix(lce,yStart,lce->text[neighbour*lgn]) >= yLength && neighbour <= lceBlock->length){
				arms.push_back(neighbour*lgn);
				neighbour++;
			}
			*/
			
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



//Funktion zur Berechnung alpha-gapped repeats mit kurzen Armen
template<typename lceDataStructure>
int calcLongArm (lceDataStructure*& lce, float alpha, vector<alphaGappedRepeat*> *grList){
	size_t n = lce->length;
	vector<int> blockRep = calcBlockRep(lce);
	stringstream ss;
    copy( blockRep.begin(), blockRep.end(), ostream_iterator<int>(ss, " ")); //blockRep in string umwandeln
    string blockString = ss.str();
	lceDataStructure* lceBlock = new lceDataStructure(blockString); //lce-ds auf string der block rep erstellen

	size_t yStart;
	size_t yLength;
	size_t y2; //Startposition von y'
	size_t suffix;
	size_t prefix;
	alphaGappedRepeat* gappedRep;
	vector<int> rightArms;
	vector<int> kBlocks;
	for ( size_t k = 0; k <= (n/log2(n)); k++){
		kBlocks = calcKBlock(lce, k);
		for ( size_t i = 0; i < kBlocks.size(); i++){
			for ( size_t start = 0; start < log2(n); start++ ){
				// linken Arm y fixieren
				yStart = kBlocks[i] + start;
				yLength = (pow(2,(k-1))) * log2(n);
				
				y2 = binarySearch(lce, lceBlock, yStart, yLength);
				if (y2 != -1){
					//TODO yLenght korrigieren, es ist nicht die eigentliche laenge gefragt, sondern die in der blockrep
					rightArms = kmpMatching (lceBlock, 0, lceBlock->length, y2, yLength);
					for ( size_t j = 0; j < rightArms.size(); j++ ){
						prefix = lcPrefix(lce, yStart, rightArms[j]);
						suffix = lcSuffix(lce, yStart, rightArms[j]);
						//TODO muessen Positionen noch um 1 verringert werden? Laenge richtig?
						gappedRep = new alphaGappedRepeat(yStart-prefix, rightArms[j]-prefix, prefix+suffix);
						//auf Gueltigkeit pruefen
						//TODO auf alpha pruefen und ob der rechte Arm im richtigen kBlock z beginnt
						if ( gappedRep->lArm + gappedRep->length < gappedRep->rArm && pow(2,(k+1)) <= gappedRep->length < pow(2,(k+2))){
							grList->push_back(gappedRep);
					}
					//TODO if ( 2^(k+1) <= yLength < 2^(k+2) ){hinzufuegen}
					}
				}
			}
		}
	}
	
	return 0;
}


