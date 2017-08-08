#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
//#include "strinalyze.cpp"
#include <iomanip>
#include <glog/logging.h>
#include "substring.hpp"
#include "checked_vector.hpp"
#include "gappedRepeats.cpp"




using namespace std;
using namespace sdsl;

//naive iterative Berechnung
int mainCalc1 ( string text, size_t alpha ){
//Stringstats erstellen
    const StringStats stats = StringStats(std::move(text));
    //Stringstatsausgabe als Test
    //stats.print(FLAGS_zeroindex);
    //cout << "gg " << stats.sa << endl;
    
    //RMQ auf LCP-Array erstellen
    rmq_succinct_sada<> rmq(&stats.lcp);
    
    //Ausgabe einer RMQ als Test
    //cout << "-----------------------------------" << endl;
    //cout << "RMQ: " << stats.lcp[rmq(5,7)] << endl;
    
    //Initialisierungen
    
    
    size_t n = stats.sa.size();
    size_t length;
    size_t l;
    size_t r;
    size_t realX;
    size_t realY;
    const char* textchar = text.c_str();
    size_t textlength = text.size();
	vector<alphaGappedRepeat*> grList;
	alphaGappedRepeat* gappedRep;
    
    for(size_t i=0; i<= n-2; i++){
        for(size_t j=i+1; j <= n-1 && (text[stats.sa[i]] == text[stats.sa[j]]); j++){
            if( stats.sa[i]==0 || stats.sa[j]==0 || text[stats.sa[i]-1] != text[stats.sa[j]-1] ){ //doppelte verhindern
                realX = abs(stats.sa[i]-stats.sa[j]);
                realY = j - i;
                if( realX < x && (realX * y) < (realY * x) ){         //wenig zeichenvergleiche noetig
                    l = min(stats.sa[i], stats.sa[j]);
                    r = max(stats.sa[i], stats.sa[j]);
                    length = ceil((r-l)/alpha);
                    length = naivComp(textchar, l, r, textlength);
                    if (length >= ceil((r-l)/alpha)){
                        //cout << "("<<l<<","<<r<<","<<length<<")"<<endl;
						gappedRep = new alphaGappedRepeat(l,r,length);
						grList.push_back(gappedRep);
                    }
                }
                else if( realY < y){                         //kurzer Abstand im Suffix-Array
                    l = min(stats.sa[i], stats.sa[j]);
                    r = max(stats.sa[i], stats.sa[j]);
                    length = lcpMin( stats, i, j );
                    if (length >= ceil((r-l)/alpha)){
                        //cout << "("<<l<<","<<r<<","<<length<<")"<<endl;
						gappedRep = new alphaGappedRepeat(l,r,length);
						grList.push_back(gappedRep);
                    }
                    
                }
                else{                                           //sonst
                    l = min(stats.sa[i], stats.sa[j]);
                    r = max(stats.sa[i], stats.sa[j]);
                    length = lcpRmqMin( stats, rmq, i, j);
                    if (length >= ceil((r-l)/alpha)){
                        //cout << "("<<l<<","<<r<<","<<length<<")"<<endl;
						gappedRep = new alphaGappedRepeat(l,r,length);
						grList.push_back(gappedRep);
                    }
                }
                
            }
        }
    }
    
	printGappedRepeat(grList);
	return 0;
}

//Berechnung von Gawrychowski, verwendet nur die Technik fuer kurze Arme und einen Superblock
int mainCalc2 ( string text, size_t alpha ){
	lceDataStructure* lce = new lceDataStructure(text);
	vector<alphaGappedRepeat*> grList;
	//calc1Arm(lce, alpha, &grList);
	calcShortArm (lce, alpha, &grList);
	cout << "alpha-gapped repeats:" << endl;
	printGappedRepeat(grList);
	return 0;
}

int mainCalc3 ( string text, size_t alpha ){
	lceDataStructure* lce = new lceDataStructure(text);
	vector<alphaGappedRepeat*> grList;
	//calc1Arm(lce, alpha, &grList);
	//calcShortArm (lce, alpha, &grList);
	//calcLongArm (lce, alpha, &grList);
	cout << "alpha-gapped repeats:" << endl;
	printGappedRepeat(grList);
	return 0;
}


int main(int argc, char *argv[]){
    
	//Test werden nur zum Implementieren verwendet, werden nach Fertigstellung entfernt
	if(argc == 2){
		istringstream ss0(argv[1]);
		string test;
    	if (!(ss0 >> test)){
        	cerr << "Ungueltiger String " << argv[1] << '\n';
        	return 1;
    	}
		if(test == "test"){
			cout << "Dies ist ein Test" << endl;

			string file0 = "data/beispiel4.txt";
  	  		ifstream inFile0;
  			inFile0.open(file0);
    		stringstream strStream0;
    		strStream0 << inFile0.rdbuf();
    		const string text0 = strStream0.str();

			cout << text0 << endl;

			string text1 = "gabbdefabbx";
			string text2 = "aaaaaaaaaa";
			string text3 = "abac";
			
			lceDataStructure* lce = new lceDataStructure(text0);
			//size_t lceabfrage;

			//cout << "length: " << lce->length << endl;
			//cout << "text: " << lce->text << endl;
			//cout << "sa:  " << lce->sa << endl;
			//cout << "isa: " << lce->isa << endl;
			//cout << "lcp: " << lce->lcp << endl;

			vector<int> leftArms = kmpMatching(lce, 0, 36, 32 , 2);
			cout << leftArms << endl;
			
			
			

			

		}
	}
	


    if(argc != 4){
		cout << "Bitte genau 3 Parameter fuer alpha, Text und Variante eingeben" << endl;
		cout << "Geben Sie z.B. ein: './build/src/gapRep2 2 beispiel.txt v1'" << endl;
        return 1;
    }
    
    //file lesen und in text speichern
    string file = argv[2];
    file = "data/" + file;
    ifstream inFile;
    inFile.open(file);
    stringstream strStream;
    strStream << inFile.rdbuf();
    string text = strStream.str();

	istringstream ss(argv[1]);
	size_t alpha;
    if (!(ss >> alpha)){
        cerr << "Ungueltige Zahl " << argv[1] << '\n';
        return 1;
    }
    
    //Textausgabe als Test
    //cout << text << endl;
    
	istringstream ss2(argv[3]);
	string variante;
    if (!(ss2 >> variante)){
        cerr << "Ungueltiger String " << argv[3] << '\n';
        return 1;
    }
	if ( variante == "v1" ){
    	mainCalc1 ( text, alpha );
    }
	else if ( variante == "v2" ){
    	mainCalc2 ( text, alpha );
    }
	else if ( variante == "v3" ){
    	mainCalc3 ( text, alpha );
    }
	else {
		cerr << "Ungueltige Variante " << argv[3] << '\n';
	}
    
    return 0;
}



