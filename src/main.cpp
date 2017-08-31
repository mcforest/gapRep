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

    
    //for(size_t i=0; i<= n-2; i++){
	for(size_t i=0; i < n; i++){
		//cout << stats.sa[i] << endl;
        //for(size_t j=i+1; j <= n-1 && (text[stats.sa[i]] == text[stats.sa[j]]); j++){
		for(size_t j=i+1; j < n && (text[stats.sa[i]] == text[stats.sa[j]]); j++){
            //if( stats.sa[i]==0 || stats.sa[j]==0 || text[stats.sa[i]-1] != text[stats.sa[j]-1] ){ //doppelte verhindern
                realX = abs(stats.sa[i]-stats.sa[j]);
                realY = j - i;
                if( realX < x && (realX * y) < (realY * x) ){         //wenig zeichenvergleiche noetig
                    l = min(stats.sa[i], stats.sa[j]);
                    r = max(stats.sa[i], stats.sa[j]);
                    length = ceil((r-l)/alpha);
                    length = naivComp(textchar, l, r, textlength);
                    if (length >= ceil((r-l)/(float)alpha)){
						if( stats.sa[i]==0 || stats.sa[j]==0 || text[stats.sa[i]-1] != text[stats.sa[j]-1] 
							|| l+length == r ){		//doppelte verhindern
                        	//cout << "("<<l<<","<<r<<","<<length<<")"<<endl;
							gappedRep = new alphaGappedRepeat(l,r,length);
							//grList.push_back(gappedRep); //auskommentiert fuer Messungen
							delete gappedRep;
							gappedRep=0;
						}
                    }
                }
                else if( realY < y){                         //kurzer Abstand im Suffix-Array
                    l = min(stats.sa[i], stats.sa[j]);
                    r = max(stats.sa[i], stats.sa[j]);
                    length = lcpMin( stats, i, j );
                    if (length >= ceil((r-l)/(float)alpha)){
                        if( stats.sa[i]==0 || stats.sa[j]==0 || text[stats.sa[i]-1] != text[stats.sa[j]-1] 
							|| l+length == r ){		//doppelte verhindern
                        	//cout << "("<<l<<","<<r<<","<<length<<")"<<endl;
							gappedRep = new alphaGappedRepeat(l,r,length);
							//grList.push_back(gappedRep);	//auskommentiert fuer Messungen
							delete gappedRep;
							gappedRep=0;
						}
                    }
                    
                }
                else{                                           //sonst
                    l = min(stats.sa[i], stats.sa[j]);
                    r = max(stats.sa[i], stats.sa[j]);
                    length = lcpRmqMin( stats, rmq, i, j);
                    if (length >= ceil((r-l)/(float)alpha)){
                        if( stats.sa[i]==0 || stats.sa[j]==0 || text[stats.sa[i]-1] != text[stats.sa[j]-1] 
							|| l+length == r ){		//doppelte verhindern
                        	//cout << "("<<l<<","<<r<<","<<length<<")"<<endl;
							gappedRep = new alphaGappedRepeat(l,r,length);
							//grList.push_back(gappedRep);	//auskommentiert fuer Messungen
							delete gappedRep;
							gappedRep=0;
						}
                    }
                }
                
            //}
        }
    }
    
	//cout << "alpha-gapped repeats:" << endl;
	//printGappedRepeat(grList);
	return 0;
}

//Berechnung von Gawrychowski, verwendet nur die Technik fuer kurze Arme und einen Superblock
int mainCalc2 ( string text, size_t alpha ){
	lceDataStructure* lce = new lceDataStructure(text);
	vector<alphaGappedRepeat*> grList;
	calc1Arm(lce, alpha, &grList);
	calcShortArm (lce, alpha, &grList, true);
	//cout << "alpha-gapped repeats:" << endl;
	//printGappedRepeat(grList);
	return 0;
}

//VARIANTE VON CALCSHORTARMS, DIE NUR SQUARES BERECHNET,
//WIRD NUR ZUM VERGLEICH VON MESSWERTEN VERWENDET
int mainCalcSquares ( string text, size_t alpha ){
	alpha = 1;
	lceDataStructure* lce = new lceDataStructure(text);
	vector<alphaGappedRepeat*> grList;
	calc1Arm(lce, alpha, &grList);
	calcSQUARES (lce, &grList, true);
	//cout << "alpha-gapped repeats:" << endl;
	//printGappedRepeat(grList);
	return 0;
}

//Berechnung von Gawrychowski, verwendet Blockrepraesentation
int mainCalc3 ( string text, size_t alpha ){
	size_t n0 = text.size();
	size_t n1 = n0;

	while ( n1/log2(n1) != (int)(n1/log2(n1)) ){
		n1++;
	}

	//string textX = text0;
	for( size_t i = n0; i < n1; i++ ){
		text = text + '\1';
	}


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

			string file0 = "data/beispiel.txt";
  	  		ifstream inFile0;
  			inFile0.open(file0);
    		stringstream strStream0;
    		strStream0 << inFile0.rdbuf();
    		const string text0 = strStream0.str();

			cout << text0 << endl;
			cout << text0.size() << endl;

			/*
			//TODO fuer lange Arme verwenden
			size_t n0 = text0.size();
			//size_t n0 = 17;
			size_t n1 = n0;

			//size_t n0 = 16;

			float logn = log2(n0);
			while ( n1/log2(n1) != (int)(n1/log2(n1)) ){
				n1++;
			}
			cout << n0 << " " << n1 << endl;

			string textX = text0;
			for( size_t i = n0; i < n1; i++ ){
				textX = textX + '\1';
			}

			cout << textX<< endl;
			cout << textX.size() << endl;

			/*
			string textX = text0+ '\1';

			cout << textX << endl;
			cout << textX.size() << endl;
			
			string text1 = "gabbdefabbx";
			string text2 = "aaaaaaaaaa";
			string text3 = "abac";
			*/
			
			
			lceDataStructure* lce = new lceDataStructure(text0);
			size_t lceabfrage;

			cout << "length: " << lce->length << endl;
			cout << "text: " << lce->text << endl;
			cout << "sa:  " << lce->sa << endl;
			cout << "isa: " << lce->isa << endl;
			cout << "lcp: " << lce->lcp << endl;
			cout << "mtext: " << lce->mtext << endl;
			cout << "msa:  " << lce->msa << endl;
			cout << "misa: " << lce->misa << endl;
			cout << "mlcp: " << lce->mlcp << endl;

			//vector<int> leftArms = kmpMatching(lce, 0, 36, 32 , 2);
			//cout << leftArms << endl;
			
			//cout << lce->text.size() << " " << lce->length  << endl;
			cout << "lce: " << lcPrefix(lce, 1, 7) << endl;
			
			/*
			vector<int> *vek;// = new vector<int>(5,0);
			vek =  new vector<int>(5,0);
			(*vek)[2]= 3;
			cout << *vek << endl;
			cout << "hi" << (*vek)[3] << endl;
			cout << (*vek).size() << endl;
			*/
			/*
			size_t test1;
			size_t test2;
			size_t test3;
			test1 = 5;
			test2 = 6;
			test3 = 1;

			cout << test1 - test2 << endl;
			cout << (int)test1 - test2*(int)test3 << endl;
			cout << test1 - (int)test2 << endl;
			cout << (int)test1 - (int)test2 << endl;
			*/


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
	else if ( variante == "vs" ){
    	mainCalcSquares ( text, alpha );
    }
	else {
		cerr << "Ungueltige Variante " << argv[3] << '\n';
	}
    
    return 0;
}



