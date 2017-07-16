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
#include "../build/src/grenzwerte.hpp"



using namespace std;
using namespace sdsl;

int mainCalc1 ( string text, float alpha ){
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
                        cout << "("<<l<<","<<r<<","<<length<<")"<<endl;
                    }
                }
                else if( realY < y){                         //kurzer Abstand im Suffix-Array
                    l = min(stats.sa[i], stats.sa[j]);
                    r = max(stats.sa[i], stats.sa[j]);
                    length = lcpMin( stats, i, j );
                    if (length >= ceil((r-l)/alpha)){
                        cout << "("<<l<<","<<r<<","<<length<<")"<<endl;
                    }
                    
                }
                else{                                           //sonst
                    l = min(stats.sa[i], stats.sa[j]);
                    r = max(stats.sa[i], stats.sa[j]);
                    length = lcpRmqMin( stats, rmq, i, j);
                    if (length >= ceil((r-l)/alpha)){
                        cout << "("<<l<<","<<r<<","<<length<<")"<<endl;
                    }
                }
                
            }
        }
    }
    
	return 0;
}


int mainCalc2 ( string text, float alpha ){
	lceDataStructure* lce = new lceDataStructure(text);
	vector<alphaGappedRepeat*> grList;
	float test = 2.0;
	calc1Arm(lce, test, &grList);
	calcShortArm (lce, alpha, &grList);
	//calcLongArm (lce, alpha, &grList);
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

			string text1 = "gabbdefabbx";
			
			lceDataStructure* lce = new lceDataStructure(text0);
			
			cout << "text: " << lce->text << endl;
			cout << "length: " << lce->length << endl;
			cout << "sa: " << lce->sa << endl;
			cout << "isa: " << lce->isa << endl;
			cout << "lcp: " << lce->lcp << endl;
			cout << "length (text+mtext): " << lce->text.size() << endl;
			cout << lce->text[10] << lce->text[13] << lce->text[14] << endl;
			//cout << "mtext: " << lce->mtext << endl;
			//cout << "msa: " << lce->msa << endl;
			//cout << "misa: " << lce->misa << endl;
			//cout << "mlcp: " << lce->mlcp << endl;
			

			/*
			vector<int> vec (4);
			cout << vec.capacity() << endl;
			vec.reserve(28);
			cout << vec.capacity() << endl;
			vec[0] = 1;
			vec.push_back (5);
			cout << vec << endl;
			cout << vec.size() << endl;
			vec.push_back (5);
			cout << vec << endl;
			cout << vec.size() << endl;
			*/
			/*
			alphaGappedRepeat* agr1 = new alphaGappedRepeat(5,29,3);
			alphaGappedRepeat* agr2 = new alphaGappedRepeat(7,10,2);
			alphaGappedRepeat* agr3 = new alphaGappedRepeat(12,25,9);
			alphaGappedRepeat* agr4 = new alphaGappedRepeat(15,22,5);
			vector<alphaGappedRepeat*> vec;
			vec.push_back(agr1);
			vec.push_back(agr2);
			vec.push_back(agr3);
			vec.push_back(agr4);
			printGappedRepeat(vec);
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
	float alpha;
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
	else {
		cerr << "Ungueltige Variante " << argv[3] << '\n';
	}
    
    return 0;
}



