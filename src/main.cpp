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

int main(int argc, char *argv[]){
    
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
	else {
		cerr << "Ungueltige Variante " << argv[3] << '\n';
	}
    
    return 0;
}



