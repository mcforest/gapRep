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
#include <time.h>

constexpr size_t kAnzahl = 50;

using namespace std;
using namespace sdsl;

int compare (const void * a, const void * b)
{
  if ( *(double*)a <  *(double*)b ) return -1;
  if ( *(double*)a == *(double*)b ) return 0;
//  if ( *(double*)a >  *(double*)b ) 
  return 1;
}


double median( double arr[]){
    //cout << arr[0] << endl;
    qsort (arr, kAnzahl, sizeof(double), compare);
    //cout << arr[0] << endl;
    return arr[kAnzahl/2];
}


int main(){

	const string text(1<<17, 'a');


    
    
    //Stringstats erstellen
    const StringStats stats = StringStats(std::move(text));
    
    //RMQ auf LCP-Array erstellen
    rmq_succinct_sada<> rmq(&stats.lcp);
    
    double timeLCPRMQ=0.0, timeLCP=0.0, timeNaiv=0.0, tstart;
    double arrLCPRMQ[kAnzahl], arrLCP[kAnzahl], arrNaiv[kAnzahl];
    
    //Dauer fuer LCPRMQ-Abfrage
    for(size_t k=0; k<kAnzahl; k++){
        tstart = clock();
        lcpRmqMin(stats, rmq, 0, 10000);    
        arrLCPRMQ[k] = 0.0 + clock() -tstart;
    }
    timeLCPRMQ = median(arrLCPRMQ);
    cout << "time lcprmq: " << timeLCPRMQ << endl;
    
    

    //Dauer fuer naiven Zeichenvergleich
    size_t x1 = 0;
    const char* textchar = text.c_str();
    size_t length = text.size();
    while(timeNaiv < timeLCPRMQ){                 //nur eine von beiden Schleifen
    //while( x1 <= 10000){
        x1++;
        for(size_t k=0; k<kAnzahl; k++){
            tstart = clock();    
            naivComp(textchar, 0, x1, length);
            arrNaiv[k] = 0.0 + clock() - tstart;
        }
        timeNaiv = median(arrNaiv);
//         if (x1 % 100 == 0){		//Ausgabe zum Messen
//             cout << x << "       "<< timeNaiv << endl;
//         }
    }
    cout << "time naiv: " << timeNaiv << endl;
    cout << "x: " << x1 << endl;
    
    
    
    //Dauer fuer iterative Minimumsuche im LCP-Array
    size_t y1 = 0;
    while(timeLCP < timeLCPRMQ){
	//while( y1 <= 10000){
        y1++;
        for(size_t k=0; k<kAnzahl; k++){
            tstart = clock();    
            lcpMin(stats, 0, y1);
            arrLCP[k] = 0.0 + clock() -tstart;   
        }
        timeLCP = median(arrLCP);
    }
//         if (y1 % 100 == 0){		//Ausgabe zum Messen
//             cout << y << "       "<< timeLCP << endl;
//         }
    cout << "time lcp: " << timeLCP << endl;
    cout << "y: " << y1 << endl;
    

    
    
    //in Datei schreiben
    ofstream fileout;
    fileout.open("grenzwerte.hpp"); 
    fileout << "constexpr size_t x = " << x1 << ";\nconstexpr size_t y = " << y << ";"; 
    fileout.close();
    
}
