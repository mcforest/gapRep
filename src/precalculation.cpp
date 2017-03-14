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

#define ANZAHL 50

using namespace std;
using namespace sdsl;

int compare (const void * a, const void * b)
{
  if ( *(double*)a <  *(double*)b ) return -1;
  if ( *(double*)a == *(double*)b ) return 0;
  if ( *(double*)a >  *(double*)b ) return 1;
}


double median( double arr[]){
    //cout << arr[0] << endl;
    qsort (arr, ANZAHL, sizeof(double), compare);
    //cout << arr[0] << endl;
    return arr[ANZAHL/2];
}


int main(){
    
    //file lesen und in text speichern
    string file = "a.txt";
    file = "data/" + file;
    ifstream inFile;
    inFile.open(file);
    stringstream strStream;
    strStream << inFile.rdbuf();
    string text = strStream.str();
    
    
    //Stringstats erstellen
    const StringStats stats = StringStats(std::move(text));
    
    //RMQ auf LCP-Array erstellen
    rmq_succinct_sada<> rmq(&stats.lcp);
    
    int test = ANZAHL;                      //Anzahl der Tests aus denen Median bestimmt wird
    double timeLCPRMQ=0.0, timeLCP=0.0, timeNaiv=0.0, tstart;
    double arrLCPRMQ[test], arrLCP[test], arrNaiv[test];
    int k;
    
    
    
    //Dauer fuer LCPRMQ-Abfrage
    for(k=0; k<test; k++){
        tstart = clock();
        lcpRmqMin(stats, rmq, 0, 10000);    
        arrLCPRMQ[k] = 0.0 + clock() -tstart;
    }
    timeLCPRMQ = median(arrLCPRMQ);
    cout << "time lcprmq: " << timeLCPRMQ << endl;
      
    
    //Dauer fuer naiven Zeichenvergleich
    size_t x = 0;
    while(timeNaiv < timeLCPRMQ){                 //nur eine von beiden Schleifen
    //while( x <= 2000){
        x++;
        for(k=0; k<test; k++){
            tstart = clock();    
            naivComp(text, 0, x);
            arrNaiv[k] = 0.0 + clock() - tstart;
        }
        timeNaiv = median(arrNaiv);
//         if (x % 100 == 0){
//             cout << x << "       "<< timeNaiv << endl;
//         }
    }
    cout << "time naiv: " << timeNaiv << endl;
    cout << "x: " << x << endl;
    
    
    
    //Dauer fuer iterative Minimumsuche im LCP-Array
    size_t y = 0;
    while(timeLCP < timeLCPRMQ){
        y++;
        for(k=0; k<test; k++){
            tstart = clock();    
            lcpMin(stats, 0, y);
            arrLCP[k] = 0.0 + clock() -tstart;   
        }
        timeLCP = median(arrLCP);
    }
    cout << "time lcp: " << timeLCP << endl;
    cout << "y: " << y << endl;
    

    
    
    //in Datei schreiben
    ofstream fileout;
    fileout.open("pizza.hpp"); 
    fileout << "constexpr size_t x = " << x << ";\nconstexpr size_t y = " << y << ";"; 
    fileout.close();
    
}
