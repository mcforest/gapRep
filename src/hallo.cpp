#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>

//kompilieren mit          g++ test.cpp -Wall -g -o test
//-I strinalyze
using namespace std;
int main(int argc, char *argv[]){
    
    
    
    //in Datei schreiben
    ofstream fileout;
    fileout.open("pizza.hpp"); //hier: "src/grenzwerte.hpp"
    fileout << "hallowelt";             //hier grenzwerte eintragen
    fileout.close();
    
    return 0;
}
 
