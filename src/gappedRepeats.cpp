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


//vergleicht die Laenge der Arme ueber Zeichenvergleiche
//i und j sind Positionen im Text, length ist bereits die geforderte Mindestlaenge
int naivCompOld(const string &text, size_t i, size_t j, size_t length){
    boost::string_ref l = text.substr(i,length);   //left arm
    boost::string_ref r = text.substr(j,length);   //right arm
    if ( l != r ){
        return 0;
    }
    else {
        while ((text[i+length+1]==text[j+length+1]) && (i + length + 1 < j)){
            length++;
        }
        return length;
    }
}


//vergleicht die Laenge der Arme ueber Zeichenvergleiche
//nutzt word packing mit __int 128
//i und j sind Positionen im Text, length ist bereits die geforderte Mindestlaenge
int naivComp128(const string &text, size_t i, size_t j){
    bool gr = 1;
    __int128 l;
    __int128 r;
    size_t k;
    size_t length = 0;
    
    for (k=0; gr==1; k++){
        if (i + 816*k <= j && j + 16*k <= text.size()){
            l = text[i+16*k]<<(16*15) | text[i+16*k+1]<<(16*14)
                                        | text[i+16*k+2]<<(16*13)
                                        | text[i+16*k+3]<<(16*12)
                                        | text[i+16*k+4]<<(16*11)
                                        | text[i+16*k+5]<<(16*10)
                                        | text[i+16*k+6]<<(16*9)
                                        | text[i+16*k+7]<<(16*8)
                                        | text[i+16*k+8]<<(16*7)
                                        | text[i+16*k+9]<<(16*6)
                                        | text[i+16*k+10]<<(16*5)
                                        | text[i+16*k+11]<<(16*4)
                                        | text[i+16*k+12]<<(16*3)
                                        | text[i+16*k+13]<<(16*2)
                                        | text[i+16*k+14]<<(16*1)
                                        | text[i+16*k+15];
            r = text[j+16*k]<<(16*15) | text[j+16*k+1]<<(16*14)
                                        | text[j+16*k+2]<<(16*13)
                                        | text[j+16*k+3]<<(16*12)
                                        | text[j+16*k+4]<<(16*11)
                                        | text[j+16*k+5]<<(16*10)
                                        | text[j+16*k+6]<<(16*9)
                                        | text[j+16*k+7]<<(16*8)
                                        | text[j+16*k+8]<<(16*7)
                                        | text[j+16*k+9]<<(16*6)
                                        | text[j+16*k+10]<<(16*5)
                                        | text[j+16*k+11]<<(16*4)
                                        | text[j+16*k+12]<<(16*3)
                                        | text[j+16*k+13]<<(16*2)
                                        | text[j+16*k+14]<<(16*1)
                                        | text[j+16*k+15];
            if (l!=r){
                gr = 0;
                length = k*16;
                while ( text[i+length] == text[j+length] && i + length <= j && j + length <= text.size() ){
                    length++;
                }
            }
        }
        else {
            gr = 0;
            length = (k-1)*16;
            while ( text[i+length+1] == text[j+length+1] && i + length+1 <= j && j + length+1 <= text.size()){
                length++;
            }
        }
        
    }
    return length;
}



//vergleicht die Laenge der Arme ueber Zeichenvergleiche
//nutzt word packing mit uint64_t
//i und j sind Positionen im Text, length ist bereits die geforderte Mindestlaenge
int naivComp(const string &text, size_t i, size_t j){
    bool gr = 1;
    uint64_t l;
    uint64_t r;
    size_t k;
    size_t length = 0;
    
    for (k=0; gr==1; k++){
        if (i + 8*k <= j && j + 8*k <= text.size()){
            l = text[i+8*k]<<(8*7) | text[i+8*k+1]<<(8*6)
                                        | text[i+8*k+2]<<(8*5)
                                        | text[i+8*k+3]<<(8*4)
                                        | text[i+8*k+4]<<(8*3)
                                        | text[i+8*k+5]<<(8*2)
                                        | text[i+8*k+6]<<(8*1)
                                        | text[i+8*k+7];
            r = text[j+8*k]<<(8*7) | text[j+8*k+1]<<(8*6)
                                        | text[j+8*k+2]<<(8*5)
                                        | text[j+8*k+3]<<(8*4)
                                        | text[j+8*k+4]<<(8*3)
                                        | text[j+8*k+5]<<(8*2)
                                        | text[j+8*k+6]<<(8*1)
                                        | text[j+8*k+7];
            if (l!=r){
                gr = 0;
                length = k*8;
                while ( text[i+length] == text[j+length] && i + length <= j && j + length <= text.size() ){
                    length++;
                }
            }
        }
        else {
            gr = 0;
            length = (k-1)*8;
            while ( text[i+length+1] == text[j+length+1] && i + length+1 <= j && j + length+1 <= text.size()){
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

//vergleich mit LCP-Array und RMQ die Laenge der Arme
//i und j sind Positionen im S/LCP-Array
int lcpRmqMin(const StringStats &stats, rmq_succinct_sada<> &rmq, size_t i, size_t j){
    
    return stats.lcp[rmq(i+1,j)];
} 
