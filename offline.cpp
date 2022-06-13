#include <thread>
#include <chrono>
#include <sndfile.h>
#include <iostream>
#include <fstream>
#include <vector>
#include<math.h>
#include <complex> 
#include <fftw3.h>
#include "mic_array.h"
#include "util.h"
#include "matplotlib-cpp/matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;



int main(void){
    unsigned long int i = 0;
    ifstream tx_file("../data/tx.bin", ios::binary | ios::in);
    ifstream rx_file("../data/rx2.bin", ios::binary | ios::in);
    
    /*
    sig_file.seekg (0, ios::end);
    int length_sig = sig_file.tellg();
    sig_file.seekg (0, ios::beg);

    pre_file.seekg (0, ios::end);
    int length_pre = pre_file.tellg();
    pre_file.seekg (0, ios::beg);
    */
    int max_max = 100000;

    unsigned long int Ns = 0; //(unsigned long int) length_sig * sizeof(char)/sizeof(short);
    unsigned long int recv_len = 0; //(unsigned long int) length_pre * sizeof(char)/sizeof(short);

    double*tx = new double[max_max];
    double*rx = new double[max_max];

    while(rx_file.read((char *) (rx + recv_len), sizeof(double))){
        //cout << rx[len1] << ' ';
        recv_len ++;
    }
    cout << endl;
    
    while(tx_file.read((char *) (tx + Ns), sizeof(double))){
        //cout << rx[len1] << ' ';
        Ns ++;
    }
    cout << recv_len <<' ' <<Ns <<endl; 
    
    rx_file.close();
    tx_file.close();
    
    unsigned long int N0 = 480; 
    int fs = 44100;
    int BW1 = 1000;
    int BW2 = 5000;
    
    
    
    int DIVIDE_FACTOR = 2;
    
    double* h = new double[Ns];
    
    cout << "--------" << endl;
    channel_estimation_freq_multiple(h,tx, Ns, rx,  recv_len,N0,  BW1,  BW2,  fs);
    
    cout << "--------" << endl;
    
    for( i=0; i<Ns; ++i) {
        cout << h[i] << ' ';
    }
    cout << endl;


    delete []h;
    /*

    std::vector<double> y(Ns);
    for( i=0; i<Ns; ++i) {
        y.at(i) = h[i];
    }
    cout << endl;
    plt::plot(y);
    plt::show();

    plt::plot(y);
    plt::show();

    
    double in1[10] = {1, 2, 1, 2, 3, 4, 5, 6, 7, 2};
    double in2[10] = {1, 2, 3, 4, 5, 6, 7, 2, 1, 2};
    unsigned long int len1 = 10;
    unsigned long int len2 = 10;
    unsigned long int max_len = 10;
    

    double *out = new double[max_len];
    
    xcorr(out, in1, in2, len1, len2, max_len);
    

    
    for(i = 0 ; i < max_len; ++i){
        cout << out[i] << ' ';
    }
    cout << endl;

    int max_idx = find_max(out, max_len);
    cout << max_idx << endl;
    
    std::vector<double> y(max_len);
    for( i=0; i<max_len; ++i) {
        y.at(i) = out[i];
    }


    delete out;
    */
    delete []tx;
    delete []rx;

    return(0);
}
