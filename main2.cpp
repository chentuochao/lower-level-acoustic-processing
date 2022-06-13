#include <thread>
#include <chrono>
#include <sndfile.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <fftw3.h>
#include "mic_array.h"

using namespace std;

const int SR=44100;


const int buffer_offset= 41520 - 38400; //41517
const int BUFFER_SIZE=1024;
const int wait_time = BUFFER_SIZE*60; //40
const int num_mics = 8;
const double threshold = 150;


int main(int argc, const char* argv[]){
    if(argc<=1){
      cout<<"Argument is needed"<<endl;
    }
    const char* out_path=argv[1];
    cout << out_path << endl;
    // init the micrphone array structure
    audio_config_t layout(BUFFER_SIZE, num_mics, 2, SR);
    mic_array new_array(layout, buffer_offset, wait_time, threshold);

    const char sound_file[] = "../raw/online1.wav";
    const char pre_file[] = "../raw/online1.bin";

    bool success = new_array.init_mic_array(10, sound_file, pre_file);
    if(! success) return 0;
    new_array.start(out_path);
    
    cout << "Begin recording ........" << endl;

    new_array.run_detect(5, 200);

    //std::this_thread::sleep_for(std::chrono::seconds(40));
    new_array.stop();
    //new_array.save_file(filename);
    return 0;
}
