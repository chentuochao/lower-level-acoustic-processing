#include <thread>
#include <chrono>
#include <sndfile.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <fftw3.h>
#include "mic_array.h"
using namespace std;

const int SR=48000;

const int BUFFER_SIZE=1024;
const int num_mics = 8;

const int record_length = 19200/960*5;
const int wait_length = 19200/960*30;

int main(int argc, const char* argv[]){
    if(argc<=1){
      cout<<"Argument is needed"<<endl;
    }
    const char* out_path=argv[1];
    cout << out_path << endl;
    // init the micrphone array structure
    audio_config_t layout(BUFFER_SIZE, num_mics, 2, SR);
    mic_array new_array(layout);
    char filename[] = "../long_120.wav";
    const char sound_file[] = "../chirp2.wav";

    bool success = new_array.init_mic_array(1, sound_file);
    if(! success) return 0;
    new_array.start(out_path);
    
    cout << "Begin recording ........" << endl;

    new_array.run(4, record_length, wait_length);

    //std::this_thread::sleep_for(std::chrono::seconds(40));
    new_array.stop();
    //new_array.save_file(filename);
    return 0;
}
