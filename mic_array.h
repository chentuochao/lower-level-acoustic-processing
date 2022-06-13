#ifndef MIC_ARRAY_H
#define MIC_ARRAY_H

#include <thread>
#include <portaudio.h>
#include <vector>
#include "readerwriterqueue.h"
#include <sndfile.h>
#include <math.h>
#include <vector> 

typedef short SAVE_FORMAT;
#define SAVE_BIT 16
#define WIN_SIZE 16384
#define CHANNELS 8

struct audio_config_t{
    int buffer_size;
    int sampling_rate;
    int num_mic;
    int num_out;
    audio_config_t(int buffer_size, int n, int out, int sr):buffer_size(buffer_size), sampling_rate(sr), num_mic(n), num_out(out){}
};

class mic_array{
protected:
    bool mode;
    const int index_offset;
    const int wait_time;
    int output_wait;
    int buffer_index;
    double threshold; 

    int sending_index_of_buffer;
    int sending_index_in_buffer;
    bool sending_status;

    audio_config_t config;
    PaStream* stream;
    SNDFILE* out_file;
    
    std::vector<double>* current_window;
    double* preamble_pattern;
    unsigned long int preamble_len;

    SAVE_FORMAT* output_signal;
    //int output_frames;
    int out_len;
    int write_pointer;

    moodycamel::ReaderWriterQueue<SAVE_FORMAT*> queue_in;


    static int paCallback(const void* inputBuffer, void* outputBuffer, 
        unsigned long frames,const PaStreamCallbackTimeInfo* timeinfo,
        PaStreamCallbackFlags statusFlags,
        void* userData);
public:
    mic_array(audio_config_t config, int index_offset, int wait_time, double threshold):
        mode(0), index_offset(index_offset), wait_time(wait_time), output_wait(0), buffer_index(0), threshold(threshold), config(config), stream(NULL), preamble_len(0), out_len(0), write_pointer(0){

        }
    virtual bool init_mic_array(int output_wait, const char* sound_file,  const char* preamble_file);

    virtual bool start(const char* filename);
    virtual bool load_preamble(const char* filename);
    virtual bool load_speaker_file(const char* sound_file);

    virtual void run(int epoch, int record_len, int wait_len);
    virtual void run_detect(int epoch, int epoch_size);
    virtual void save_file(const char* out_path);
    
    virtual bool stop();
    ~mic_array(){
            if(stream){
                Pa_AbortStream(stream);
                Pa_CloseStream(stream);
            }
            sf_close(out_file);
            Pa_Terminate();
            //for(int i = 0; i < output_frames; i++){
            //    delete [] output_signal[i];
            //}
            delete []output_signal;
    
            delete [] current_window;
            delete [] preamble_pattern;
        }

};

#endif