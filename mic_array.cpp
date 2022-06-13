#include <thread>
#include <chrono>
#include <portaudio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include "mic_array.h"
#include <assert.h> 
#include <fstream>
#include "util.h"
#include "matplotlib-cpp/matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

int mic_array::paCallback(const void* inputBuffer, void* outputBuffer, 
    unsigned long frames,const PaStreamCallbackTimeInfo* timeinfo,
    PaStreamCallbackFlags statusFlags,
    void* userData){

    mic_array* context=(mic_array*)userData;
    auto config=context->config;
    assert(frames==(unsigned long)(config.buffer_size));

    if(context->output_wait>0){
        SAVE_FORMAT* copy=new SAVE_FORMAT[config.num_mic*frames];
        memcpy(copy, inputBuffer, sizeof(SAVE_FORMAT)*config.num_mic*frames);
        delete [] copy;

        context->output_wait--;
        memset(outputBuffer, 0, sizeof(SAVE_FORMAT)*config.num_out*frames);

    }
    else{
        SAVE_FORMAT* copy=new SAVE_FORMAT[config.num_mic*frames];
        memcpy(copy, inputBuffer, sizeof(SAVE_FORMAT)*config.num_mic*frames);
        context->queue_in.enqueue(copy);

        if(statusFlags == 1) cout << "Warning! Underflow" << endl;
        else if(statusFlags == 2) cout << "Warning! Overflow" << endl;
        // feed the sound file to the speaker

        memset(outputBuffer, 0, sizeof(SAVE_FORMAT)*config.num_out*frames);

        
        if(context -> sending_status && context -> write_pointer <  context ->out_len - 1 && context ->buffer_index >= context -> sending_index_of_buffer){
            int end_pointer = context ->out_len - 1;
            int now_pointer = context -> write_pointer;
            int index_of_buffer = context -> sending_index_of_buffer;
            int index_in_buffer =  context -> sending_index_in_buffer;
            int channels = config.num_out;
            int frames2 = (int)frames;

            if( now_pointer == 0){
                if(context ->buffer_index == index_of_buffer){
                    //memcpy
                    //cout << "---1: " << now_pointer << ' '  << end_pointer << ' ' <<  index_of_buffer << ' ' << context -> buffer_index << endl;
                    memcpy(  (void*)((SAVE_FORMAT*)outputBuffer + channels*index_in_buffer), 
                        context->output_signal, sizeof(SAVE_FORMAT)*channels*(frames2 - index_in_buffer));
                    context -> write_pointer += (frames2 - index_in_buffer);
                }
                else{
                    cout << "Warning! miss the sending ddl" <<endl;
                    context -> sending_status = false;
                }
            }
            else{
                if(now_pointer + frames2 < end_pointer){
                    //memcpy
                    //cout << "---2: " << now_pointer << ' '  << end_pointer << ' ' <<  index_of_buffer << ' ' << context -> buffer_index << endl;
                    memcpy(outputBuffer, context->output_signal + channels*now_pointer, sizeof(SAVE_FORMAT)*channels*frames2);
                    context -> write_pointer += frames2;
                }
                else{
                    //memcpy
                    //cout << "---3: " << now_pointer << ' '  << end_pointer << ' ' <<  index_of_buffer << ' ' << context -> buffer_index << endl;
                    int remaining_length = end_pointer - now_pointer + 1;
                    memcpy(outputBuffer, context->output_signal + channels*now_pointer , sizeof(SAVE_FORMAT)*channels*remaining_length);
                    if(context ->mode == 1){                  
                        context -> write_pointer = end_pointer;
                        context -> sending_status = false;
                    }
                    else{
                        context -> write_pointer = 0;
                        context ->sending_index_of_buffer += 50;
                        context -> sending_status = true;                        
                    }
                }
            }
        }

        /*
        int out_index = context->output_index;
        out_index ++;
        if(out_index>= context->output_frames){
            out_index = 0;
        }
        //cout << out_index << ' ' << config.num_out*frames<< endl;
        memcpy(outputBuffer, context->output_signal[out_index], sizeof(SAVE_FORMAT)*config.num_out*frames);
        context->output_index  = out_index;
        */
        context -> buffer_index ++;
    }

    return 0;
}


bool mic_array::load_preamble(const char* filename){
    int max_max = 100000;
    unsigned long int i = 0;

    ifstream pre_file(filename, ios::binary | ios::in);
    short *pre = new short[max_max];
    preamble_len = 0;
    while(pre_file.read((char *) (pre + preamble_len),sizeof(short))){
        preamble_len ++;
    }
    preamble_pattern = new double[preamble_len];
    pre_file.close();

    for(i = 0; i < preamble_len ; ++i){
        preamble_pattern[i] = (double)pre[i]/30000.;
    }
    cout << "Load preamble length: " << preamble_len << endl;

    delete []pre;

    return true;
}

bool mic_array::load_speaker_file(const char* sound_file){
    SF_INFO in_file_info;
    SNDFILE* in_file=sf_open(sound_file, SFM_READ, &in_file_info);
    int channel = in_file_info.channels;
    int length0 = in_file_info.frames;
    int sr=in_file_info.samplerate;

    assert(channel==config.num_out);
    assert(sr==config.sampling_rate);

    cout << "successfully load the audio file with channels: " << channel << " and frame number: " << ' ' <<length0<<endl;
    out_len = length0;

    int total_num = out_len*channel;
    
    output_signal = new SAVE_FORMAT[total_num]; 
    memset(output_signal, 0, sizeof(SAVE_FORMAT)*total_num);
    
#if SAVE_BIT==32       
    sf_read_float(in_file, output_signal, total_num);
#else
    sf_read_short(in_file, output_signal, total_num);

#endif
    
    sf_close(in_file);

    return true;
}
/*
bool mic_array::load_speaker_file(const char* sound_file){
    int i = 0;
    SF_INFO in_file_info;
    SNDFILE* in_file=sf_open(sound_file, SFM_READ, &in_file_info);
    int channel = in_file_info.channels;
    int length0 = in_file_info.frames;
    int sr=in_file_info.samplerate;

    assert(channel==config.num_out);
    assert(sr==config.sampling_rate);
    int buffer_size = config.buffer_size*config.num_out;
    output_frames = (int)ceil((float)length0/(float)config.buffer_size);

    cout << "successfully load the audio file with channels: " << channel << " and frame number: " <<  length0 << ' ' <<output_frames<<endl;
    int write_pointer = length0*channel;
    
    output_signal = new SAVE_FORMAT*[output_frames]; 
    
    for(i = 0; i < output_frames; ++i){
        //cout << i << " write: " <<  write_pointer << ' ' << buffer_size <<endl;
        output_signal[i] = new SAVE_FORMAT[buffer_size];
        memset(output_signal[i], 0, sizeof(SAVE_FORMAT)*buffer_size);
        

#if SAVE_BIT==32       
        if(write_pointer > buffer_size){
            sf_read_float(in_file, output_signal[i], buffer_size);
        }
        else{
            sf_read_float(in_file, output_signal[i], write_pointer);
        }
#else
        if(write_pointer > buffer_size){
            sf_read_short(in_file, output_signal[i], buffer_size);
        }
        else{
            sf_read_short(in_file, output_signal[i], write_pointer);
        }
#endif
        write_pointer -= buffer_size;
    }
    sf_close(in_file);

    return true;
}
*/

bool mic_array::init_mic_array(int output_wait, const char* sound_file, const char* preamble_file){
    /*
        intialize the vector array for sliding window
    */
    int i = 0;
    int channel_num = config.num_mic;
    current_window = new vector<double>[channel_num];

    // preset the memory size for the vector, it will do not need to 
    for(i = 0; i < channel_num ; ++i){
        current_window[i].reserve(2*WIN_SIZE+1025);
    }

    load_preamble(preamble_file);
    load_speaker_file(sound_file);

    write_pointer = 0;
    sending_index_of_buffer= 39;
    sending_index_in_buffer = 0;
    if(mode == 1) // reply mode
      sending_status = false;
    else // sending type
      sending_status = true;
    
    /*
        intialize the input audio file and load preamble to the program  
    */
    //load_speaker_file(sound_file);

    /*
        intialize the portaudio stream
    */
    PaStreamParameters inputPara, outputPara;
    PaStream *stream=NULL;

    PaError err;
    err=Pa_Initialize();
    if(err!=paNoError) return false;
    printf("to be intial finish ...........!\n");

    // scan all the available audio devices 
    
    PaDeviceIndex valid_num = Pa_GetDeviceCount();
    //cout << valid_num << endl;
    PaDeviceIndex input_index = -1;
    PaDeviceIndex output_index = -1;

    for(i =0; i< valid_num ; ++i){
        auto info = Pa_GetDeviceInfo(i);
        //cout << i << ' '<< info->name << ' '  << info->maxInputChannels << endl;
        if(info->name[0]=='a' && info->name[1] =='c' && info->name[2] == '1' && info->name[3] == '0'&& info->name[4] == '8'){
            input_index = i;
        }
        if(info->name[0]=='a' && info->name[1] =='c' && info->name[2] == '1' && info->name[3] == '0'&& info->name[4] == '1'){
            output_index = i;
        }
    }
    if(input_index == -1 || output_index == -1){
        cout << "No seeed device found in audio" << endl;
        return false;
    }
    /*
    printf("aaaa: %d, out %d", inputPara.device, outputPara.device);
    auto inputInfo=Pa_GetDeviceInfo(Pa_GetDefaultInputDevice());
    auto outputInfo=Pa_GetDeviceInfo(Pa_GetDefaultOutputDevice());
    */
    
    auto inputInfo=Pa_GetDeviceInfo(input_index);
    auto outputInfo=Pa_GetDeviceInfo(output_index); 
    cout << input_index << ' ' << inputInfo->name << " input "  << inputInfo->maxInputChannels << endl;
    cout << output_index << ' ' << outputInfo->name << " output "  << outputInfo->maxOutputChannels << endl;

    inputPara.device=input_index;
    outputPara.device=output_index;
    inputPara.channelCount=config.num_mic;
    inputPara.suggestedLatency=inputInfo->defaultLowInputLatency;
    inputPara.hostApiSpecificStreamInfo=NULL;
    outputPara.channelCount=config.num_out;

#if SAVE_BIT==32
    cout << "setting to float during init" << endl;
    inputPara.sampleFormat=paFloat32;
    outputPara.sampleFormat=paFloat32;
#else
    inputPara.sampleFormat=paInt16;
    outputPara.sampleFormat=paInt16;
#endif
    
    outputPara.suggestedLatency=outputInfo->defaultLowOutputLatency;
    outputPara.hostApiSpecificStreamInfo=NULL;

    err=Pa_OpenStream(
        &stream,
        &inputPara,
        &outputPara,
        config.sampling_rate,
        config.buffer_size,
        paClipOff,
        paCallback,
        this);

    if(err!=paNoError) return false;
    this->stream=stream;

    this->output_wait=output_wait;
    return true;
}


bool mic_array::start(const char* filename){
    // open the port audio stream
    cout <<"try to open" << endl;
    auto err=Pa_StartStream(stream);
    if(err!=paNoError){
        cout << "stream open error" << endl;
        return err;
    }

    cout << "success configure stream" << endl;
    
    // configure the saving file
    SF_INFO out_file_info;
    out_file_info.channels=config.num_mic;
    out_file_info.samplerate=config.sampling_rate;

#if SAVE_BIT==32
    out_file_info.format=SF_FORMAT_WAV|SF_FORMAT_FLOAT;
#else
    out_file_info.format=SF_FORMAT_WAV|SF_FORMAT_PCM_16;
#endif
    out_file=sf_open(filename, SFM_WRITE, &out_file_info);

    return true;
}


void mic_array::run(int epoch, int record_len, int wait_len){
    SAVE_FORMAT* result=NULL;

    int buffer_size = config.buffer_size*config.num_mic;
    int total_num = 0;

    while(epoch > 0){
        if(queue_in.try_dequeue(result)){
            total_num ++;

#if SAVE_BIT==32
            sf_write_float(out_file, result, buffer_size);
#else
            sf_write_short(out_file, result, buffer_size);
#endif  
        
            if(total_num >= 500){
                sf_write_sync(out_file);
                epoch --;
                cout << "savefile, remaining items: " << queue_in.size_approx() << ' ' <<total_num<<endl;
                total_num = 0;
            }


            delete[] result;
        }
        else{
            Pa_Sleep(5);
        }
    }
}

void mic_array::run_detect(int epoch, int epoch_size){
    bool report_time = false;
    bool check_next = false;
    double previous_max = -1;
    int previous_max_index = -1;
    int sleep_loop = 0; // sleep to avoid recv the signal from our own speaker
    // intial the 
    /*
    ofstream  f1("../data/channel1.bin", ios::binary | ios::out);
    ofstream  f2("../data/channel2.bin", ios::binary | ios::out);
    ofstream  f3("../data/channel3.bin", ios::binary | ios::out);
    ofstream  f4("../data/channel4.bin", ios::binary | ios::out);
    */
    SAVE_FORMAT* result=NULL;
    int i = 0;
    int j = 0;
    int begin_idx = 0;
    bool previous_valid = false;

    int buffer_size = config.buffer_size*config.num_mic;
    int total_num = 0;
    int read_index = 0;
    int Loop_index = 0;

    while(epoch > 0){
        if(queue_in.try_dequeue(result)){
            total_num ++; 
            read_index ++;
#if SAVE_BIT==32
            sf_write_float(out_file, result, buffer_size);
#else
            sf_write_short(out_file, result, buffer_size);
#endif  
            if(total_num >= epoch_size){
                
                epoch --;
                //cout << "savefile, remaining items: " << queue_in.size_approx() << ' ' <<total_num<<endl;
                total_num = 0;
            }

            for(i = 0; i < config.buffer_size; ++i){
                for(j = 0; j < config.num_mic; ++j){
                    current_window[j].push_back((double) result[i*config.num_mic + j] );
                }
            }

            if(current_window[0].size() >= 2*WIN_SIZE){
                double** window_data = new double*[config.num_mic];

                //cout << "indexes: " << read_index << ' ' << buffer_index << ' ' << queue_in.size_approx() << endl;

                //if(report_time) auto t_start = std::chrono::high_resolution_clock::now();
                for(j = 0; j < config.num_mic; ++j){
                    assert (current_window[j].size() == 2*WIN_SIZE);
                    window_data[j] = new double[WIN_SIZE*2];
                    copy(current_window[j].begin(), current_window[j].end(), window_data[j]);

                    current_window[j].erase(current_window[j].begin(), current_window[j].begin()+WIN_SIZE);
                    assert (current_window[j].size() == WIN_SIZE);
                }


                //if(report_time) auto t_end1 = std::chrono::high_resolution_clock::now();
                /* process the recving data here */
                unsigned long int data_len = (unsigned long int ) WIN_SIZE*2;
                unsigned long int max_len = max(data_len, preamble_len);
                
                double *out = new double[max_len];

                if(!previous_valid){
                    begin_idx = calibrate_channel_4mic(window_data, config.num_mic, 2000, 10000);
                    previous_valid = true;
                    //auto t_start = std::chrono::high_resolution_clock::now();
                    //auto t_end = std::chrono::high_resolution_clock::now();
                    //double time1 = std::chrono::duration<double, std::milli>(t_end-t_start).count();
                    //cout << time1 <<endl; 
                }
                //if(report_time) auto t_end2 = std::chrono::high_resolution_clock::now();
                //f1.write((const char*)window_data[(begin_idx + 0)%config.num_mic], sizeof(double) * WIN_SIZE);
                //f2.write((const char*)window_data[(begin_idx + 1)%config.num_mic], sizeof(double) * WIN_SIZE);
                //f3.write((const char*)window_data[(begin_idx + 2)%config.num_mic], sizeof(double) * WIN_SIZE);
                //f4.write((const char*)window_data[(begin_idx + 3)%config.num_mic], sizeof(double) * WIN_SIZE);
                int channel_index = (begin_idx + config.num_mic - 1)%config.num_mic;

                xcorr(out, window_data[begin_idx], preamble_pattern, data_len, preamble_len, max_len);
                int max_idx = find_max(out, max_len);
                double current_max = out[max_idx]/10000;


                // checking the valid of the x-corr peak
                int preamble_begin_index = 0;
                bool valid_preamble = false;
                
                if(check_next){
                    if(current_max >= previous_max){
                        preamble_begin_index = max_idx;
                    }
                    else{
                        preamble_begin_index = previous_max_index - WIN_SIZE;
                    }
                    valid_preamble = true;
                    check_next = false;
                }
                else{
                    if(current_max > threshold){
                        //cout << "suppass threshold: " << current_max << ' ' << max_idx << endl;
                        previous_max = current_max;
                        previous_max_index = max_idx;
                        check_next = true;
                    }else{
                        check_next = false;
                    }
                }

                if(sleep_loop > 0 ) sleep_loop --;
                cout << "Loop: " << Loop_index << ' ' << read_index<< ' ' << current_max << ' ' << max_idx << endl;

                if(mode == 1 && valid_preamble && sending_status == false && sleep_loop <= 0){
                    int temp_index =  
                        (read_index*config.buffer_size - 2*WIN_SIZE)
                        + preamble_begin_index
                        + (wait_time - index_offset);

                    //cout << (read_index*config.buffer_size - 2*WIN_SIZE) + preamble_begin_index <<  ' ' << (wait_time - index_offset) << ' ' <<wait_time << ' ' <<index_offset<<' '<<temp_index <<endl;

                    sending_index_in_buffer = temp_index%(config.buffer_size);
                    sending_index_of_buffer = temp_index/(config.buffer_size) - 1;
                    write_pointer = 0;
                    sending_status = true;

                    sleep_loop = ceil(  (1.2*wait_time + out_len + preamble_begin_index )   / (double)WIN_SIZE ) + 1;
                    cout << "detect valid preamble and reply: " << sending_index_of_buffer << "----" << buffer_index  << ' ' << sleep_loop << endl;
                }

                //if(report_time) auto t_end3 = std::chrono::high_resolution_clock::now();
                for(j = 0; j < config.num_mic; ++j){
                    delete [] window_data[j];
                }
                delete []window_data;
                delete []out;

                sf_write_sync(out_file);
                Loop_index ++;
                /*
                if(report_time){
                    auto t_end4 = std::chrono::high_resolution_clock::now();
                    double time1 = std::chrono::duration<double, std::milli>(t_end1-t_start).count();
                    double time2 = std::chrono::duration<double, std::milli>(t_end2-t_end1).count();
                    double time3 = std::chrono::duration<double, std::milli>(t_end3-t_end2).count();
                    double time4 = std::chrono::duration<double, std::milli>(t_end4-t_end3).count();

                    cout << "Time consuming: " << time1 << ' ' << time2 <<' ' << time3 << ' ' << time4 <<endl; 
                }
                */

            }

            delete[] result;
        }
        else{
            Pa_Sleep(5);
        }
    }
    /*
    f1.close();
    f2.close();
    f3.close();
    f4.close();
    */
}


void mic_array::save_file(const char* out_path){
    SF_INFO out_file_info;
    out_file_info.channels=config.num_mic;
    out_file_info.samplerate=config.sampling_rate;

#if SAVE_BIT==32
    cout << "setting to float during output" << endl;
    out_file_info.format=SF_FORMAT_WAV|SF_FORMAT_FLOAT;
#else
    out_file_info.format=SF_FORMAT_WAV|SF_FORMAT_PCM_16;
#endif

    int buffer_size = config.buffer_size*config.num_mic;

    SNDFILE* out_file=sf_open(out_path, SFM_WRITE, &out_file_info);
    SAVE_FORMAT* result=NULL;
    while(queue_in.try_dequeue(result)){

#if SAVE_BIT==32
        sf_write_float(out_file, result, buffer_size);
#else
        sf_write_short(out_file, result, buffer_size);
#endif
        delete[] result;
        result=NULL;
    }
    sf_write_sync(out_file);
    sf_close(out_file);
}

bool mic_array::stop(){
    auto err=Pa_StopStream(stream);
    if(err==paNoError){
        //worker_stop=true;
        //working_thread.join();
        return true;
    } else {
        //TODO
        return false;
    }
}
