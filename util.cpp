#include <iostream>
#include <string.h>
#include <math.h>
#include "util.h"
using namespace std;


void multiple(fftw_complex * cross_vector, fftw_complex * in1, fftw_complex * in2, unsigned long int len_in){
    for(unsigned long int i = 0; i < len_in; ++i){
        cross_vector[i][0] = ((in1[i][0] * in2[i][0]) - (in1[i][1] * (-in2[i][1])) ) / len_in;
		cross_vector[i][1] = ((in1[i][0] * (-in2[i][1])) + (in1[i][1] * in2[i][0]) ) / len_in;
    }
}

void xcorr(double* corr_result, double *sig, double *preamble, unsigned long int sig_len, unsigned long int preamble_len, unsigned long int max_len){  
    unsigned long int i = 0;
    double* in1 = new double[max_len]; 
    double* in2 = new double[max_len];
    
    memset(in1, 0, sizeof(double)*max_len);
    memset(in2, 0, sizeof(double)*max_len);
    memcpy(in1, sig, sizeof(double)*sig_len);
    memcpy(in2, preamble, sizeof(double)*preamble_len);
    
    //for(i = 0; i < max_len; ++i){
    //    cout << in2[i] << ' ';
    //}
    //cout << endl;

    fftw_complex *out1=(fftw_complex*)fftw_malloc(max_len*sizeof(fftw_complex));
    if(out1==NULL){fprintf(stderr,"malloc failed\n");exit(1);} 
    fftw_complex *out2=(fftw_complex*)fftw_malloc(max_len*sizeof(fftw_complex));
    if(out2==NULL){fprintf(stderr,"malloc failed\n");exit(1);}
    
    
    fftw_plan p1=fftw_plan_dft_r2c_1d(max_len, in1, out1, FFTW_ESTIMATE);
    if (p1==NULL){fprintf(stderr,"plan creation failed\n");exit(1);}

    fftw_plan p2=fftw_plan_dft_r2c_1d(max_len, in2, out2, FFTW_ESTIMATE);
    if (p2==NULL){fprintf(stderr,"plan creation failed\n");exit(1);}
    
    fftw_execute(p1);
    fftw_execute(p2);

    //for(i = 0; i < 1200; ++i){
    //    cout << out2[i][0] << ' ' << out2[i][1] << endl;
    //}
    //cout << endl;

    fftw_complex *cross_vector = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*max_len);
    multiple(cross_vector, out1, out2, max_len);

    fftw_plan ifft=fftw_plan_dft_c2r_1d(max_len, cross_vector, corr_result, FFTW_ESTIMATE);
    fftw_execute(ifft);


    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(ifft);

    fftw_free(out1);
    fftw_free(out2);
    fftw_free(cross_vector);

    delete [] in1;
    delete [] in2;
}


int find_max(double* in, unsigned long int len_in){
    int maxPosition = max_element(in,in+len_in) - in; 
    return maxPosition;
}

int find_min(double* in, unsigned long int len_in){
    int minPosition = min_element(in,in+len_in) - in; 
    return minPosition;
}

int calibrate_channel_4mic(double** data, int channel_num, int begin_id, int end_id){
    double *channel_mean = new double[channel_num]; 

    memset(channel_mean, 0,  sizeof(double)*channel_num);

    for(int i = 0; i < channel_num; ++i){
        for(int j = begin_id; j < end_id; ++j){
            channel_mean[i] += abs(data[i][j]);
        }
        channel_mean[i] = channel_mean[i]/(double)(end_id - begin_id);
    }

    int max_value = -1;
    int max_index = -1;
    for(int i = 0; i < channel_num; ++i) {
        int value = 0;     
        for(int j = i; j < i+4; ++j){
            value += channel_mean[j%channel_num];
        }
        if(value > max_value){
            max_value = value;
            max_index = i;
        }
    }

    cout << "------calibrate channel mean value ----" << endl;
    for(int i = 0; i < channel_num; ++i)    cout << channel_mean[i] << ' ';
    cout << "begin index: "<< max_index << endl;
    delete []channel_mean;

    return max_index;
}

int calibrate_channel_6mic(double** data, int channel_num, int begin_id, int end_id){
    double *channel_mean = new double[channel_num]; 

    memset(channel_mean, 0,  sizeof(double)*channel_num);

    for(int i = 0; i < channel_num; ++i){
        for(int j = begin_id; j < end_id; ++j){
            channel_mean[i] += abs(data[i][j]);
        }
        channel_mean[i] = channel_mean[i]/(double)(end_id - begin_id);
    }

    int max_value = -1;
    int max_index = -1;
    for(int i = 0; i < channel_num; ++i) {
        int value = 0;     
        for(int j = i; j < i+6; ++j){
            value += channel_mean[j%channel_num];
        }
        if(value > max_value){
            max_value = value;
            max_index = i;
        }
    }

    cout << "------calibrate channel mean value ----" << endl;
    for(int i = 0; i < channel_num; ++i)    cout << channel_mean[i] << ' ';
    cout << "begin index: "<< max_index << endl;
    delete []channel_mean;

    return max_index;
}



double multiptle_sum(double* in1, int begin1, int end1, double * in2 ,int begin2, int end2){
    assert(end1-begin1 == end2 - begin2);
    double sum0 = 0;
    for(int i = 0; i < end1 - begin1 ; ++i){
        sum0 += in1[begin1 + i]*in2[begin2+i];
    }

    return sum0;
}

double sum_sqr(double* in, int begin, int end){
    double sum0 = 0;
    for(int i = begin; i < end; ++i){
        sum0 += in[i]*in[i];
    }
    return sum0;
}


int naiser_corr(double* signal, int total_length , int Nu, int N0, int DIVIDE_FACTOR, bool include_zero, double threshold){
    short PN_seq[8] = {1, -1, -1, -1, -1, -1, 1, -1};


    int N_both = Nu + N0;
    int preamble_L = 8*N_both;

    
    if(total_length-preamble_L < 0){
        cout << "input signal too short" << endl;
        return -1;
    }

    int len_corr = (total_length - preamble_L)/DIVIDE_FACTOR + 2;
    int num = 0;
    double* Mn = new double[len_corr];
    memset(Mn, 0, sizeof(double)*len_corr);

    for(int i = 0; i < total_length - preamble_L - 1; i += DIVIDE_FACTOR){
        //cout << i << ' ' << total_length - preamble_L - 1 << ' '  << len_corr  << endl;
        double Pd = 0;
        for(int k = 0; k < 8 - 1; ++k){
            int bk = PN_seq[k]*PN_seq[k+1];
            double temp_P = multiptle_sum(signal + i, k*N_both + N0, (k + 1)*N_both, signal + i, (k + 1)*N_both + N0, (k+2)*N_both);
            Pd +=  temp_P*bk;
        
        }
        double Rd = 0;
        if(include_zero){
            Rd = (sum_sqr(signal + i, 0, preamble_L)*Nu)/N_both;
        }
        else{
            for(int k = 0; k < 8; ++k){
                Rd += sum_sqr(signal + i, k*N_both + N0, (k + 1)*N_both);
            }
        }

        Mn[num] = Pd/Rd;

        num ++;
    }

    int max_idx = find_max(Mn, (unsigned long int)len_corr);
    double max_value = Mn[max_idx];

    if(max_value < threshold){
        delete [] Mn;
        return -1;
    }

    double shoulder = 0.85*max_value;
    int right = -1;
    int left = -1;


    // find the right shoulder point
    for(int i = max_idx; i < len_corr - 1; ++i){
        if(Mn[i] >= shoulder && Mn[i+1] <= shoulder){
            right = i;
            break;
        }
    }

    // find the left shoulder point
    for(int i = max_idx; i > 0; --i){
        if(Mn[i] >= shoulder && Mn[i-1] <= shoulder){
            left = i;
            break;
        }
    }

    if(left == -1 || right == -1){
        delete [] Mn;
        return max_idx*DIVIDE_FACTOR;
    }
    else{
        double middle = (right+left)*DIVIDE_FACTOR/2;
        delete [] Mn;
        return lround(middle);
    }
    
    

}




// z = abs(c)
double Complex_abs(fftw_complex c){
    return sqrt(pow(c[0], 2 ) + pow(c[1], 2 ) );
}

// z = power(c, 2)
double Complex_power(fftw_complex c){
    return pow(c[0], 2 ) + pow(c[1], 2 );
}

// z = x/y
void Complex_divide(fftw_complex x, fftw_complex y, fftw_complex z){ 
    double y_pow = Complex_power(y);
    z[0] = (x[0]*y[0] + x[1]*y[1])/y_pow;
    z[1] = (x[1]*y[0] - x[0]*y[1])/y_pow;
}

void abs_complex(double * out, fftw_complex* in, unsigned long int Nu){
    unsigned long int  i = 0;
    for( i = 0; i < Nu ; ++i){
        out[i] = Complex_abs(in[i]);
    }
}

void real_complex(double * out, fftw_complex* in, unsigned long int Nu){
    unsigned long int  i = 0;
    for( i = 0; i < Nu ; ++i){
        out[i] = in[i][0];
    }
}

void imag_complex(double * out, fftw_complex* in, unsigned long int Nu){
    unsigned long int  i = 0;
    for( i = 0; i < Nu ; ++i){
        out[i] = in[i][1];
    }
}


void estimate_H(fftw_complex* H, fftw_complex* X, fftw_complex * Y, int fft_begin, int fft_end, unsigned long int Nu){
    unsigned long int  i = 0;
    for( i = 0; i < Nu ; ++i){
        if(i >= fft_begin && i < fft_end){
            double y_pow = Complex_power(X[i]);
            H[i][0] = (Y[i][0]*X[i][0] + X[i][1]*Y[i][1])/y_pow;
            H[i][1] = (Y[i][1]*X[i][0] - Y[i][0]*X[i][1])/y_pow;
        }
        else{
            H[i][0] = 0;
            H[i][1] = 0;
        }
    }
}

void channel_estimation_freq(double* result, double * tx, double * rx, unsigned long int sig_len, int BW1, int BW2, int fs){
    unsigned long int i = 0;

    fftw_complex *rx_fft=(fftw_complex*)fftw_malloc(sig_len*sizeof(fftw_complex));
    if(rx_fft==NULL){fprintf(stderr,"malloc failed\n");exit(1);} 

    fftw_complex *tx_fft=(fftw_complex*)fftw_malloc(sig_len*sizeof(fftw_complex));
    if(tx_fft==NULL){fprintf(stderr,"malloc failed\n");exit(1);}
    
    
    fftw_plan p1=fftw_plan_dft_r2c_1d(sig_len, tx, tx_fft, FFTW_ESTIMATE);
    if (p1==NULL){fprintf(stderr,"plan creation failed\n");exit(1);}

    fftw_plan p2=fftw_plan_dft_r2c_1d(sig_len, rx, rx_fft, FFTW_ESTIMATE);
    if (p2==NULL){fprintf(stderr,"plan creation failed\n");exit(1);}


    fftw_execute(p1);
    fftw_execute(p2);
    //cout <<sig_len <<endl;

    fftw_complex *H=(fftw_complex*)fftw_malloc(sig_len*sizeof(fftw_complex));
    if(H==NULL){fprintf(stderr,"malloc failed\n");exit(1);}


    double delta_f = (double)fs/(double)sig_len;
    int fft_begin = (int) round(BW1/delta_f);
    int fft_end = (int) round(BW2/delta_f);
    /*
    for(i = fft_begin; i < fft_end; ++i){
        cout << tx_fft[i][0] << " + i" <<  tx_fft[i][1]<< endl; 
    }

    cout << "-------------" << endl;

    for(i = fft_begin; i < fft_end; ++i){
        cout << rx_fft[i][0] << " + i" <<  rx_fft[i][1]<< endl; 
    }
    cout << "-------------" << endl;
    */
    estimate_H(H, tx_fft, rx_fft, fft_begin, fft_end, sig_len);
    
    fftw_complex *h=(fftw_complex*)fftw_malloc(sig_len*sizeof(fftw_complex));
    if(h==NULL){fprintf(stderr,"malloc failed\n");exit(1);}

    fftw_plan ifft=fftw_plan_dft_1d(sig_len, H, h, FFTW_BACKWARD, FFTW_ESTIMATE);

    //double* abs_h = new double[Nu];

    fftw_execute(ifft);

    abs_complex(result, h, sig_len);

    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(ifft);

    fftw_free(tx_fft);
    fftw_free(rx_fft);
    fftw_free(H);
    fftw_free(h);

    //delete []abs_h;
}


void channel_estimation_freq_single(fftw_complex* H_result, double * tx, double * rx, unsigned long int sig_len, int BW1, int BW2, int fs){

    fftw_complex *rx_fft=(fftw_complex*)fftw_malloc(sig_len*sizeof(fftw_complex));
    if(rx_fft==NULL){fprintf(stderr,"malloc failed\n");exit(1);} 

    fftw_complex *tx_fft=(fftw_complex*)fftw_malloc(sig_len*sizeof(fftw_complex));
    if(tx_fft==NULL){fprintf(stderr,"malloc failed\n");exit(1);}
    
    
    fftw_plan p1=fftw_plan_dft_r2c_1d(sig_len, tx, tx_fft, FFTW_ESTIMATE);
    if (p1==NULL){fprintf(stderr,"plan creation failed\n");exit(1);}

    fftw_plan p2=fftw_plan_dft_r2c_1d(sig_len, rx, rx_fft, FFTW_ESTIMATE);
    if (p2==NULL){fprintf(stderr,"plan creation failed\n");exit(1);}


    fftw_execute(p1);
    fftw_execute(p2);
    //cout <<sig_len <<endl;

    double delta_f = (double)fs/(double)sig_len;
    int fft_begin = (int) round(BW1/delta_f);
    int fft_end = (int) round(BW2/delta_f);

    estimate_H(H_result, tx_fft, rx_fft, fft_begin, fft_end, sig_len);

    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);

    fftw_free(tx_fft);
    fftw_free(rx_fft);

    //delete []abs_h;
}

// result, tx -> channel length = Ns
// rx -> recv_len = 8*(Ns+N0) 
void channel_estimation_freq_multiple(double* result, double * tx, unsigned long int Ns, double * rx, unsigned long int recv_len, unsigned long int N0, int BW1, int BW2, int fs){
    assert(  recv_len == 8*(Ns+N0) );
    short PN_seq[8] = {1, -1, -1, -1, -1, -1, 1, -1};

    fftw_complex *H_avg =(fftw_complex*)fftw_malloc(Ns*sizeof(fftw_complex)); 
    if(H_avg==NULL){fprintf(stderr,"malloc failed\n");exit(1);}

    for(int i = 0; i<8; ++i){
        fftw_complex *H_result =(fftw_complex*)fftw_malloc(Ns*sizeof(fftw_complex)); 
        if(H_result==NULL){fprintf(stderr,"malloc failed\n");exit(1);}

        channel_estimation_freq_single(H_result, tx, rx + i*(N0+Ns) + N0, Ns, BW1, BW2, fs);
        if(i == 0){
            for(unsigned long int j = 0; j < Ns; ++j ){
                H_avg[j][0] = PN_seq[i]*H_result[j][0];
                H_avg[j][1] = PN_seq[i]*H_result[j][1];
            }
        }
        else{
            for(unsigned long int j = 0; j < Ns; ++j ){
                H_avg[j][0] += PN_seq[i]*H_result[j][0];
                H_avg[j][1] += PN_seq[i]*H_result[j][1];
            }            
        }
        
        fftw_free(H_result);
        
    }

    for(unsigned long int j = 0; j < Ns; ++j ){
        H_avg[j][0] /= (8*Ns); //H_result[j][0];
        H_avg[j][1] /= (8*Ns); //H_result[j][1];
    }         


    fftw_complex *h=(fftw_complex*)fftw_malloc(Ns*sizeof(fftw_complex));
    if(h==NULL){fprintf(stderr,"malloc failed\n");exit(1);}

    fftw_plan ifft=fftw_plan_dft_1d(Ns, H_avg, h, FFTW_BACKWARD, FFTW_ESTIMATE);

    //double* abs_h = new double[Nu];

    fftw_execute(ifft);

    abs_complex(result, h, Ns);

    fftw_destroy_plan(ifft);

    fftw_free(H_avg);
    fftw_free(h);
}   