#ifndef UTIL_H
#define UTIL_H

#include <math.h>
#include <vector>
#include "readerwriterqueue.h"
#include <fftw3.h>
#include <iostream>
#include<algorithm>

//void convert_queue_to_fftcomplex();

void multiple(fftw_complex * cross_vector, fftw_complex * in1, fftw_complex * in2, unsigned long int len_in);

void xcorr(double* corr_result, double *sig, double *preamble, unsigned long int sig_len, unsigned long int preamble_len, unsigned long int max_len);

int find_max(double* in, unsigned long int len_in);

int find_min(double* in, unsigned long int len_in);

int calibrate_channel_4mic(double** data, int channel_num, int begin_id, int end_id);
int calibrate_channel_6mic(double** data, int channel_num, int begin_id, int end_id);

int  naiser_corr(double* signal, int total_length , int Nu, int N0, int DIVIDE_FACTOR, bool include_zero, double threshold);
double sum_sqr(double* in, int begin, int end);
double multiptle_sum(double* in1, int begin1, int end1, double * in2 ,int begin2, int end2);


double Complex_abs(fftw_complex c);
double Complex_power(fftw_complex c);
void Complex_divide(fftw_complex x, fftw_complex y, fftw_complex z);
void abs_complex(double * out, fftw_complex* in, unsigned long int Nu);
void real_complex(double * out, fftw_complex* in, unsigned long int Nu);
void imag_complex(double * out, fftw_complex* in, unsigned long int Nu);
void estimate_H(fftw_complex* H, fftw_complex* X, fftw_complex * Y, int fft_begin, int fft_end, unsigned long int Nu);
void channel_estimation_freq(double* result, double * tx, double * rx, unsigned long int sig_len, int BW1, int BW2, int fs);

void channel_estimation_freq_single(fftw_complex* H_result, double * tx, double * rx, unsigned long int sig_len, int BW1, int BW2, int fs);

void channel_estimation_freq_multiple(double* result, double * tx, unsigned long int Ns, double * rx, unsigned long int recv_len, unsigned long int N0, int BW1, int BW2, int fs);

#endif
