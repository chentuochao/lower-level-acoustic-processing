cmake_minimum_required(VERSION 3.10)
project(beamform_mic_array)

set(CMAKE_CXX_STANDARD 17)
link_directories("/usr/local/lib")
link_directories("/usr/lib")
include_directories("/usr/local/include")
include_directories("/usr/include/python3.9")

add_executable(offline offline.cpp mic_array.cpp util.cpp)

SET(CMAKE_CXX_FLAGS "-DASIO_STANDALONE -Wall -O3")
target_link_libraries(offline pthread portaudio sndfile fftw3 python3.9)
