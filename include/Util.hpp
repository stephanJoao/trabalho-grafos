// For compiling with Mersenne twister random number generator include MTWISTER
// compiling directive

#ifndef UTIL_H
#define UTIL_H


#define MTWISTER

#include <ctime>
#include <unistd.h>
#include <sys/times.h>

#include <vector>

using namespace std;
//#ifdef MTWISTER
//#define MTWISTER

void init_genrand(unsigned long s);
unsigned long genrand_int32(void);
unsigned int intRandom(const unsigned int maxValue);

#endif
