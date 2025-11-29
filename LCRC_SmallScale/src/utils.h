#ifndef _UTILS_H
#define _UTILS_H
#include <cstdlib>
#include <cstdint>

#define NAME2STR(name) (#name)
#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

void reverseComplement(char *seq, int len);

void PRNG(unsigned *randVal, int order, unsigned num, size_t seed);
void feedbackCoefficients_PRNG(char *feedCoefficient, int order, char id);


int getMaxIndex(const int *val, unsigned len);
int getMinIndex(const int *val, unsigned len);

void createFolder(const char *folderName);
bool isPrime(unsigned val);
bool isCopirme(unsigned *len, int num);
int gcd(int a, int b);
int quickpow(int a, int b, int c);

void *alloc_init(size_t n, size_t n_size);

#endif
