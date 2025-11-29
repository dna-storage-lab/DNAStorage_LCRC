#ifndef _COMPOSITECODE_H
#define _COMPOSITECODE_H
#include <cstdint>
#include <cstdlib>

int createCompositeCode_majority(char *compositeCode,
                                 size_t codeLen,
                                 unsigned cNum,
                                 const int *cPhase,
                                 const unsigned *cLen,
                                 int **c);
int PNSequenceGenerator(int *sequence, unsigned period, char P_id = 0);
void LegendreCode(int *sequence, unsigned period);
void twoPrimeCode(int *sequence, unsigned period);
void mSequenceGenerator(int *mSeq, int mOrder, char P_id = 0);
void feedbackCoefficients(char *feedCoefficient, int mOrder, char P_id = 0);

void probeCorrelation(int *corrVal, const char *compositeSeq, unsigned seqLen, const int *component, unsigned cLen);
int autoCorrelation(const char *rvdSeq, unsigned rvdSeqLen, const char *refSeq, size_t refSeqLen, size_t offset = 0);

void getCRTparam(size_t *M, size_t *M_sup, const unsigned *cLen, unsigned num, size_t period);
size_t performCRT(int *phase, const size_t *M, const size_t *M_sup, size_t period, unsigned num);

#endif
