#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>

#include "compositeCode.h"
#include "utils.h"

#define OLIGO_NUM 299700
#define PAYLOAD_LEN 160

const int n_component = 5;
const char *refer_file_PNseq = "config/PNseq_refer.fa";

int main()
{

    unsigned cLen[n_component] = {31, 35, 43, 47, 59};
    int cPhase_init[n_component] = {0};

    if (!isCopirme(cLen, n_component))
    {
        fprintf(stderr, "\nPeriods of components must be coprime!\n\n");
        return 1;
    }

    time_t start_time, end_time;
    char timebuff[100];
    time(&start_time);
    strftime(timebuff, 100, "%Y-%m-%d %H:%M:%S ", localtime(&start_time));
    printf("\nStart running: %s\n\n", timebuff);

    const size_t data_size = OLIGO_NUM * PAYLOAD_LEN;

    uint64_t period = 1;
    for (int i = 0; i < n_component; i++)
    {
        period *= *(cLen + i);
        // printf("C%d = %d, ", i + 1, *(cLen + i));
    }
    // printf("Period = %lu\n", period);

    size_t limit_len_composite = data_size;

    int **c = (int **)malloc(n_component * sizeof(int *)); // Component sequences
    for (int i = 0; i < n_component; i++)
    {
        *(c + i) = (int *)calloc(*(cLen + i), sizeof(int));
        int ret = PNSequenceGenerator(*(c + i), *(cLen + i), 0);
        if (ret)
            return 1;
    }

    uint64_t M[n_component] = {0}; // Parameters of Chinese Remainder Theorem
    uint64_t M_sup[n_component] = {0};
    getCRTparam(M, M_sup, cLen, n_component, period);

    uint64_t start_phase = performCRT(cPhase_init, M, M_sup, period, n_component);

    char *compositeCode_refer_AC = (char *)alloc_init(limit_len_composite, sizeof(char));

    printf("Generating LCRC...\n");
    if (createCompositeCode_majority(compositeCode_refer_AC, limit_len_composite, n_component, cPhase_init, cLen, c))
        return 1;
    for (size_t ti = 0; ti < limit_len_composite; ++ti)
        *(compositeCode_refer_AC + ti) = *(compositeCode_refer_AC + ti) == 1 ? 'A' : 'C';
    FILE *fpw_refPN = NULL;
    if ((fpw_refPN = fopen(refer_file_PNseq, "w")) == NULL)
    {
        fprintf(stderr, "Could not open file %s!\n", refer_file_PNseq);
        return 1;
    }
    fprintf(fpw_refPN, ">1\n");
    fprintf(fpw_refPN, "%s\n", compositeCode_refer_AC);
    fclose(fpw_refPN);

    int status;
    char cmdLine[1000] = {0};

    printf("\nBuild bwa index...\n\n");

    sprintf(cmdLine, "bwa index %s", refer_file_PNseq);
    status = system(cmdLine);
    if (status == -1)
    {
        fprintf(stderr, "Errors occurred while BWA index!\n");
        return 1;
    }

    printf("\nThe generated LCRC is written to %s\n", refer_file_PNseq);

    time(&end_time);
    strftime(timebuff, 100, "%Y-%m-%d %H:%M:%S ", localtime(&end_time));
    printf("\nEnd running: %s\n", timebuff);
    printf("\nElapsed time: %.2lf sec\n\n", difftime(end_time, start_time));
    return 0;
}