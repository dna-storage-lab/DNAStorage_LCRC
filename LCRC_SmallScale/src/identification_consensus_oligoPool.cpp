/************************************
 * Copyright (C) 2025 Tianjin University (TJU)
 *
 * @brief: Read identification and bit-wise majority voting
 *  for DNA data storage with small-scale oligo pools
 *
 *************************************/
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cstdint>
#include <cassert>
#include <cctype>
#include <cmath>

#include <pthread.h>

#include "compositeCode.h"
#include "utils.h"
#include "edlib.h"

#define MAX_COMPONENTS 5 // Number of component codes
#define THREAD_MAX 128   // Thread limits

#define OLIGO_NUM 11745
#define OLIGO_LEN 200
#define PAYLOAD_LEN 160
#define PRIMER_LEN 20

#define LDPC_N 64800
#define BLOCKCODE_NUM 29

#define MAX_READ_LEN 1000

#define N_TYPE 2
const char element[N_TYPE] = {'0', '1'};

enum Filetype
{
    TXT = 0,
    FASTQ = 1
};

typedef struct
{
    int n_component;
    uint64_t period;
    unsigned start_phase;

    uint64_t *M;
    uint64_t *M_sup;

    unsigned *cLen;
    int *cPhase;
    int **c;

} Configuration;

typedef struct
{
    int tid;

    unsigned limit_len_composite;

    unsigned window_len;

    bool mode;
    double threshold_primer;
    double threshold_autoCorr;
    double threshold_ED;

    size_t n_reads;
    size_t n_reads_primer_valid;
    size_t n_reads_valid;
    size_t n_reads_valid_corr;
    size_t tot_base;
    size_t reads_BWA;

    Filetype filetype;
    char reads_file[500];
    char *result_dir;

    char *shorten_compositeCode_refer;

    unsigned short *unique_pos[N_TYPE];

    Configuration *config;

} WorkerHandler;

const char *primerL = "ATAATTGGCTCCTGCTTGCA";
const char *primerR = "AATGTAGGCGGAAAGTGCAA";
const char *primerL_revComp = "TGCAAGCAGGAGCCAATTAT";
const char *primerR_revComp = "TTGCACTTTCCGCCTACATT";

const char *refer_file_PNseq = "PNseq_refer.fa";
const char *tmp_ofile_data_BWA = "tmp_read_BWA";
const char *tmp_ofile_PNseq_BWA = "tmp_PNseq_BWA";
const char *ofile_data_BWA = "data_BWA.txt";
const char *ofile_PNseq_BWA = "PNseq_BWA.fq";
const char *ofile_SAM_BWA = "alignment.sam";

void *pipelineWorker(void *args);

int main(int argc, char *argv[])
{
    setvbuf(stdout, NULL, _IONBF, 0); // none buffer stdout

    char result_dir[500] = {0};
    char reads_file_prefix[500] = {0};
    char refer_file[500] = {0};
    char consensus_file[500] = {0};

    int n_thread = 1;

    const int n_component = 5;
    int mode = 1;
    int bwa_idx_flag = 0;

    double threshold_primer = 0.2;
    double threshold_autoCorr = 0.8;

    unsigned cLen[MAX_COMPONENTS] = {31, 35, 43, 47, 59};
    int cPhase_init[MAX_COMPONENTS] = {0};

    int status;
    char cmdLine[1000] = {0};

    if (argc != 8)
    {
        fprintf(stderr, "Usage: "
                        "%s [Result dir] [Fastq file prefix] [Consensus file]"
                        "[Thread num] [Threshold] [Mode] [BWA index flag]\n",
                argv[0]);
        return 1;
    }
    sscanf(argv[1], "%s", result_dir);
    sscanf(argv[2], "%s", reads_file_prefix);
    sscanf(argv[3], "%s", consensus_file);
    sscanf(argv[4], "%d", &n_thread);
    sscanf(argv[5], "%lf", &threshold_autoCorr);
    sscanf(argv[6], "%d", &mode);
    sscanf(argv[7], "%d", &bwa_idx_flag);
    assert(n_component <= MAX_COMPONENTS);

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
    size_t limit_len_composite = data_size;

    uint64_t period = 1;
    for (int i = 0; i < n_component; i++)
    {
        period *= *(cLen + i);
        printf("C%d = %d, ", i + 1, *(cLen + i));
    }
    printf("\nFull-length composite code period = %lu\n", period);

    int **c = (int **)malloc(n_component * sizeof(int *));
    for (int i = 0; i < n_component; i++)
    {
        *(c + i) = (int *)alloc_init(*(cLen + i), sizeof(int));

        if (PNSequenceGenerator(*(c + i), *(cLen + i), 0))
            return 1;
    }

    uint64_t M[MAX_COMPONENTS] = {0}; // Parameters of Chinese Remainder Theorem
    uint64_t M_sup[MAX_COMPONENTS] = {0};
    getCRTparam(M, M_sup, cLen, n_component, period);

    uint64_t start_phase = performCRT(cPhase_init, M, M_sup, period, n_component);

    char *shorten_compositeCode_refer = (char *)alloc_init(limit_len_composite, sizeof(char));
    if (createCompositeCode_majority(shorten_compositeCode_refer, limit_len_composite, n_component, cPhase_init, cLen, c))
        return 1;
    for (size_t ti = 0; ti < limit_len_composite; ++ti)
        *(shorten_compositeCode_refer + ti) = *(shorten_compositeCode_refer + ti) == 1 ? '0' : '1';

#if 1 // Run it once for enhanced mode
    char *shorten_compositeCode_refer_AC = NULL;
    FILE *fpw_refPN = NULL;
    char refer_path_PNseq[500] = {0};
    if (mode)
    {
        sprintf(refer_path_PNseq, "%s/%s", result_dir, refer_file_PNseq);

        if (bwa_idx_flag)
        {
            shorten_compositeCode_refer_AC = (char *)alloc_init(limit_len_composite, sizeof(char));
            for (size_t ti = 0; ti < limit_len_composite; ++ti)
                *(shorten_compositeCode_refer_AC + ti) = *(shorten_compositeCode_refer + ti) == '0' ? 'A' : 'C';

            if ((fpw_refPN = fopen(refer_path_PNseq, "w")) == NULL)
            {
                fprintf(stderr, "Could not open file %s!\n", refer_path_PNseq);
                return 1;
            }
            fprintf(fpw_refPN, ">1\n");
            fprintf(fpw_refPN, "%s\n", shorten_compositeCode_refer_AC);
            fclose(fpw_refPN);

            sprintf(cmdLine, "bwa index %s", refer_path_PNseq);
            status = system(cmdLine);
            if (status == -1)
            {
                fprintf(stderr, "Errors occurred while BWA index!\n");
                return 1;
            }
            return 0;
        }
    }
#endif

    Configuration config;
    config.n_component = n_component;
    config.cLen = cLen;
    config.c = c;
    config.period = period;
    config.M = M;
    config.M_sup = M_sup;
    config.start_phase = start_phase;
    config.cPhase = cPhase_init;

    printf("\n----------Step 1: Read identification----------\n\n");

    printf("----------Component code correlation---------\n\n");

    //---- Multi-threaded processing
    pthread_t tid[THREAD_MAX];
    WorkerHandler worker_hdlr[THREAD_MAX];
    for (int ni = 0; ni < n_thread; ++ni)
    {
        (worker_hdlr + ni)->tid = ni;

        (worker_hdlr + ni)->limit_len_composite = limit_len_composite;

        (worker_hdlr + ni)->mode = mode;
        (worker_hdlr + ni)->threshold_primer = threshold_primer;
        (worker_hdlr + ni)->threshold_autoCorr = threshold_autoCorr;

        (worker_hdlr + ni)->filetype = Filetype::FASTQ;
        (worker_hdlr + ni)->result_dir = result_dir;
        sprintf((worker_hdlr + ni)->reads_file, "%s%03d", reads_file_prefix, ni);

        (worker_hdlr + ni)->shorten_compositeCode_refer = shorten_compositeCode_refer;

        for (int j = 0; j < N_TYPE; ++j)
            *((worker_hdlr + ni)->unique_pos + j) = (unsigned short *)alloc_init(data_size, sizeof(unsigned short));

        (worker_hdlr + ni)->config = &config;

        if (pthread_create(tid + ni, NULL, &pipelineWorker, (void *)(worker_hdlr + ni)) != 0)
        {
            fprintf(stderr, "Failed to create thread-%d!\n", ni);
            exit(1);
        }
    }

    for (int ni = 0; ni < n_thread; ++ni)
    {
        if (pthread_join(*(tid + ni), NULL) != 0)
        {
            fprintf(stderr, "Thread-%d failed!\n", ni);
        }
    }

    char *consensus_codewords = (char *)alloc_init(data_size, sizeof(char));
    unsigned short *unique_pos[N_TYPE];
    for (int j = 0; j < N_TYPE; ++j)
        *(unique_pos + j) = (unsigned short *)alloc_init(data_size, sizeof(unsigned short));

    double *lratio_codewords = (double *)alloc_init(data_size, sizeof(double));

    FILE *fpr_dat = NULL, *fpr_sam = NULL;
    char filepath_data[500] = {0};
    char filepath_PNseq[500] = {0};
    char filepath_SAM[500] = {0};
    if (mode)
    {
        sprintf(filepath_data, "%s/%s", result_dir, ofile_data_BWA);
        sprintf(filepath_PNseq, "%s/%s", result_dir, ofile_PNseq_BWA);
        sprintf(filepath_SAM, "%s/%s", result_dir, ofile_SAM_BWA);
    }

    size_t reads_BWA = 0;
    for (int ni = 0; ni < n_thread; ++ni)
    {
        for (int sb = 0; sb < data_size; ++sb)
        {
            for (int j = 0; j < N_TYPE; ++j)
            {
                *(*(unique_pos + j) + sb) += *(*((worker_hdlr + ni)->unique_pos + j) + sb);
            }
        }

        if (mode)
        {
            reads_BWA += (worker_hdlr + ni)->reads_BWA;

            if (ni == 0)
                sprintf(cmdLine, "cat %s/%s%03d > %s", result_dir, tmp_ofile_data_BWA, ni, filepath_data);
            else
                sprintf(cmdLine, "cat %s/%s%03d >> %s", result_dir, tmp_ofile_data_BWA, ni, filepath_data);
            status = system(cmdLine);
            if (status == -1)
            {
                fprintf(stderr, "Errors occurred while merge files %s*!\n", tmp_ofile_data_BWA);
                return 1;
            }

            if (ni == 0)
                sprintf(cmdLine, "cat %s/%s%03d > %s", result_dir, tmp_ofile_PNseq_BWA, ni, filepath_PNseq);
            else
                sprintf(cmdLine, "cat %s/%s%03d >> %s", result_dir, tmp_ofile_PNseq_BWA, ni, filepath_PNseq);
            status = system(cmdLine);
            if (status == -1)
            {
                fprintf(stderr, "Errors occurred while merge files %s*!\n", tmp_ofile_PNseq_BWA);
                return 1;
            }
        }
    }

    if (mode)
    {
        sprintf(cmdLine, "rm %s/%s*", result_dir, tmp_ofile_data_BWA);
        status = system(cmdLine);
        sprintf(cmdLine, "rm %s/%s*", result_dir, tmp_ofile_PNseq_BWA);
        status = system(cmdLine);
    }

    int flag;
    uint64_t startPos;
    char header[500] = {0}, pre_header[500] = {0};
    char cigar[1000] = {0};
    char tmp_symbol;
    int op_num, n_bp;

    unsigned rec_payload_len, rec_payload_polished_len;
    char payload_noisy[MAX_READ_LEN] = {0};
    char payload_polished[MAX_READ_LEN] = {0};

    size_t n_reads_valid_align = 0;

    // -----------Composite code alignment----------- //
    if (mode && reads_BWA > 0)
    {
        printf("----------Composite code alignment---------\n\n");

        printf("Remain reads for BWA align: %lu\n\n", reads_BWA);

        sprintf(cmdLine, "bwa mem -t %d %s %s -o %s", n_thread, refer_path_PNseq, filepath_PNseq, filepath_SAM);
        status = system(cmdLine);
        if (status == -1)
        {
            fprintf(stderr, "Errors occurred while BWA mem!\n");
            return 1;
        }

        if ((fpr_dat = fopen(filepath_data, "r")) == NULL)
        {
            fprintf(stderr, "Could not open file %s!\n", filepath_data);
            return 1;
        }

        if ((fpr_sam = fopen(filepath_SAM, "r")) == NULL)
        {
            fprintf(stderr, "Could not open file %s!\n", filepath_SAM);
            return 1;
        }
        while (fscanf(fpr_sam, "%s", header) != EOF)
        {
            if (strlen(header) > 0)
            {
                if (*(header + 0) == '@')
                    fscanf(fpr_sam, "%*[^\n]%*c");
                else
                {
                    fscanf(fpr_sam, "%d", &flag);
                    if (flag == 4)
                    {
                        fscanf(fpr_sam, "%*[^\n]%*c");
                        fscanf(fpr_dat, "%*[^\n]%*c");
                    }
                    else if (flag == 0 || flag == 16)
                    {
                        ++n_reads_valid_align;

                        fscanf(fpr_sam, "%*s %lu %*s %s %*[^\n]s", &startPos, cigar);
                        startPos -= 1;
                        fscanf(fpr_dat, "%s\n", payload_noisy);
                        rec_payload_len = strlen(payload_noisy);

                        n_bp = 0;
                        op_num = 0;
                        rec_payload_polished_len = 0;
                        for (int ni = 0; ni < strlen(cigar); ++ni)
                        {
                            tmp_symbol = *(cigar + ni);
                            if (isdigit(tmp_symbol))
                            {
                                op_num = op_num * 10 + (tmp_symbol - '0');
                            }
                            else
                            {
                                if (tmp_symbol == 'I' || tmp_symbol == 'S' || tmp_symbol == 'H')
                                {
                                    n_bp += op_num;
                                }
                                else if (tmp_symbol == 'D')
                                {
                                    for (int ti = 0; ti < op_num; ++ti)
                                    {
                                        *(payload_polished + rec_payload_polished_len) = 'E';
                                        rec_payload_polished_len++;
                                    }
                                }
                                else if (tmp_symbol == 'M')
                                {
                                    for (int ti = 0; ti < op_num; ++ti)
                                    {
                                        *(payload_polished + rec_payload_polished_len) = *(payload_noisy + n_bp);
                                        rec_payload_polished_len++;
                                        n_bp++;
                                    }
                                }
                                op_num = 0;
                            }
                        }

                        for (int ti = 0; ti < rec_payload_polished_len; ++ti)
                        {
                            if (startPos + ti >= limit_len_composite || startPos + ti < 0)
                                continue;

                            int value = (*(payload_polished + ti) == 'E') ? 2 : *(payload_polished + ti) - 48;
                            if (value != 2)
                                ++*(*(unique_pos + value) + startPos + ti);
                        }
                    }
                    else
                        fscanf(fpr_sam, "%*[^\n]%*c");
                }
            }
        }
        fclose(fpr_dat);
        fclose(fpr_sam);

        printf("\n");
    }

    size_t n_reads = 0;
    size_t n_reads_primer_valid = 0;
    size_t n_reads_valid = 0;
    size_t n_reads_valid_corr = 0;
    size_t tot_base = 0;
    for (int ni = 0; ni < n_thread; ++ni)
    {
        n_reads += (worker_hdlr + ni)->n_reads;
        n_reads_primer_valid += (worker_hdlr + ni)->n_reads_primer_valid;

        n_reads_valid += (worker_hdlr + ni)->n_reads_valid;
        n_reads_valid_corr += (worker_hdlr + ni)->n_reads_valid_corr;
        tot_base += (worker_hdlr + ni)->tot_base;
    }
    n_reads_valid += n_reads_valid_align;
    printf("Total reads = %lu\n", n_reads);
    printf("Coverage = %.2fx\n", 1.0 * tot_base / data_size);
    printf("Reads with valid primers = %lu\n", n_reads_primer_valid);
    printf("Reads with valid indices = %lu\n", n_reads_valid);
    printf("  1) correlation = %lu (%.2f)\n", n_reads_valid_corr, 1.0 * n_reads_valid_corr / n_reads_valid * 100);
    printf("  2) alignment   = %lu (%.2f)\n\n", n_reads_valid_align, 1.0 * n_reads_valid_align / n_reads_valid * 100);

    printf("----------Step 2: Bit-wise majority voting---------");
    //------------------------ Bit-wise majority voting -------------------------//
    int confidence;
    size_t n_era_consensus = 0;
    double era_rate_consensus;
    int tmp_val;
    int bin_fre[N_TYPE] = {0};
    int arr_sub[N_TYPE] = {0};

    for (size_t sb = 0; sb < data_size; ++sb)
    {
        if (*(*(unique_pos + 0) + sb) + *(*(unique_pos + 1) + sb) == 0)
        {
            *(consensus_codewords + sb) = 'E';
            confidence = 0;
        }
        else
        {
            if (*(*(unique_pos + 0) + sb) == *(*(unique_pos + 1) + sb))
            {
                *(consensus_codewords + sb) = 'E';
                confidence = 0;
                ++n_era_consensus;
            }
            else if (*(*(unique_pos + 0) + sb) > *(*(unique_pos + 1) + sb))
            {
                *(consensus_codewords + sb) = '0';
                confidence = (int)*(*(unique_pos + 0) + sb) - (int)*(*(unique_pos + 1) + sb);
            }
            else
            {
                *(consensus_codewords + sb) = '1';
                confidence = (int)(*(*(unique_pos + 0) + sb) - (int)*(*(unique_pos + 1) + sb));
            }
        }

        *(lratio_codewords + sb) = exp(-2.0 * confidence);
    }
    era_rate_consensus = 1.0 * n_era_consensus / data_size;

    FILE *fpw_cons = NULL;
    if ((fpw_cons = fopen(consensus_file, "wb")) == NULL)
    {
        fprintf(stderr, "Could not open file %s!\n", consensus_file);
        return 1;
    }
    fwrite(lratio_codewords, sizeof(double), data_size, fpw_cons);
    fclose(fpw_cons);

    ///-------------------- Evaluate performance -----------------///

    // Error evaluation. 0: off, 1: on
#if 0
    const char *dna_refer = "reference/DNA_oligoPool/oligoPool_11.7K.fa";
    char *dna_refer = (char *)alloc_init(data_size, sizeof(char));
    char seq_buff[MAX_READ_LEN] = {0};
    FILE *fpr_refer = NULL;
    if ((fpr_refer = fopen(refer_file, "r")) == NULL)
    {
        fprintf(stderr, "Could not open file %s!\n", refer_file);
        return 1;
    }
    size_t nn = 0;
    while (fscanf(fpr_refer, "%*[^\n]%*c") != EOF)
    {
        if (fscanf(fpr_refer, "%s\n", seq_buff) == EOF)
            break;
        int cur_len = strlen(seq_buff);
        assert(cur_len == PAYLOAD_LEN + 2 * PRIMER_LEN);

        memcpy(dna_refer + PAYLOAD_LEN * (nn++), seq_buff + PRIMER_LEN, PAYLOAD_LEN);
    }
    fclose(fpr_refer);

    size_t erasure_num = 0, error_num = 0;
    double erasure_rate = 0.0, error_rate = 0.0;

    char bit1, bit2;
    for (size_t nn = 0; nn < data_size; ++nn)
    {
        bit1 = *(consensus_codewords + nn);
        bit2 = (*(dna_refer + nn) == 'A' || *(dna_refer + nn) == 'T') ? '0' : '1';

        if (bit1 != bit2)
        {
            if (bit1 == 'E' || bit2 == 'E')
            {
                ++erasure_num;
            }
            else
            {
                ++error_num;
            }
        }
    }
    error_rate = 1.0 * error_num / data_size;
    erasure_rate = 1.0 * erasure_num / data_size;

    printf("\n\n--------------------Summary--------------------\n");
    printf("Bit error num = %lu, bit error rate = %.8f\n", error_num, error_rate);
    printf("Bit erasure num = %lu, bit erasure rate = %.8f\n", erasure_num, erasure_rate);

#endif

    time(&end_time);
    strftime(timebuff, 100, "%Y-%m-%d %H:%M:%S ", localtime(&end_time));
    printf("\nEnd running: %s\n", timebuff);
    printf("Elapsed time: %.2lf sec\n\n", difftime(end_time, start_time));
    return 0;
}

void *pipelineWorker(void *args)
{
    WorkerHandler *worker_hdlr = (WorkerHandler *)args;

    FILE *fpr_read = NULL;
    if ((fpr_read = fopen(worker_hdlr->reads_file, "r")) == NULL)
    {
        fprintf(stderr, "Could not open file %s!\n", worker_hdlr->reads_file);
        return NULL;
    }

    FILE *fpw_data = NULL;
    char ofile_data[500] = {0};
    FILE *fpw_PNseq = NULL;
    char ofile_PNseq[500] = {0};
    if (worker_hdlr->mode)
    {
        sprintf(ofile_data, "%s/%s%03d", worker_hdlr->result_dir, tmp_ofile_data_BWA, worker_hdlr->tid);
        if ((fpw_data = fopen(ofile_data, "w")) == NULL)
        {
            fprintf(stderr, "Could not open file %s!\n", ofile_data);
        }

        sprintf(ofile_PNseq, "%s/%s%03d", worker_hdlr->result_dir, tmp_ofile_PNseq_BWA, worker_hdlr->tid);
        if ((fpw_PNseq = fopen(ofile_PNseq, "w")) == NULL)
        {
            fprintf(stderr, "Could not open file %s!\n", ofile_PNseq);
        }
    }

    int n_component = worker_hdlr->config->n_component;
    unsigned limit_len_composite = worker_hdlr->limit_len_composite;
    uint64_t period = worker_hdlr->config->period;
    unsigned start_phase = worker_hdlr->config->start_phase;
    unsigned *cLen = worker_hdlr->config->cLen;
    int **c = worker_hdlr->config->c;
    uint64_t *M = worker_hdlr->config->M;
    uint64_t *M_sup = worker_hdlr->config->M_sup;

    double threshold_primer = worker_hdlr->threshold_primer;
    double threshold_autoCorr = worker_hdlr->threshold_autoCorr;

    char *shorten_compositeCode_refer = worker_hdlr->shorten_compositeCode_refer;

    int forward_flag = 1;
    int read_len;
    int payload_len, payload_polished_len;

    char *read_base = (char *)alloc_init(1000000, sizeof(char));
    char payload_base[MAX_READ_LEN] = {0};
    char payload_noisy[MAX_READ_LEN] = {0};
    char compositeCode_noisy[MAX_READ_LEN] = {0};
    char compositeCode_noisy_char[MAX_READ_LEN] = {0};

    char compositeCode_regen[MAX_READ_LEN] = {0};
    char payload_polished[MAX_READ_LEN] = {0};

    char compositeCode_AC[MAX_READ_LEN] = {0};
    char quality[MAX_READ_LEN] = {0};

    int primer_valid_flag = 0;
    int start_loc, end_loc;
    int minED_l_op, minED_r_op;
    int minED_l[2], minED_r[2];
    EdlibEqualityPair additionalEqualities[4] = {{'N', 'A'}, {'N', 'T'}, {'N', 'G'}, {'N', 'C'}};

    int **crossCorrVal = (int **)malloc(n_component * sizeof(int *));
    for (int i = 0; i < n_component; ++i)
        crossCorrVal[i] = (int *)alloc_init(*(cLen + i), sizeof(int));

    int recCorrVal, idealCorrVal;
    uint64_t recLoc;
    int cPhase_rec[MAX_COMPONENTS] = {0};

    int tid_flag;
    int minED_align;
    uint64_t recLoc_align;

    int valid_flag = 0;
    int start_offset;
    int64_t recIndex;

    size_t n_reads = 0;
    size_t n_reads_primer_valid = 0;
    size_t n_reads_valid = 0;
    size_t n_reads_valid_corr = 0;
    size_t n_reads_valid_align = 0;
    size_t tot_base = 0;
    size_t reads_BWA = 0;
    while (fgets(read_base, MAX_READ_LEN, fpr_read))
    {
        if (worker_hdlr->filetype == Filetype::FASTQ)
        {
            fscanf(fpr_read, "%s\n", read_base);
            fscanf(fpr_read, "%*[^\n]%*c");
            fscanf(fpr_read, "%*[^\n]%*c");
            read_len = strlen(read_base);
        }
        else if (worker_hdlr->filetype == Filetype::TXT)
            read_len = strlen(read_base) - 1;

        if (read_len == 0 || read_len > MAX_READ_LEN)
            continue;
        ++n_reads;

        // ----------Paired-end perime alignment---------- //
        forward_flag = 1;
        primer_valid_flag = 0;
        EdlibAlignResult align_rslt_left[2];
        *(align_rslt_left + 0) = edlibAlign(primerL, PRIMER_LEN, read_base, read_len,
                                            edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 4));
        *(align_rslt_left + 1) = edlibAlign(primerR_revComp, PRIMER_LEN, read_base, read_len,
                                            edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 4));
        *(minED_l + 0) = (align_rslt_left + 0)->editDistance;
        *(minED_l + 1) = (align_rslt_left + 1)->editDistance;
        if (*(minED_l + 0) <= *(minED_l + 1))
        {
            forward_flag = 1;
            minED_l_op = *(minED_l + 0);
            start_loc = *((align_rslt_left + 0)->endLocations) + 1;
        }
        else
        {
            forward_flag = 0;
            minED_l_op = *(minED_l + 1);
            start_loc = *((align_rslt_left + 1)->endLocations) + 1;
        }
        if (minED_l_op <= threshold_primer * PRIMER_LEN)
        {
            primer_valid_flag = 1;

            EdlibAlignResult align_rslt_right1;
            if (forward_flag)
                align_rslt_right1 = edlibAlign(primerR, PRIMER_LEN, read_base, read_len,
                                               edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 4));
            else
                align_rslt_right1 = edlibAlign(primerL_revComp, PRIMER_LEN, read_base, read_len,
                                               edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 4));
            minED_r_op = align_rslt_right1.editDistance;
            if (minED_r_op <= threshold_primer * PRIMER_LEN)
            {
                end_loc = *align_rslt_right1.startLocations - 1;
            }
            else
            {
                end_loc = read_len - 1;
            }
            edlibFreeAlignResult(align_rslt_right1);
        }
        else
        {
            start_loc = 0;

            EdlibAlignResult align_rslt_right[2];
            *(align_rslt_right + 0) = edlibAlign(primerR, PRIMER_LEN, read_base, read_len,
                                                 edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 4));

            *(align_rslt_right + 1) = edlibAlign(primerL_revComp, PRIMER_LEN, read_base, read_len,
                                                 edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 4));
            *(minED_r + 0) = (align_rslt_right + 0)->editDistance;
            *(minED_r + 1) = (align_rslt_right + 1)->editDistance;
            if (*(minED_r + 0) <= *(minED_r + 1))
            {
                forward_flag = 1;
                minED_r_op = *(minED_r + 0);
                end_loc = *((align_rslt_right + 0)->startLocations) - 1;
            }
            else
            {
                forward_flag = 0;
                minED_r_op = *(minED_r + 1);
                end_loc = *((align_rslt_right + 1)->startLocations) - 1;
            }
            if (minED_r_op <= threshold_primer * PRIMER_LEN)
            {
                primer_valid_flag = 1;
            }
            for (int i = 0; i < 2; ++i)
                edlibFreeAlignResult(align_rslt_right[i]);
        }
        for (int i = 0; i < 2; ++i)
            edlibFreeAlignResult(align_rslt_left[i]);

        if (primer_valid_flag)
        {
            payload_len = end_loc - start_loc;
            if (payload_len <= 0)
                continue;
            tot_base += payload_len;

            ++n_reads_primer_valid;

            memset(payload_base, 0, MAX_READ_LEN);
            memset(payload_noisy, 0, MAX_READ_LEN);
            memset(compositeCode_noisy, 0, MAX_READ_LEN);
            memset(compositeCode_noisy_char, 0, MAX_READ_LEN);
            memcpy(payload_base, read_base + start_loc, payload_len);
            if (!forward_flag)
                reverseComplement(payload_base, payload_len);

            memset(compositeCode_AC, 0, MAX_READ_LEN);
            memset(quality, 0, MAX_READ_LEN);
            for (int ti = 0; ti < payload_len; ++ti)
            {
                if (*(payload_base + ti) == 'A')
                {
                    *(payload_noisy + ti) = '0';
                    *(compositeCode_noisy + ti) = 0;
                }
                else if (*(payload_base + ti) == 'T')
                {
                    *(payload_noisy + ti) = '0';
                    *(compositeCode_noisy + ti) = 1;
                }
                else if (*(payload_base + ti) == 'G')
                {
                    *(payload_noisy + ti) = '1';
                    *(compositeCode_noisy + ti) = 0;
                }
                else if (*(payload_base + ti) == 'C')
                {
                    *(payload_noisy + ti) = '1';
                    *(compositeCode_noisy + ti) = 1;
                }
                else
                {
                    *(payload_noisy + ti) = 'E';
                    *(compositeCode_noisy + ti) = 0;
                }
                *(compositeCode_noisy + ti) = 1 - 2 * *(compositeCode_noisy + ti);
                *(compositeCode_noisy_char + ti) = *(compositeCode_noisy + ti) == 1 ? '0' : '1';

                *(compositeCode_AC + ti) = *(compositeCode_noisy + ti) == 1 ? 'A' : 'C';
                *(quality + ti) = 'F';
            }

            // -----------component code correlation----------- //

            if (threshold_autoCorr <= 1)
            {
                for (int nc = 0; nc < n_component; ++nc)
                {
                    probeCorrelation(*(crossCorrVal + nc), compositeCode_noisy, payload_len, *(c + nc), *(cLen + nc));

                    *(cPhase_rec + nc) = getMaxIndex(*(crossCorrVal + nc), *(cLen + nc));
                }

                recLoc = performCRT(cPhase_rec, M, M_sup, period, n_component);

                createCompositeCode_majority(compositeCode_regen, payload_len, n_component, cPhase_rec, cLen, c);
                idealCorrVal = payload_len;
                recCorrVal = autoCorrelation(compositeCode_noisy, payload_len, compositeCode_regen, payload_len, 0);

                recIndex = (int64_t)recLoc - (int64_t)start_phase;
            }
            else
            {
                recCorrVal = 0;
                idealCorrVal = 1;
            }

            valid_flag = 0;
            start_offset = 0;
            payload_polished_len = payload_len;
            memcpy(payload_polished, payload_noisy, payload_len);
            *(payload_polished + payload_len) = '\0';
            if (recCorrVal >= threshold_autoCorr * idealCorrVal)
            {
                valid_flag = 1;

                for (int ni = 0; ni < payload_len; ++ni)
                {
                    if (*(compositeCode_noisy + ni) != *(compositeCode_regen + ni))
                    {
                        *(payload_polished + ni) = 'E';
                    }
                }
            }
            else
            {
                if (worker_hdlr->mode)
                {
                    fprintf(fpw_data, "%s\n", payload_noisy);

                    fprintf(fpw_PNseq, "@%d%lu\n", worker_hdlr->tid, reads_BWA);
                    fprintf(fpw_PNseq, "%s\n", compositeCode_AC);
                    fprintf(fpw_PNseq, "+\n");
                    fprintf(fpw_PNseq, "%s\n", quality);
                    ++reads_BWA;
                }
            }

            if (valid_flag)
            {
                if (valid_flag == 1)
                    ++n_reads_valid_corr;
                ++n_reads_valid;

                for (int ti = 0; ti < payload_polished_len; ++ti)
                {
                    if (recIndex + start_offset + ti >= limit_len_composite || recIndex + start_offset + ti < 0)
                        continue;

                    int value = (*(payload_polished + ti) == 'E') ? 2 : *(payload_polished + ti) - 48;
                    if (value != 2)
                        ++*(*(worker_hdlr->unique_pos + value) + recIndex + start_offset + ti);
                }
            }
        }
    }
    fclose(fpr_read);

    if (worker_hdlr->mode)
    {
        fclose(fpw_data);
        fclose(fpw_PNseq);
    }

    worker_hdlr->n_reads = n_reads;
    worker_hdlr->n_reads_primer_valid = n_reads_primer_valid;
    worker_hdlr->n_reads_valid = n_reads_valid;
    worker_hdlr->n_reads_valid_corr = n_reads_valid_corr;
    worker_hdlr->tot_base = tot_base;
    worker_hdlr->reads_BWA = reads_BWA;

    return NULL;
}
