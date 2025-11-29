/************************************
 * Copyright (C) 2025 Tianjin University (TJU)
 *
 * @brief: Real-time readout of small-scale oligo pools via nanopore sequencing
 *
 *************************************/
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cassert>

#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>

#include "edlib.h"

#include "compositeCode.h"
#include "utils.h"

#include "ldpc_decode.h"

#include <pthread.h>
#include "thpool.h"

#define DATA_RECOVERY
#ifdef DATA_RECOVERY
#define IMG1_BYTE 100161
#define IMG2_BYTE 95558
#endif

#define OLIGO_NUM 11745
#define OLIGO_LEN 200
#define PAYLOAD_LEN 160
#define PRIMER_LEN 20

#define MSEQ1 65535
#define MSEQ2 31

#define BLOCKCODE_NUM 29

#define LDPC_N 64800
#define LDPC_M 10800
#define LDPC_K 54000

#define MAX_COMPONENTS 5 // Number of component sequences
#define THREAD_MAX 60    // Thread limits

#define MAX_READ_LEN 1000
#define MAX_PAYLOAD_LEN 1000

#define N_TYPE 2
const char element[N_TYPE] = {'0', '1'};

const char *ldpc_pchk_file = "config/ldpc.pchk";

typedef struct
{
    int tid;
    int n_jobs;

    int process_code_num;

    int n_corrupted_code;

    char *mSeq1;
    char *mSeq2;

    unsigned *permutation;
    double *lratio_codewords;
    char *codeword_status;

    char **rec_payload;

    pthread_mutex_t mutex;

} DecoderContext;

typedef struct
{
    int tid;
    int n_jobs;

    unsigned data_size;
    int process_bits;

    unsigned avail_base_num;
    unsigned n_era_consensus;

    char *consensus_codewords; // '0'/'1'/'E'
    double *lratio_codewords;  // soft information

    unsigned short **unique_pos;

    pthread_mutex_t mutex;

} ConsensuContext;

const char *primerL = "ATAATTGGCTCCTGCTTGCA";
const char *primerR = "AATGTAGGCGGAAAGTGCAA";
const char *primerL_revComp = "TGCAAGCAGGAGCCAATTAT";
const char *primerR_revComp = "TTGCACTTTCCGCCTACATT";

const char *refer_file_PNseq = "config/PNseq_refer.fa";
const char *tmp_ofile_data_BWA = "tmp_read_BWA";
const char *tmp_ofile_PNseq_BWA = "tmp_PNseq_BWA";
const char *ofile_data_BWA = "data_BWA.txt";
const char *ofile_PNseq_BWA = "PNseq_BWA.fq";
const char *ofile_SAM_BWA = "alignment.sam";

// Reference files for recovery evaluation
const char *refer_file = "reference/DNA_oligoPool/oligoPool_11.7K.fa"; // reference FASTA
const char *srcdata_file = "reference/DNA_oligoPool/srcdata.txt";

char bit2ascii(const char *inBit)
{
    char symbol = 0;
    for (int i = 0; i < 8; i++)
    {
        symbol += *(inBit + i) * (1 << (7 - i));
    }
    return symbol;
}

void ConsensusWorker(void *args);
void LDPCDecodingWorker(void *args);

size_t current_time_ms()
{
    struct timeval time_now;
    gettimeofday(&time_now, NULL);
    return (size_t)time_now.tv_sec * 1000 + time_now.tv_usec / 1000;
}

int isFileReady(const char *filename, size_t set_lines);

// const double available_ratio = 0.938; // offline
// const double available_ratio = 0.955;    // online
// Max iteration for ldpc decoding
const int maxIteration_ldpc = 40;

int main(int argc, char *argv[])
{
    setvbuf(stdout, NULL, _IONBF, 0); // none buffer stdout

    char result_dir[500] = {0};
    char reads_file_prefix[500] = {0};
    char consensus_file[500] = {0};

    int status;
    char cmdLine[1000] = {0};

    const int n_component = 5;

    int target_reads = -1;

    int n_thread = 1;
    int mode = 1; // 0: correlation-based, 1: alignment-based

    double threshold_primer = 0.2;
    double threshold_autoCorr = 0.8;

    double available_ratio = 0.938;

    unsigned cLen[MAX_COMPONENTS] = {31, 35, 43, 47, 59};
    int cPhase_init[MAX_COMPONENTS] = {0};

    if (argc < 7)
    {
        fprintf(stderr, "Usage: "
                        "%s [result dir] [Fastq file] [Threads] [Correlation threshold] [Mode] [Base available threshold]"
                        "[Target reads per file]\n",
                argv[0]);
        return 1;
    }
    sscanf(argv[1], "%s", result_dir);
    sscanf(argv[2], "%s", reads_file_prefix);
    sscanf(argv[3], "%d", &n_thread);
    sscanf(argv[4], "%lf", &threshold_autoCorr);
    sscanf(argv[5], "%d", &mode);
    sscanf(argv[6], "%lf", &available_ratio);
    assert(n_component <= MAX_COMPONENTS);

    if (argc == 8)
        sscanf(argv[7], "%d", &target_reads);

    time_t start_time, end_time;
    char timebuff[100];
    time(&start_time);
    strftime(timebuff, 100, "%Y-%m-%d %H:%M:%S ", localtime(&start_time));
    printf("\nStart running: %s\n\n", timebuff);

    const unsigned data_size = OLIGO_NUM * PAYLOAD_LEN;
    const unsigned limit_len = data_size;

    int read_len;
    int tot_block = 0;
    int n_read = 0;
    int n_reads_BWA = 0;

    size_t tot_reads = 0;
    size_t tot_read_BWA = 0;

    int forward_flag = 1;
    int primer_valid_flag = 0;
    int start_loc, end_loc;
    int minED_l_op, minED_r_op;
    int minED_l[2], minED_r[2];
    EdlibEqualityPair additionalEqualities[4] = {{'N', 'A'}, {'N', 'T'}, {'N', 'G'}, {'N', 'C'}};

    FILE *fpr_read = NULL;
    char reads_file[500] = {0};

    FILE *fpw_PN = NULL, *fpw_dat = NULL;
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

    int flag;
    int64_t startPos;
    char tmp_symbol;
    int op_num, n_bp;
    int payload_len, payload_polished_len;
    char header[500] = {0};
    char cigar[1000] = {0};

    int confidence;
    int tmp_val;
    int bin_fre[N_TYPE] = {0};
    int arr_sub[N_TYPE] = {0};

    double read_coverage = 0.0;
    double payload_coverage = 0.0;

    unsigned avail_base_num = 0;
    double avail_base_ratio = 0.0;

    unsigned n_era_consensus = 0;
    double era_rate_consensus;

    size_t tot_read_base = 0;
    size_t tot_payload_base = 0;

    size_t tot_reads_primer_valid = 0;
    size_t tot_reads_valid = 0;
    size_t tot_reads_valid_corr = 0;
    size_t tot_reads_valid_align = 0;

    char payload_noisy_base[MAX_PAYLOAD_LEN] = {0};  // A, T, G, C
    char payload_noisy[MAX_PAYLOAD_LEN] = {0};       // '0', '1'
    char compositeCode_noisy[MAX_PAYLOAD_LEN] = {0}; // +1, -1

    char compositeCode_regen[MAX_PAYLOAD_LEN] = {0}; // +1, -1
    char payload_polished[MAX_PAYLOAD_LEN] = {0};    // '0', '1'

    char compositeCode_AC[MAX_PAYLOAD_LEN] = {0}; // '0' -> 'A', '1' -> 'C'
    char quality[MAX_PAYLOAD_LEN] = {0};
    memset(quality, 'F', MAX_PAYLOAD_LEN);

    char *read_base = (char *)alloc_init(1000000, sizeof(char)); // A, T, G, C

    int **crossCorrVal = (int **)malloc(n_component * sizeof(int *));
    for (int i = 0; i < n_component; ++i)
        *(crossCorrVal + i) = (int *)alloc_init(*(cLen + i), sizeof(int));

    int recCorrVal, idealCorrVal;
    uint64_t recLoc;
    int cPhase_rec[MAX_COMPONENTS] = {0}; // 捕获恢复的子码相位偏移量

    int valid_flag = 0; // 0: unvalid, 1: correlation-based, 2: alignment-based
    int64_t recIndex;

    unsigned short *unique_pos[N_TYPE];
    for (int j = 0; j < N_TYPE; ++j)
        *(unique_pos + j) = (unsigned short *)alloc_init(data_size, sizeof(unsigned short));

    char *consensus_codewords = (char *)alloc_init(data_size, sizeof(char));    // '0'/'1'/'E'
    double *lratio_codewords = (double *)alloc_init(data_size, sizeof(double)); // soft information

    char **rec_payload = (char **)malloc(BLOCKCODE_NUM * sizeof(char *));
    for (int i = 0; i < BLOCKCODE_NUM; ++i)
        *(rec_payload + i) = (char *)alloc_init(LDPC_K, sizeof(char));

    char codeword_status[BLOCKCODE_NUM] = {0}; // 0: success, 1:fail

    //----- Load m sequence (scrambler)
    unsigned mSeq1Len = MSEQ1;
    unsigned mSeq2Len = MSEQ2;
    char mSeq1[MSEQ1];
    char mSeq2[MSEQ2];
    FILE *fpr = NULL;
    if ((fpr = fopen("config/m65535.txt", "rb")) == NULL)
    {
        fprintf(stderr, "Could not open file config/m65535.txt!\n");
        return 1;
    }
    fread(mSeq1, sizeof(char), mSeq1Len, fpr);
    fclose(fpr);
    for (int i = 0; i < mSeq1Len; i++)
    {
        *(mSeq1 + i) -= '0';
    }

    if ((fpr = fopen("config/m31.txt", "rb")) == NULL)
    {
        fprintf(stderr, "Could not open file config/m31.txt!\n");
        return 1;
    }
    fread(mSeq2, sizeof(char), mSeq2Len, fpr);
    fclose(fpr);
    for (int i = 0; i < mSeq2Len; i++)
    {
        *(mSeq2 + i) -= '0';
    }

    //----- Generate permutation table
    unsigned *randVal = (unsigned *)alloc_init(65535, sizeof(unsigned));
    unsigned *permutation = (unsigned *)alloc_init(LDPC_N, sizeof(unsigned));
    PRNG(randVal, 16, 65535, 64); // seed = 64
    for (unsigned i = 0, nn = 0; i < 65535; ++i)
    {
        if (*(randVal + i) <= LDPC_N)
        {
            *(permutation + nn) = *(randVal + i) - 1;
            ++nn;
        }
    }

    threadpool thpool = thpool_init(n_thread + 2);

    uint64_t period = 1;
    for (int i = 0; i < n_component; i++)
    {
        period *= *(cLen + i);
    }

    int **c = (int **)malloc(n_component * sizeof(int *)); // Component sequences
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

    DecoderContext decoder_ctx;
    decoder_ctx.mSeq1 = mSeq1;
    decoder_ctx.mSeq2 = mSeq2;
    decoder_ctx.lratio_codewords = lratio_codewords;
    decoder_ctx.permutation = permutation;
    decoder_ctx.codeword_status = codeword_status;
    decoder_ctx.rec_payload = rec_payload;
    decoder_ctx.process_code_num = ceil(1.0 * BLOCKCODE_NUM / n_thread);
    int n_thread_on = n_thread * decoder_ctx.process_code_num > BLOCKCODE_NUM
                          ? ceil(1.0 * BLOCKCODE_NUM / decoder_ctx.process_code_num)
                          : n_thread;
    decoder_ctx.n_jobs = n_thread_on;
    pthread_mutex_init(&decoder_ctx.mutex, NULL);

    ConsensuContext cons_ctx;
    cons_ctx.data_size = data_size;
    cons_ctx.consensus_codewords = consensus_codewords;
    cons_ctx.lratio_codewords = lratio_codewords;
    cons_ctx.unique_pos = unique_pos;
    cons_ctx.process_bits = decoder_ctx.process_code_num * LDPC_N;
    cons_ctx.n_jobs = n_thread_on;
    pthread_mutex_init(&cons_ctx.mutex, NULL);

#ifdef DATA_RECOVERY
    //---- Data recovery ----//
    size_t tot_bytes = BLOCKCODE_NUM * LDPC_K / 8;
    char *pbuff_ptr = NULL;
    char *rec_data_bytes = (char *)alloc_init(tot_bytes, sizeof(char));
    int n_byte = 0;

    FILE *fpw = NULL;
    char img1_filename[500] = {0};
    char img2_filename[500] = {0};

    sprintf(img1_filename, "%s/img1.jpg", result_dir);
    sprintf(img2_filename, "%s/img2.jpg", result_dir);

    int open_flag = 0;
#endif

    time_t pre_primer_time, post_primer_time;
    time_t pre_corr_time, post_corr_time;
    time_t pre_align_time, post_align_time;
    double primer_used_time = 0.0;
    double corr_used_time = 0.0;
    double align_used_time = 0.0;

    printf("\n----------------Start read-by-read processing------------------\n");
    time_t pre_time, post_time;
    pre_time = current_time_ms();

    while (1)
    {
        sprintf(reads_file, "%s%d.fastq", reads_file_prefix, tot_block);

        if (target_reads > 0)
            while (!isFileReady(reads_file, target_reads))
                ;

        if ((fpr_read = fopen(reads_file, "r")) == NULL)
        {
            fprintf(stderr, "Could not open file %s!\n", reads_file);
            return 1;
        }

        if (mode)
        {
            n_reads_BWA = 0;
            if ((fpw_dat = fopen(filepath_data, "w")) == NULL) // writer for data
            {
                fprintf(stderr, "Could not open file %s!\n", filepath_data);
                return 1;
            }
            if ((fpw_PN = fopen(filepath_PNseq, "w")) == NULL) // writer for composite code
            {
                fprintf(stderr, "Could not open file %s!\n", filepath_PNseq);
                return 1;
            }
        }

        printf("\n--------------Step1: Primer & index identification-------------\n");

        n_read = 0;
        while (fscanf(fpr_read, "%*[^\n]%*c") != EOF)
        {
            fscanf(fpr_read, "%s\n", read_base); // reads
            fscanf(fpr_read, "%*[^\n]%*c");
            fscanf(fpr_read, "%*[^\n]%*c");
            read_len = strlen(read_base);

            if (read_len > MAX_READ_LEN)
                read_len = MAX_READ_LEN;

            tot_read_base += read_len;
            ++n_read;

            pre_primer_time = current_time_ms();

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
                // Forward
                forward_flag = 1;
                minED_l_op = *(minED_l + 0);
                start_loc = *((align_rslt_left + 0)->endLocations) + 1;
            }
            else
            {
                // Reverse
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
                    end_loc = start_loc + 160 - 1;
                    if (read_len < 160)
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
                    // Forward
                    forward_flag = 1;
                    minED_r_op = *(minED_r + 0);
                    end_loc = *((align_rslt_right + 0)->startLocations) - 1;
                }
                else
                {
                    // Reverse
                    forward_flag = 0;
                    minED_r_op = *(minED_r + 1);
                    end_loc = *((align_rslt_right + 1)->startLocations) - 1;
                }
                if (minED_r_op <= threshold_primer * PRIMER_LEN)
                {
                    if (end_loc >= 160)
                        start_loc = end_loc - 160 + 1;

                    primer_valid_flag = 1;
                }
                for (int i = 0; i < 2; ++i)
                    edlibFreeAlignResult(align_rslt_right[i]);
            }
            for (int i = 0; i < 2; ++i)
                edlibFreeAlignResult(align_rslt_left[i]);

            post_primer_time = current_time_ms();
            primer_used_time += (post_primer_time - pre_primer_time);

            if (primer_valid_flag)
            {
                pre_corr_time = current_time_ms();

                payload_len = end_loc - start_loc + 1;
                if (payload_len <= 0)
                    continue;
                tot_payload_base += payload_len;

                ++tot_reads_primer_valid;

                memset(payload_noisy_base, 0, MAX_PAYLOAD_LEN);
                memset(payload_noisy, 0, MAX_PAYLOAD_LEN);
                memset(compositeCode_noisy, 0, MAX_PAYLOAD_LEN);
                memset(compositeCode_AC, 0, MAX_PAYLOAD_LEN);
                memcpy(payload_noisy_base, read_base + start_loc, payload_len);
                if (!forward_flag)
                    reverseComplement(payload_noisy_base, payload_len); // reverse complementary

                for (int ti = 0; ti < payload_len; ++ti)
                {
                    if (*(payload_noisy_base + ti) == 'A')
                    {
                        *(payload_noisy + ti) = '0';
                        *(compositeCode_noisy + ti) = 1;
                    }
                    else if (*(payload_noisy_base + ti) == 'T')
                    {
                        *(payload_noisy + ti) = '0';
                        *(compositeCode_noisy + ti) = -1;
                    }
                    else if (*(payload_noisy_base + ti) == 'G')
                    {
                        *(payload_noisy + ti) = '1';
                        *(compositeCode_noisy + ti) = 1;
                    }
                    else if (*(payload_noisy_base + ti) == 'C')
                    {
                        *(payload_noisy + ti) = '1';
                        *(compositeCode_noisy + ti) = -1;
                    }
                    else // 'N'
                    {
                        *(payload_noisy + ti) = 'E';
                        *(compositeCode_noisy + ti) = 0;
                    }
                    *(compositeCode_AC + ti) = *(compositeCode_noisy + ti) == 1 ? 'A' : 'C';
                }

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

                if (valid_flag == 1)
                {
                    ++tot_reads_valid_corr;

                    for (int ti = 0; ti < payload_polished_len; ++ti)
                    {
                        if (recIndex + ti >= data_size || recIndex + ti < 0 || *(payload_polished + ti) == 'E')
                            continue;

                        ++*(*(unique_pos + (*(payload_polished + ti) - 48)) + recIndex + ti);
                    }
                }

                post_corr_time = current_time_ms();
                corr_used_time += (post_corr_time - pre_corr_time);

                if (mode && valid_flag == 0)
                {
                    //--- Payload DNA
                    fprintf(fpw_dat, "%s\n", payload_noisy);

                    //--- Noisy composite code ('A'+'C')
                    fprintf(fpw_PN, "@%d\n", n_reads_BWA);
                    fprintf(fpw_PN, "%s\n", compositeCode_AC);
                    fprintf(fpw_PN, "+\n");
                    fwrite(quality, 1, payload_len, fpw_PN);
                    fprintf(fpw_PN, "\n");

                    ++n_reads_BWA;
                }
            }
        }
        fclose(fpr_read);
        printf("Loaded %d usable reads from file-%d!\n", n_read, tot_block + 1);

        tot_reads += (size_t)n_read;

        if (mode)
        {
            fclose(fpw_dat);
            fclose(fpw_PN);

            tot_read_BWA += (size_t)n_reads_BWA;
        }

        pre_align_time = current_time_ms();
        if (mode && n_reads_BWA > 0)
        {
            printf("Remain reads for BWA align: %d\n", n_reads_BWA);

            sprintf(cmdLine, "bwa mem -v 1 -t %d %s %s -o %s", n_thread, refer_file_PNseq, filepath_PNseq, filepath_SAM);
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
                            ++tot_reads_valid_align;

                            fscanf(fpr_sam, "%*s %ld %*s %s %*[^\n]s", &startPos, cigar);
                            startPos -= 1;
                            startPos -= start_phase;
                            fscanf(fpr_dat, "%s\n", payload_noisy);
                            payload_len = strlen(payload_noisy);

                            n_bp = 0;
                            op_num = 0;
                            payload_polished_len = 0;
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
                                            *(payload_polished + payload_polished_len) = 'E';
                                            payload_polished_len++;
                                        }
                                    }
                                    else if (tmp_symbol == 'M')
                                    {
                                        for (int ti = 0; ti < op_num; ++ti)
                                        {
                                            *(payload_polished + payload_polished_len) = *(payload_noisy + n_bp);
                                            payload_polished_len++;
                                            n_bp++;
                                        }
                                    }
                                    op_num = 0; // reset
                                }
                            }

                            for (int ti = 0; ti < payload_polished_len; ++ti)
                            {
                                if (startPos + ti >= data_size || startPos + ti < 0 || *(payload_polished + ti) == 'E')
                                    continue;

                                ++*(*(unique_pos + (*(payload_polished + ti) - 48)) + startPos + ti);
                            }
                        }
                        else
                            fscanf(fpr_sam, "%*[^\n]%*c");
                    }
                }
            }
            fclose(fpr_dat);
            fclose(fpr_sam);
        }
        post_align_time = current_time_ms();
        align_used_time += post_align_time - pre_align_time;

        tot_reads_valid = tot_reads_valid_corr + tot_reads_valid_align;
        read_coverage = 1.0 * tot_reads / OLIGO_NUM;
        payload_coverage = 1.0 * tot_payload_base / data_size;
        printf("Total processed reads = %lu\n", tot_reads);
        printf("Reads with valid primers = %lu (%.2fx)\n", tot_reads_primer_valid, payload_coverage);
        printf("Reads with valid indices = %lu\n", tot_reads_valid);
        printf("  1) correlation = %lu (%.2f)\n", tot_reads_valid_corr, 1.0 * tot_reads_valid_corr / tot_reads_primer_valid * 100);
        printf("  2) alignment   = %lu (%.2f)\n", tot_reads_valid_align, 1.0 * tot_reads_valid_align / tot_reads_primer_valid * 100);

        printf("\n----------------Step2: Bit-wise majority voting----------------\n");
        //---- Bit-wise majority voting
        cons_ctx.tid = 0;
        cons_ctx.n_era_consensus = 0;
        cons_ctx.avail_base_num = 0;
        for (int tid = 0; tid < n_thread_on; ++tid)
            thpool_add_work(thpool, ConsensusWorker, (void *)&cons_ctx);
        thpool_wait(thpool);
        n_era_consensus = cons_ctx.n_era_consensus;
        avail_base_num = cons_ctx.avail_base_num;

        era_rate_consensus = 1.0 * n_era_consensus / data_size;
        avail_base_ratio = 1.0 * avail_base_num / data_size;

        printf("Available base ratio: %.2f\n", avail_base_ratio * 100);

        //---- Error correction using LDPC
        if (avail_base_ratio >= available_ratio)
        {
            printf("\n--------------------Step3: Error correction--------------------\n");
            decoder_ctx.tid = 0;
            decoder_ctx.n_corrupted_code = 0;

            for (int tid = 0; tid < n_thread_on; ++tid)
                thpool_add_work(thpool, LDPCDecodingWorker, (void *)&decoder_ctx);
            thpool_wait(thpool);

#ifdef DATA_RECOVERY
            n_byte = 0;
            for (int blk = 0; blk < BLOCKCODE_NUM; ++blk)
            {
                for (int ni = 0; ni < LDPC_K / 8; ++ni)
                {
                    pbuff_ptr = *(rec_payload + blk) + ni * 8;
                    *(rec_data_bytes + n_byte) = bit2ascii(pbuff_ptr);
                    n_byte++;
                }
            }

            if ((fpw = fopen(img1_filename, "wb")) == NULL)
            {
                fprintf(stderr, "Faile to open %s\n", img1_filename);
                return 1;
            }
            fwrite(rec_data_bytes, 1, IMG1_BYTE, fpw);
            fclose(fpw);

            if ((fpw = fopen(img2_filename, "wb")) == NULL)
            {
                fprintf(stderr, "Faile to open %s\n", img2_filename);
                return 1;
            }
            fwrite(rec_data_bytes + IMG1_BYTE, 1, IMG2_BYTE, fpw);
            fclose(fpw);

#if 1 // !!!!!!Note: show image
            if (!open_flag)
            {
                sprintf(cmdLine, "eog %s 2> /dev/null &"
                                 "eog %s 2> /dev/null &",
                        img1_filename, img2_filename);
                status = system(cmdLine);

                sprintf(cmdLine,
                        "sleep 0.5;"
                        "wmctrl -r \"img1.jpg\" -e 0,900,640,540,360;"
                        "wmctrl -r \"img2.jpg\" -e 0,1350,690,540,360");
                status = system(cmdLine);

                open_flag = 1;
            }
#endif
#endif

            if (!decoder_ctx.n_corrupted_code)
                break;
        }

        ++tot_block;
    }

    post_time = current_time_ms();
    double duration_run = difftime(post_time, pre_time) / 1000;

    printf("\n++++++++++++++++++++++++++++Summary++++++++++++++++++++++++++++\n");
    printf("Total processed reads = %lu\n", tot_reads);
    printf("Reads with valid primers = %lu (%.2fx)\n", tot_reads_primer_valid, payload_coverage);
    printf("Reads with valid indices = %lu\n", tot_reads_valid);
    printf("  1) correlation = %lu (%.2f)\n", tot_reads_valid_corr, 1.0 * tot_reads_valid_corr / tot_reads_primer_valid * 100);
    printf("  2) alignment   = %lu (%.2f)\n", tot_reads_valid_align, 1.0 * tot_reads_valid_align / tot_reads_primer_valid * 100);
    printf("Available base ratio: %.2f\n", avail_base_ratio * 100);

    // Error evaluation. 0: off, 1: on
#if 0
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
        assert(cur_len == OLIGO_LEN);

        memcpy(dna_refer + PAYLOAD_LEN * (nn++), seq_buff + PRIMER_LEN, PAYLOAD_LEN);
    }
    fclose(fpr_refer);

    unsigned erasure_num = 0, error_num = 0;
    double erasure_rate = 0.0, error_rate = 0.0;
    unsigned error_num_block[BLOCKCODE_NUM] = {0};
    unsigned erasure_num_block[BLOCKCODE_NUM] = {0};
    double error_rate_block[BLOCKCODE_NUM] = {0};
    double erasure_rate_block[BLOCKCODE_NUM] = {0};
    int block_id = 0;

    char bit1, bit2;
    for (size_t nn = 0; nn < data_size; ++nn)
    {
        bit1 = *(consensus_codewords + nn);
        bit2 = (*(dna_refer + nn) == 'A' || *(dna_refer + nn) == 'T') ? '0' : '1';

        if (bit1 != bit2)
        {
            if (bit1 == 'E' || bit2 == 'E')
            {
                ++erasure_num; // erasure
                ++*(erasure_num_block + block_id);
            }
            else
            {
                ++error_num; // error
                ++*(error_num_block + block_id);
            }
        }

        if ((nn + 1) % LDPC_N == 0)
        {
            *(error_rate_block + block_id) = 1.0 * *(error_num_block + block_id) / LDPC_N;
            *(erasure_rate_block + block_id) = 1.0 * *(erasure_num_block + block_id) / LDPC_N;
            block_id++;
        }
    }
    error_rate = 1.0 * error_num / data_size;
    erasure_rate = 1.0 * erasure_num / data_size;

    printf("-------------------------After consensus-------------------------\n");
    printf("Bit error num = %u, bit error rate = %.8f\n", error_num, error_rate);
    printf("Bit erasure num = %u, bit erasure rate = %.8f\n", erasure_num, erasure_rate);

    char *srcdata_refer = (char *)alloc_init(BLOCKCODE_NUM * LDPC_K, sizeof(char));
    if ((fpr_refer = fopen(srcdata_file, "rb")) == NULL)
    {
        fprintf(stderr, "Could not open file %s!\n", srcdata_file);
        return 1;
    }
    size_t rt1 = fread(srcdata_refer, 1, BLOCKCODE_NUM * LDPC_K, fpr_refer);
    assert(rt1 == BLOCKCODE_NUM * LDPC_K);
    fclose(fpr_refer);

    unsigned bit_error_num = 0;
    unsigned bit_error_num_block[BLOCKCODE_NUM] = {0};
    double bit_error_rate_block[BLOCKCODE_NUM] = {0.0};

    for (int blk = 0; blk < BLOCKCODE_NUM; ++blk)
    {
        for (int nb = 0; nb < LDPC_K; ++nb)
        {
            if (*(*(rec_payload + blk) + nb) != *(srcdata_refer + LDPC_K * blk + nb) - 48)
            {
                ++bit_error_num;
                ++*(bit_error_num_block + blk);
            }
        }

        *(bit_error_rate_block + blk) = 1.0 * *(bit_error_num_block + blk) / LDPC_K;
    }
    double bit_error_rate = 1.0 * bit_error_num / (BLOCKCODE_NUM * LDPC_K);
    int dec_flag = bit_error_num == 0 ? 0 : bit_error_num;

    printf("----------------------------After ECC----------------------------\n");
    printf("Bit error num = %u, Bit error rate = %.6f\n", bit_error_num, bit_error_rate);
    printf("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");

#endif

    thpool_destroy(thpool);
    pthread_mutex_destroy(&decoder_ctx.mutex);
    pthread_mutex_destroy(&cons_ctx.mutex);

    time(&end_time);
    strftime(timebuff, 100, "%Y-%m-%d %H:%M:%S ", localtime(&end_time));
    printf("\nEnd running: %s\n", timebuff);
    printf("Elapsed time: %.2lf sec\n\n", duration_run);
    return 0;
}

void ConsensusWorker(void *args)
{
    ConsensuContext *cons_ctx = (ConsensuContext *)args;
    pthread_mutex_lock(&cons_ctx->mutex);
    int tid = cons_ctx->tid;
    ++cons_ctx->tid;
    pthread_mutex_unlock(&cons_ctx->mutex);

    unsigned bits_offset = cons_ctx->process_bits * tid;
    unsigned process_bits = (cons_ctx->n_jobs == tid + 1) ? cons_ctx->data_size - bits_offset
                                                          : cons_ctx->process_bits;
    int cof0, cof1;
    int confidence;

    double *lratio_codewords = cons_ctx->lratio_codewords + bits_offset;
    char *consensus_codewords = cons_ctx->consensus_codewords + bits_offset;
    unsigned short **unique_pos = cons_ctx->unique_pos;

    unsigned avail_base_num = 0;
    unsigned n_era_consensus = 0;

    if (process_bits <= 0)
        return;

    for (unsigned sb = 0; sb < process_bits; ++sb)
    {
        cof0 = (int)*(*(unique_pos + 0) + bits_offset + sb);
        cof1 = (int)*(*(unique_pos + 1) + bits_offset + sb);

        if (cof0 + cof1 == 0)
        {
            *(consensus_codewords + sb) = 'E';
            confidence = 0;
        }
        else
        {
            if (cof0 == cof1)
            {
                *(consensus_codewords + sb) = 'E';
                confidence = 0;
                ++n_era_consensus;
            }
            else if (cof0 > cof1)
            {
                ++avail_base_num;

                *(consensus_codewords + sb) = '0';
                confidence = cof0 - cof1;
            }
            else
            {
                ++avail_base_num;

                *(consensus_codewords + sb) = '1';
                confidence = cof0 - cof1;
            }
        }

        *(lratio_codewords + sb) = exp(-2.0 * confidence);
    }

    pthread_mutex_lock(&cons_ctx->mutex);
    cons_ctx->n_era_consensus += n_era_consensus;
    cons_ctx->avail_base_num += avail_base_num;
    pthread_mutex_unlock(&cons_ctx->mutex);
}

void LDPCDecodingWorker(void *args)
{
    DecoderContext *decoder_ctx = (DecoderContext *)args;

    pthread_mutex_lock(&decoder_ctx->mutex);
    int tid = decoder_ctx->tid;
    ++decoder_ctx->tid;

    pthread_mutex_unlock(&decoder_ctx->mutex);

    unsigned code_offset = decoder_ctx->process_code_num * tid;
    int process_code_num = (decoder_ctx->n_jobs == tid + 1) ? BLOCKCODE_NUM - code_offset
                                                            : decoder_ctx->process_code_num;
    if (process_code_num <= 0)
        return;

    char *mSeq1 = decoder_ctx->mSeq1;
    char *mSeq2 = decoder_ctx->mSeq2;

    double *lratio_codewords = decoder_ctx->lratio_codewords + LDPC_N * code_offset;
    unsigned *permutation = decoder_ctx->permutation;
    char **rec_payload = decoder_ctx->rec_payload + code_offset;

    double *lratio = (double *)alloc_init(LDPC_N, sizeof(double));

    LDPCDecode *ldpc_decode = new LDPCDecode();
    if (ldpc_decode->decode_init(ldpc_pchk_file))
    {
        fprintf(stderr, "Failed to initialize LDPC decoder!\n");
        return;
    }

    int ldpc_chk;
    int n_corrupted_code = 0;
    for (int blk = 0; blk < process_code_num; ++blk)
    {
        if (*(decoder_ctx->codeword_status + code_offset + blk) == 0)
        {
            //// de-permuation
            for (int nb = 0; nb < LDPC_N; ++nb)
            {
                *(lratio + *(permutation + nb)) = *(lratio_codewords + LDPC_N * blk + nb);
            }

            //---- LDPC(64800, 54000) decoding
            ldpc_chk = ldpc_decode->decode_LDPC_soft(*(rec_payload + blk), lratio, maxIteration_ldpc);

            if (ldpc_chk)
                ++n_corrupted_code;
            else
                *(decoder_ctx->codeword_status + code_offset + blk) = 1;

            //---- descramble
            for (int nb = 0; nb < LDPC_K; ++nb)
            {
                *(*(rec_payload + blk) + nb) = (*(*(rec_payload + blk) + nb) +
                                                *(mSeq1 + ((unsigned)LDPC_K * (code_offset + blk) + nb) % MSEQ1) +
                                                *(mSeq2 + ((unsigned)LDPC_K * (code_offset + blk) + nb) % MSEQ2)) %
                                               2;
            }
        }
    }

    if (n_corrupted_code)
    {
        pthread_mutex_lock(&decoder_ctx->mutex);
        decoder_ctx->n_corrupted_code += n_corrupted_code;
        pthread_mutex_unlock(&decoder_ctx->mutex);
    }

    free(lratio);
    ldpc_decode->mem_free_dec();
    delete ldpc_decode;
    ldpc_decode = NULL;
}

int isFileReady(const char *filename, size_t set_lines)
{
    int status;
    size_t lines = 0;
    char cmdLine[1000] = {0};
    char buff[100] = {0};

    struct stat file_stat;
    if (stat(filename, &file_stat) == -1)
        return 0;

    sprintf(cmdLine, "awk 'END{print NR}' %s", filename);

    FILE *pipe;
    pipe = popen(cmdLine, "r");
    if (pipe == NULL)
        return 0;

    if (fgets(buff, sizeof(buff), pipe) != NULL)
        lines = atoi(buff);
    // printf("%lu\n", lines);

    status = pclose(pipe);

    if (lines == set_lines * 4) // !!! NOTE FSATQ
        return 1;
    else
        return 0;
}
