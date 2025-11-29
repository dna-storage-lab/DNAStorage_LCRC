/************************************
 * Copyright (C) 2025 Tianjin University (TJU)
 *
 * @brief: Real-time readout of medium-scale oligo pools via nanopore sequencing
 *
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
#include "rsCodec.h"

#include <pthread.h>
#include "thpool.h"

#define DATA_RECOVERY
#define FILE_NUM 4

#define OLIGO_NUM 299700
#define OLIGO_LEN 200
#define PAYLOAD_LEN 160
#define PRIMER_LEN 20

#define LDPC_K 54000
#define LDPC_N 64800
#define LDPC_M 10800
#define CHECK_DEGREE 20

#define RS_ORDER 10
#define RS_N 1023
#define RS_K 1018
#define RS_N_SHORTEN 740
#define RS_K_SHORTEN 735

#define MSEQ 65535

#define MAX_COMPONENTS 5 // Number of component sequences
#define MAX_THREAD_FILE 100

#define MAX_READ_LEN 1000
#define MAX_PAYLOAD_LEN 1000

#define N_TYPE 2
const char element[N_TYPE] = {'0', '1'};

typedef struct
{
    int n_component;
    uint64_t period;
    uint64_t start_phase;

    uint64_t *M;
    uint64_t *M_sup;

    unsigned *cLen;
    int *cPhase;
    int **c;

} Configuration;

typedef struct
{
    unsigned n_reads;
    unsigned n_reads_BWA;

    unsigned n_reads_primer_valid;
    unsigned n_reads_valid_corr;
    unsigned n_reads_valid_align;
    unsigned n_read_base;
    unsigned n_payload_base;

    double primer_used_time;
    double corr_used_time;
    double align_used_time;

    unsigned short *unique_pos[N_TYPE];
} CorrelationContext;

typedef struct
{
    int tid;

    int n_threads_bwa;

    int chip1_cnt_batch;
    int chip2_cnt_batch;

    unsigned data_size;

    bool mode;
    double threshold_primer;
    double threshold_autoCorr;

    char *log_dir;
    char *chip1_reads_file_prefix;
    char *chip2_reads_file_prefix;

    CorrelationContext *corrCtx;

    pthread_mutex_t mutex;

    Configuration *config;

} WorkerHandler;

typedef struct
{
    int tid_ldpc;    
    int n_jobs_ldpc; 
    int tot_ldpc;
    int n_ldpc_per_tid; // LDPC num per thread

    int tid_rs;   
    int n_jobs_rs; 
    int tot_rs;
    int n_rs_per_tid; // RS num per thread

    unsigned *n_corrupted_ldpc; // corrupted ldpc per tid
    unsigned *n_corrupted_rs;   // corrupted rs per tid

    char *codeword_status;
    char **bit_status;

    double **lratio_codewords;
    char **rec_codeword_errfree;
    char **product_codeword;

    pthread_mutex_t mutex;

} DecoderContext;

typedef struct
{
    int tid;
    int n_jobs;

    int n_tid_file;

    unsigned data_size;
    int process_bits;

    unsigned avail_base_num;
    unsigned n_era_consensus;

    char *consensus_codewords;
    double *lratio_codewords; 

    CorrelationContext *corrCtx;

    pthread_mutex_t mutex;

} ConsensuContext;

const char *primerL = "AATCATGGCCTTCAAACCGT"; 
const char *primerR = "AACGCTCCGAAAGTCTTGTT"; 

const char *primerL_revComp = "ACGGTTTGAAGGCCATGATT";
const char *primerR_revComp = "AACAAGACTTTCGGAGCGTT";

const char *ldpc_pchk_file = "config/ldpc.pchk";

const char *refer_file_PNseq = "config/PNseq_refer.fa";
const char *tmp_ofile_data_BWA = "tmp_read_BWA";
const char *tmp_ofile_PNseq_BWA = "tmp_PNseq_BWA";
const char *tmp_ofile_SAM_BWA = "tmp_SAM_BWA";

const char *imagfile_list[] = {"Tundra.jpg", "Lake.jpg", "Galaxy.jpg", "Poem_EN.txt"};
const unsigned filesize_bytes[] = {1657329, 1710060, 1591676, 2159};

// Reference files for recovery evaluation
const char *refer_file = "reference/oligoPool.fa";  // reference FASTA

char bit2ascii(const char *inBit)
{
    char symbol = 0;
    for (int i = 0; i < 8; i++)
    {
        symbol += *(inBit + i) * (1 << (7 - i));
    }
    return symbol;
}

void pipelineWorker_chip1(void *args);
void pipelineWorker_chip2(void *args);
void correlation_alignment(WorkerHandler *worker_hdlr, int tid, const char *reads_file);

void ConsensusWorker(void *args);
void LDPCDecodingWorker(void *args);
void RSDecodingWorker(void *args);

size_t current_time_ms()
{
    struct timeval time_now;
    gettimeofday(&time_now, NULL);
    return (size_t)time_now.tv_sec * 1000 + time_now.tv_usec / 1000;
}

int isFileReady(const char *filename, size_t set_lines);


// Max iteration for ldpc decoding
const int maxIteration_ldpc = 40;

const double threshold_coverage = 3.0;

const int n_threads_bwa = 10;
const int target_files = 15;

int main(int argc, char *argv[])
{
    setvbuf(stdout, NULL, _IONBF, 0); // none buffer stdout

    char result_dir[500] = {0};
    char log_dir[500] = {0};
    char chip1_reads_file_prefix[500] = {0};
    char chip2_reads_file_prefix[500] = {0};

    int status;
    char cmdLine[1000] = {0};

    const int n_component = 5;

    int target_reads = -1;

    int n_threads = 1;
    int n_threads_file = 4;
    int mode = 1;           // 0: correlation-based, 1: alignment-based

    double threshold_primer = 0.2;
    double threshold_autoCorr = 0.8;

    double available_ratio = 0.945;

    unsigned cLen[MAX_COMPONENTS] = {31, 35, 43, 47, 59};
    int cPhase_init[MAX_COMPONENTS] = {0};

    if (argc < 10)
    {
        fprintf(stderr, "Usage: "
                        "%s [result dir] [log dir] [Chip 1 fastq file] [Chip 2 fastq file] [File threads] [Threads]"
                        "[Correlation threshold] [Mode] [Base available threshold] [Target reads per file]\n",
                argv[0]);
        return 1;
    }
    sscanf(argv[1], "%s", result_dir);
    sscanf(argv[2], "%s", log_dir);
    sscanf(argv[3], "%s", chip1_reads_file_prefix);
    sscanf(argv[4], "%s", chip2_reads_file_prefix);
    sscanf(argv[5], "%d", &n_threads_file);
    sscanf(argv[6], "%d", &n_threads);
    sscanf(argv[7], "%lf", &threshold_autoCorr);
    sscanf(argv[8], "%d", &mode);
    sscanf(argv[9], "%lf", &available_ratio);
    assert(n_component <= MAX_COMPONENTS);

    if (argc >= 11)
        sscanf(argv[10], "%d", &target_reads);
    assert(target_reads > 0);

    time_t start_time, end_time;
    char timebuff[100];
    time(&start_time);
    strftime(timebuff, 100, "%Y-%m-%d %H:%M:%S ", localtime(&start_time));
    printf("\nStart running: %s\n\n", timebuff);

    /// RS(740, 735), GF(2^10)
    int rs_2t = RS_N - RS_K;
    int rs_n_shorten = RS_N_SHORTEN;
    int rs_k_shorten = rs_n_shorten - rs_2t;

    unsigned ldpc_rs_k = rs_k_shorten * LDPC_K;
    unsigned ldpc_rs_n = rs_n_shorten * LDPC_N;
    unsigned num_ldpc = rs_n_shorten;   
    unsigned num_rs = LDPC_K / RS_ORDER;

    unsigned data_size = ldpc_rs_n;

    printf("++++++++++++Configuration++++++++++++\n");
    printf("Thread number = %d\n", n_threads);
    if (mode == 1)
        printf("Mode = Alignment-based\n");
    else
        printf("Mode = Correlation-based\n");
    printf("Double check threshold = %.2f\n", threshold_autoCorr);
    printf("Available base ratio threshold = %.4f\n", available_ratio);

    int tot_batch = 0;
    int chip1_cnt_batch = 0;
    int chip2_cnt_batch = 0;
    int num_working = 0;
    int num_working_chip1 = 0, num_working_chip2 = 0;
    int num_decoding_file = 0;

    int n_reads = 0;
    int n_reads_BWA = 0;

    unsigned tot_reads = 0;
    unsigned tot_read_BWA = 0;

    double read_coverage = 0.0;
    double payload_coverage = 0.0;

    unsigned avail_base_num = 0;
    double avail_base_ratio = 0.0;

    unsigned n_era_consensus = 0;
    double era_rate_consensus = 0.0;

    size_t tot_read_base = 0;
    size_t tot_payload_base = 0;

    unsigned tot_reads_primer_valid = 0;
    unsigned tot_reads_valid = 0;
    unsigned tot_reads_valid_corr = 0;
    unsigned tot_reads_valid_align = 0;

    char reads_file[500] = {0};

    char codeword_status[RS_N_SHORTEN] = {0}; // 0: success, 1:fail

    char *consensus_codewords = (char *)alloc_init(data_size, sizeof(char));
    double *lratio_codewords_permutate = (double *)alloc_init(data_size, sizeof(double));

    double **lratio_codewords = (double **)alloc_init(RS_N_SHORTEN, sizeof(double *));
    char **bit_status = (char **)malloc(RS_N_SHORTEN * sizeof(char *));           // 0: success, 1:fail
    char **rec_codeword_errfree = (char **)malloc(RS_N_SHORTEN * sizeof(char *)); 
    char **product_codeword = (char **)malloc(RS_N_SHORTEN * sizeof(char *));
    for (int i = 0; i < RS_N_SHORTEN; ++i)
    {
        *(lratio_codewords + i) = (double *)alloc_init(LDPC_N, sizeof(double));
        *(bit_status + i) = (char *)alloc_init(LDPC_N, sizeof(char));
        *(rec_codeword_errfree + i) = (char *)alloc_init(LDPC_N, sizeof(char));
        *(product_codeword + i) = (char *)alloc_init(LDPC_N, sizeof(char));
    }

    char *rec_payload = (char *)alloc_init(data_size, sizeof(char));

    CorrelationContext corrCtx[MAX_THREAD_FILE];
    for (int i = 0; i < n_threads_file; ++i)
    {
        for (int j = 0; j < N_TYPE; ++j)
            *((corrCtx + i)->unique_pos + j) = (unsigned short *)alloc_init(data_size, sizeof(unsigned short));
    }

    //----- Load m sequence (scrambler)
    unsigned mSeqLen = MSEQ;
    char mSeq[MSEQ];
    FILE *fpr = NULL;
    if ((fpr = fopen("config/m65535.txt", "rb")) == NULL)
    {
        fprintf(stderr, "Could not open file config/m65535.txt!\n");
        return 1;
    }
    fread(mSeq, sizeof(char), mSeqLen, fpr);
    fclose(fpr);
    for (int i = 0; i < mSeqLen; i++)
        *(mSeq + i) -= '0';

    threadpool thpool = thpool_init(n_threads);

    uint64_t period = 1;
    for (int i = 0; i < n_component; i++)
        period *= *(cLen + i);

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

    Configuration config;
    config.n_component = n_component;
    config.cLen = cLen;
    config.c = c;
    config.period = period;
    config.M = M;
    config.M_sup = M_sup;
    config.start_phase = start_phase;
    config.cPhase = cPhase_init;

    WorkerHandler worker_hdlr;
    worker_hdlr.data_size = data_size;

    worker_hdlr.mode = mode;
    worker_hdlr.threshold_primer = threshold_primer;
    worker_hdlr.threshold_autoCorr = threshold_autoCorr;

    worker_hdlr.chip1_reads_file_prefix = chip1_reads_file_prefix;
    worker_hdlr.chip2_reads_file_prefix = chip2_reads_file_prefix;
    worker_hdlr.log_dir = log_dir;

    worker_hdlr.corrCtx = corrCtx;
    worker_hdlr.config = &config;

    worker_hdlr.n_threads_bwa = n_threads_bwa;

    pthread_mutex_init(&worker_hdlr.mutex, NULL);

    DecoderContext decoder_ctx;
    decoder_ctx.lratio_codewords = lratio_codewords;
    decoder_ctx.codeword_status = codeword_status;
    decoder_ctx.bit_status = bit_status;
    decoder_ctx.rec_codeword_errfree = rec_codeword_errfree;
    decoder_ctx.product_codeword = product_codeword;
    decoder_ctx.tot_ldpc = num_ldpc;
    decoder_ctx.n_ldpc_per_tid = ceil(1.0 * num_ldpc / n_threads);
    int n_threads_on = n_threads * decoder_ctx.n_ldpc_per_tid > num_ldpc
                           ? ceil(1.0 * num_ldpc / decoder_ctx.n_ldpc_per_tid)
                           : n_threads;
    decoder_ctx.n_jobs_ldpc = n_threads_on;
    decoder_ctx.tot_rs = num_rs;
    decoder_ctx.n_rs_per_tid = ceil(1.0 * num_rs / n_threads);
    decoder_ctx.n_jobs_rs = n_threads * decoder_ctx.n_rs_per_tid > num_rs
                                ? ceil(1.0 * num_rs / decoder_ctx.n_rs_per_tid)
                                : n_threads;
    pthread_mutex_init(&decoder_ctx.mutex, NULL);

    unsigned tot_corrupted_ldpc = 0;
    unsigned tot_corrupted_rs = 0;
    unsigned *n_corrupted_ldpc = (unsigned *)alloc_init(decoder_ctx.n_jobs_ldpc, sizeof(unsigned));
    unsigned *n_corrupted_rs = (unsigned *)alloc_init(decoder_ctx.n_jobs_rs, sizeof(unsigned));
    decoder_ctx.n_corrupted_ldpc = n_corrupted_ldpc;
    decoder_ctx.n_corrupted_rs = n_corrupted_rs;

    ConsensuContext cons_ctx;
    cons_ctx.data_size = data_size;
    cons_ctx.consensus_codewords = consensus_codewords;
    cons_ctx.lratio_codewords = lratio_codewords_permutate;
    cons_ctx.process_bits = decoder_ctx.n_ldpc_per_tid * LDPC_N;
    cons_ctx.n_jobs = n_threads_on;

    cons_ctx.n_tid_file = n_threads_file;
    cons_ctx.corrCtx = corrCtx;
    pthread_mutex_init(&cons_ctx.mutex, NULL);

#ifdef DATA_RECOVERY
    //---- Data recovery ----//
    size_t tot_bytes = ldpc_rs_k / 8;
    char *pbuff_ptr = NULL;
    char *rec_data_bytes = (char *)alloc_init(tot_bytes, sizeof(char));
    unsigned offset_byte = 0;

    FILE *fpw_rec = NULL;
    char img_filename[500] = {0};

    int open_flag = 0;
#endif

    time_t pre_corr_time, post_corr_time;
    time_t pre_align_time, post_align_time;
    time_t pre_cons_time, post_cons_time;
    time_t pre_dec_time, post_dec_time;

    double primer_used_time = 0.0;
    double corr_used_time = 0.0;
    double align_used_time = 0.0;
    double cons_used_time = 0.0;
    double dec_used_time = 0.0;

    printf("\n----------------Start read-by-read processing------------------\n");
    time_t pre_time, post_time;
    time_t init_time, now_time;
    pre_time = current_time_ms();

    while (1)
    {
        printf("\n--------------Step1: Primer & index identification-------------\n");
        init_time = current_time_ms();
        num_working = 0;
        num_working_chip1 = 0;
        num_working_chip2 = 0;
        worker_hdlr.tid = 0;
        while (num_working < n_threads_file)
        {
            sprintf(reads_file, "%s%d.fastq", chip1_reads_file_prefix, chip1_cnt_batch);
            if (isFileReady(reads_file, target_reads))
            {
                thpool_add_work(thpool, pipelineWorker_chip1, (void *)&worker_hdlr);
                ++chip1_cnt_batch;
                ++num_working;
                ++num_working_chip1;
            }
            if (num_working == n_threads_file)
                break;


            sprintf(reads_file, "%s%d.fastq", chip2_reads_file_prefix, chip2_cnt_batch);
            if (isFileReady(reads_file, target_reads))
            {
                thpool_add_work(thpool, pipelineWorker_chip2, (void *)&worker_hdlr);
                ++chip2_cnt_batch;
                ++num_working;
                ++num_working_chip2;
            }

            now_time = current_time_ms();
            if (now_time - init_time > 5000) // 5 sec detection
            {
                init_time = current_time_ms();
                if (num_working > 0)
                {
                    printf("Timeout for 5 sec, start processing!\n");
                    break;
                }
                else
                    printf("Timeout for 5 sec, waiting for sequencing data!\n");
            }
        }
        tot_batch += num_working;
        if (num_working_chip1)
            printf("Loaded %d fastq files from Chip-1: file-%d to file-%d!\n",
                   num_working_chip1, chip1_cnt_batch - num_working_chip1 + 1, chip1_cnt_batch);
        if (num_working_chip2)
            printf("Loaded %d fastq files from Chip-2: file-%d to file-%d!\n",
                   num_working_chip2, chip2_cnt_batch - num_working_chip2 + 1, chip2_cnt_batch);
        thpool_wait(thpool);

        n_reads = 0;
        n_reads_BWA = 0;
        for (int ni = 0; ni < num_working; ++ni)
        {
            n_reads += (worker_hdlr.corrCtx + ni)->n_reads;
            n_reads_BWA += (worker_hdlr.corrCtx + ni)->n_reads_BWA;

            tot_reads_primer_valid += (worker_hdlr.corrCtx + ni)->n_reads_primer_valid;
            tot_reads_valid_corr += (worker_hdlr.corrCtx + ni)->n_reads_valid_corr;
            tot_reads_valid_align += (worker_hdlr.corrCtx + ni)->n_reads_valid_align;

            tot_read_base += (worker_hdlr.corrCtx + ni)->n_read_base;
            tot_payload_base += (worker_hdlr.corrCtx + ni)->n_payload_base;

            primer_used_time += (worker_hdlr.corrCtx + ni)->primer_used_time;
            corr_used_time += (worker_hdlr.corrCtx + ni)->corr_used_time;
            align_used_time += (worker_hdlr.corrCtx + ni)->align_used_time;
        }
        tot_reads += n_reads;
        tot_read_BWA += n_reads_BWA;

        tot_reads_valid = tot_reads_valid_corr + tot_reads_valid_align;
        read_coverage = 1.0 * tot_reads / OLIGO_NUM;
        payload_coverage = 1.0 * tot_payload_base / data_size;
        printf("Processed %d usable reads from %d fastq files!\n", n_reads, num_working);
        printf("Total processed reads = %u (%d files)\n", tot_reads, tot_batch);
        printf("Reads with valid primers = %u (%.2fx)\n", tot_reads_primer_valid, payload_coverage);
        printf("Reads with valid indices = %u\n", tot_reads_valid);
        printf("  1) correlation = %u (%.2f)\n", tot_reads_valid_corr, 1.0 * tot_reads_valid_corr / tot_reads_primer_valid * 100);
        printf("  2) alignment   = %u (%.2f)\n", tot_reads_valid_align, 1.0 * tot_reads_valid_align / tot_reads_primer_valid * 100);

        if (payload_coverage >= threshold_coverage)
        {
            printf("\n----------------Step2: Bit-wise majority voting----------------\n");

            pre_cons_time = current_time_ms();

            cons_ctx.tid = 0;
            cons_ctx.n_era_consensus = 0;
            cons_ctx.avail_base_num = 0;
            for (int tid = 0; tid < n_threads_on; ++tid)
                thpool_add_work(thpool, ConsensusWorker, (void *)&cons_ctx);
            thpool_wait(thpool);

            post_cons_time = current_time_ms();
            cons_used_time += difftime(post_cons_time, pre_cons_time);

            n_era_consensus = cons_ctx.n_era_consensus;
            avail_base_num = cons_ctx.avail_base_num;

            era_rate_consensus = 1.0 * n_era_consensus / data_size;
            avail_base_ratio = 1.0 * avail_base_num / data_size;

            printf("Available base ratio: %.2f\n", avail_base_ratio * 100);

            if (avail_base_ratio >= available_ratio && (tot_batch - num_decoding_file) >= target_files)
            {
                num_decoding_file = tot_batch;

                printf("\n--------------------Step3: Error correction--------------------\n");
                pre_dec_time = current_time_ms();

                // de-interleaving
                for (int row = 0; row < RS_N_SHORTEN; ++row)
                    for (int col = 0; col < LDPC_N; ++col)
                    {
                        *(*(lratio_codewords + (row + col) % RS_N_SHORTEN) + col) = *(lratio_codewords_permutate + row * LDPC_N + col);
                    }

                // LDPC decoding
                decoder_ctx.tid_ldpc = 0;
                for (int tid = 0; tid < decoder_ctx.n_jobs_ldpc; ++tid)
                {
                    *(decoder_ctx.n_corrupted_ldpc + tid) = 0;
                    thpool_add_work(thpool, LDPCDecodingWorker, (void *)&decoder_ctx);
                }
                thpool_wait(thpool);

                tot_corrupted_ldpc = 0;
                for (int tid = 0; tid < decoder_ctx.n_jobs_ldpc; ++tid)
                    tot_corrupted_ldpc += *(decoder_ctx.n_corrupted_ldpc + tid);

                // RS decoding
                decoder_ctx.tid_rs = 0;
                for (int tid = 0; tid < decoder_ctx.n_jobs_rs; ++tid)
                {
                    *(decoder_ctx.n_corrupted_rs + tid) = 0;
                    thpool_add_work(thpool, RSDecodingWorker, (void *)&decoder_ctx);
                }
                thpool_wait(thpool);

                tot_corrupted_rs = 0;
                for (int tid = 0; tid < decoder_ctx.n_jobs_rs; ++tid)
                    tot_corrupted_rs += *(decoder_ctx.n_corrupted_rs + tid);

                fprintf(stderr, "\nResidual corrupted LDPC = %u/%u, corrupted RS = %u/%u\n",
                        tot_corrupted_ldpc, num_ldpc, tot_corrupted_rs, num_rs);

                // De-scrambling
                for (int row = 0; row < RS_K_SHORTEN; ++row)
                {
                    pbuff_ptr = *(product_codeword + row);
                    for (int col = 0; col < LDPC_K; ++col)
                        *(rec_payload + row * LDPC_K + col) = (*(pbuff_ptr + col) + *(mSeq + col)) % 2;
                }

#ifdef DATA_RECOVERY
                for (unsigned ni = 0; ni < tot_bytes; ++ni)
                    *(rec_data_bytes + ni) = bit2ascii(rec_payload + ni * 8);

                offset_byte = 0;
                for (int n = 0; n < FILE_NUM; ++n)
                {
                    sprintf(img_filename, "%s/%s", result_dir, imagfile_list[n]);
                    if ((fpw_rec = fopen(img_filename, "wb")) == NULL)
                    {
                        fprintf(stderr, "Failed to open %s\n", img_filename);
                        return 1;
                    }

                    fwrite(rec_data_bytes + offset_byte, 1, *(filesize_bytes + n), fpw_rec);
                    offset_byte += *(filesize_bytes + n);
                    fclose(fpw_rec);
                }

                post_dec_time = current_time_ms();
                dec_used_time += difftime(post_dec_time, pre_dec_time);
#endif

                // if (!tot_corrupted_ldpc && !tot_corrupted_rs)
                if (!tot_corrupted_rs)
                    break;
            }
        }
    }

    post_time = current_time_ms();
    double duration_run = difftime(post_time, pre_time) / 1000;

    printf("\n++++++++++++++++++++++++++++Summary++++++++++++++++++++++++++++\n");
    printf("Total processed fastq files = %d\n", tot_batch);
    printf("Total processed reads = %u\n", tot_reads);
    printf("Reads with valid primers = %u (%.2fx)\n", tot_reads_primer_valid, payload_coverage);
    printf("Reads with valid indices = %u\n", tot_reads_valid);
    printf("  1) correlation = %u (%.2f)\n", tot_reads_valid_corr, 1.0 * tot_reads_valid_corr / tot_reads_primer_valid * 100);
    printf("  2) alignment   = %u (%.2f)\n", tot_reads_valid_align, 1.0 * tot_reads_valid_align / tot_reads_primer_valid * 100);
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

        for (int i = 0; i < PAYLOAD_LEN; ++i)
            *(dna_refer + PAYLOAD_LEN * nn + i) = *(seq_buff + PRIMER_LEN + i);
        // memcpy(dna_refer + PAYLOAD_LEN * (nn++), seq_buff + PRIMER_LEN, PAYLOAD_LEN);
        ++nn;
    }
    fclose(fpr_refer);

    unsigned erasure_num = 0, error_num = 0;
    double erasure_rate = 0.0, error_rate = 0.0;

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
            }
            else
            {
                ++error_num; // error
            }
        }

        if ((nn + 1) % LDPC_N == 0)
        {
            block_id++;
        }
    }
    error_rate = 1.0 * error_num / (data_size - erasure_num);
    erasure_rate = 1.0 * erasure_num / data_size;

    printf("-------------------------After consensus-------------------------\n");
    printf("Bit error num = %u, bit error rate = %.8f\n", error_num, error_rate);
    printf("Bit erasure num = %u, bit erasure rate = %.8f\n", erasure_num, erasure_rate);
#endif

    printf("\nPrimer ident. time: %.2lf sec\n", primer_used_time / n_threads_file / 1000);
    printf("Correlation time: %.2lf sec\n", corr_used_time / n_threads_file / 1000);
    printf("Alignment time: %.2lf sec\n", align_used_time / n_threads_file / 1000);
    printf("Consensus time: %.2lf sec\n", cons_used_time / 1000);
    printf("Decoding time: %.2lf sec\n", dec_used_time / 1000);

    free(n_corrupted_ldpc);
    free(n_corrupted_rs);
    free(rec_payload);
    free(consensus_codewords);
    free(lratio_codewords_permutate);
    for (int i = 0; i < RS_N_SHORTEN; ++i)
    {
        free(*(lratio_codewords + i));
        free(*(product_codeword + i));
        free(*(rec_codeword_errfree + i));
    }
    free(lratio_codewords);
    free(product_codeword);
    free(rec_codeword_errfree);

    for (int i = 0; i < n_threads_file; ++i)
    {
        for (int j = 0; j < N_TYPE; ++j)
            free(*((corrCtx + i)->unique_pos + j));
    }

    thpool_destroy(thpool);
    pthread_mutex_destroy(&decoder_ctx.mutex);
    pthread_mutex_destroy(&cons_ctx.mutex);

    time(&end_time);
    strftime(timebuff, 100, "%Y-%m-%d %H:%M:%S ", localtime(&end_time));
    printf("\nEnd running: %s\n", timebuff);
    printf("Elapsed time: %.2lf sec\n\n", duration_run);
    return 0;
}

void pipelineWorker_chip1(void *args)
{
    WorkerHandler *worker_hdlr = (WorkerHandler *)args;

    pthread_mutex_lock(&worker_hdlr->mutex);
    int tid = worker_hdlr->tid;
    ++worker_hdlr->tid;

    int chip1_cnt = worker_hdlr->chip1_cnt_batch;
    ++worker_hdlr->chip1_cnt_batch;

    pthread_mutex_unlock(&worker_hdlr->mutex);

    char reads_file[500] = {0};
    sprintf(reads_file, "%s%d.fastq", worker_hdlr->chip1_reads_file_prefix, chip1_cnt);

    correlation_alignment(worker_hdlr, tid, reads_file);
}

void pipelineWorker_chip2(void *args)
{
    WorkerHandler *worker_hdlr = (WorkerHandler *)args;

    pthread_mutex_lock(&worker_hdlr->mutex);
    int tid = worker_hdlr->tid;
    ++worker_hdlr->tid;

    int chip2_cnt = worker_hdlr->chip2_cnt_batch;
    ++worker_hdlr->chip2_cnt_batch;

    pthread_mutex_unlock(&worker_hdlr->mutex);

    char reads_file[500] = {0};
    sprintf(reads_file, "%s%d.fastq", worker_hdlr->chip2_reads_file_prefix, chip2_cnt);

    correlation_alignment(worker_hdlr, tid, reads_file);
}

void correlation_alignment(WorkerHandler *worker_hdlr, int tid, const char *reads_file)
{
    unsigned data_size = worker_hdlr->data_size;
    double threshold_primer = worker_hdlr->threshold_primer;
    double threshold_autoCorr = worker_hdlr->threshold_autoCorr;
    bool mode = worker_hdlr->mode;

    int n_component = worker_hdlr->config->n_component;      // number of component code
    uint64_t period = worker_hdlr->config->period;           // period of composite code
    uint64_t start_phase = worker_hdlr->config->start_phase; // start phase of composite code
    unsigned *cLen = worker_hdlr->config->cLen;
    int **c = worker_hdlr->config->c;
    uint64_t *M = worker_hdlr->config->M;
    uint64_t *M_sup = worker_hdlr->config->M_sup;

    char *read_base = (char *)alloc_init(1000000, sizeof(char)); // A, T, G, C

    int read_len;
    int payload_len, payload_polished_len;

    int forward_flag = 1; // 1: forward, 0: reverse
    int primer_valid_flag = 0;
    int start_loc, end_loc;
    int minED_l_op, minED_r_op;
    int minED_l[2], minED_r[2];
    EdlibEqualityPair additionalEqualities[4] = {{'N', 'A'}, {'N', 'T'}, {'N', 'G'}, {'N', 'C'}};

    int **crossCorrVal = (int **)malloc(n_component * sizeof(int *));
    for (int i = 0; i < n_component; ++i)
        *(crossCorrVal + i) = (int *)alloc_init(*(cLen + i), sizeof(int));

    int recCorrVal, idealCorrVal;
    uint64_t recLoc;
    int cPhase_rec[MAX_COMPONENTS] = {0}; 

    int valid_flag = 0; 
    int64_t recIndex;  

    char payload_noisy_base[MAX_PAYLOAD_LEN] = {0};
    char payload_noisy[MAX_PAYLOAD_LEN] = {0};      
    char compositeCode_noisy[MAX_PAYLOAD_LEN] = {0};

    char compositeCode_regen[MAX_PAYLOAD_LEN] = {0};
    char payload_polished[MAX_PAYLOAD_LEN] = {0};    

    char compositeCode_AC[MAX_PAYLOAD_LEN] = {0};
    char quality[MAX_PAYLOAD_LEN] = {0};
    memset(quality, 'I', MAX_PAYLOAD_LEN);

    unsigned short **unique_pos = (worker_hdlr->corrCtx + tid)->unique_pos;

    unsigned n_reads = 0;
    unsigned n_reads_primer_valid = 0;
    unsigned n_reads_valid_corr = 0;
    unsigned n_reads_valid_align = 0;
    unsigned n_reads_BWA = 0;

    unsigned tot_read_base = 0;
    unsigned tot_payload_base = 0;

    time_t pre_primer_time, post_primer_time;
    time_t pre_corr_time, post_corr_time;
    time_t pre_align_time, post_align_time;
    double primer_used_time = 0.0;
    double corr_used_time = 0.0;
    double align_used_time = 0.0;

    FILE *fp_data = NULL;
    char ofile_data[500] = {0};
    FILE *fpw_PNseq = NULL;
    char ofile_PNseq[500] = {0};
    char filepath_SAM[500] = {0};
    if (worker_hdlr->mode)
    {
        sprintf(ofile_data, "%s/%s%03d.txt", worker_hdlr->log_dir, tmp_ofile_data_BWA, tid);
        if ((fp_data = fopen(ofile_data, "w+")) == NULL) 
        {
            fprintf(stderr, "Could not open file %s!\n", ofile_data);
        }

        sprintf(ofile_PNseq, "%s/%s%03d.fq", worker_hdlr->log_dir, tmp_ofile_PNseq_BWA, tid);
        if ((fpw_PNseq = fopen(ofile_PNseq, "w")) == NULL)
        {
            fprintf(stderr, "Could not open file %s!\n", ofile_PNseq);
        }

        sprintf(filepath_SAM, "%s/%s%03d.sam", worker_hdlr->log_dir, tmp_ofile_SAM_BWA, tid);
    }

    FILE *fpr_read = NULL;
    if ((fpr_read = fopen(reads_file, "r")) == NULL)
        fprintf(stderr, "Could not open file %s!\n", reads_file);

    while (fscanf(fpr_read, "%*[^\n]%*c") != EOF)
    {
        fscanf(fpr_read, "%s\n", read_base); // reads
        fscanf(fpr_read, "%*[^\n]%*c");
        fscanf(fpr_read, "%*[^\n]%*c");
        read_len = strlen(read_base);

        if (read_len > MAX_READ_LEN)
            read_len = MAX_READ_LEN;

        tot_read_base += read_len;
        ++n_reads;

        // Step.1 Paired-end Primers alignment & Determine boundary
        // --------------------------------------------------------
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

            start_loc = 0;
            if (minED_r_op <= threshold_primer * PRIMER_LEN)
            {
                primer_valid_flag = 1;

                if (end_loc >= 160)
                    start_loc = end_loc - 160 + 1;
            }
            for (int i = 0; i < 2; ++i)
                edlibFreeAlignResult(align_rslt_right[i]);
        }
        for (int i = 0; i < 2; ++i)
            edlibFreeAlignResult(align_rslt_left[i]);

        post_primer_time = current_time_ms();
        primer_used_time += difftime(post_primer_time, pre_primer_time);

        if (primer_valid_flag)
        {
            payload_len = end_loc - start_loc + 1;
            if (payload_len <= 0)
                continue;
            tot_payload_base += payload_len;

            ++n_reads_primer_valid;

            memset(payload_noisy_base, 0, MAX_PAYLOAD_LEN);
            memset(payload_noisy, 0, MAX_PAYLOAD_LEN);
            memset(compositeCode_noisy, 0, MAX_PAYLOAD_LEN);
            memset(compositeCode_AC, 0, MAX_PAYLOAD_LEN);
            for (int ni = 0; ni < payload_len; ++ni)
                *(payload_noisy_base + ni) = *(read_base + start_loc + ni);
            if (!forward_flag)
                reverseComplement(payload_noisy_base, payload_len);

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

            if (threshold_autoCorr <= 1.0)
            {
                pre_corr_time = current_time_ms();

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

                post_corr_time = current_time_ms();
                corr_used_time += difftime(post_corr_time, pre_corr_time);
            }
            else
            {
                recCorrVal = 0;
                idealCorrVal = 1;
            }

            valid_flag = 0;
            payload_polished_len = payload_len;
            for (int ni = 0; ni < payload_len; ++ni)
                *(payload_polished + ni) = *(payload_noisy + ni);
            *(payload_polished + payload_len) = '\0';
            if (1.0 * recCorrVal >= threshold_autoCorr * idealCorrVal)
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

            if (recIndex >= data_size || recIndex < 0)
                valid_flag = 0;
            if (worker_hdlr->mode && valid_flag == 0)
            {
                fprintf(fp_data, "%s\n", payload_noisy);


                fprintf(fpw_PNseq, "@%d%u\n", tid, n_reads_BWA);
                fprintf(fpw_PNseq, "%s\n", compositeCode_AC);
                fprintf(fpw_PNseq, "+\n");
                fwrite(quality, 1, payload_len, fpw_PNseq);
                fprintf(fpw_PNseq, "\n");
                ++n_reads_BWA;
            }

            if (valid_flag == 1)
            {
                ++n_reads_valid_corr;

                for (int ti = 0; ti < payload_polished_len; ++ti)
                {
                    if (recIndex + ti >= data_size || recIndex + ti < 0 || *(payload_polished + ti) == 'E')
                        continue;

                    ++*(*(unique_pos + (*(payload_polished + ti) - 48)) + recIndex + ti);
                }
            }
        }
    }
    fclose(fpr_read);

    int status;
    char cmdLine[1000] = {0};

    int flag;
    char tmp_symbol;
    int op_num, n_bp;
    int64_t startPos;
    char header[500] = {0};
    char cigar[1000] = {0};

    FILE *fpr_sam = NULL;

    // ================================= BWA workflow ==================================
    if (mode)
    {
        fclose(fpw_PNseq);

        if (n_reads_BWA > 0)
        {
            pre_align_time = current_time_ms();

            sprintf(cmdLine, "bwa mem -v 1 -t %d %s %s -o %s > /dev/null 2>&1", 
            worker_hdlr->n_threads_bwa, refer_file_PNseq, ofile_PNseq, filepath_SAM);
            status = system(cmdLine);
            if (status == -1)
            {
                fprintf(stderr, "Errors occurred while BWA mem!\n");
            }

            rewind(fp_data);

            if ((fpr_sam = fopen(filepath_SAM, "r")) == NULL)
            {
                fprintf(stderr, "Could not open file %s!\n", filepath_SAM);
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
                            fscanf(fp_data, "%*[^\n]%*c");
                        }
                        else if (flag == 0 || flag == 16)
                        {
                            ++n_reads_valid_align;

                            fscanf(fpr_sam, "%*s %ld %*s %s %*[^\n]s", &startPos, cigar);
                            startPos -= 1;
                            startPos -= start_phase;
                            fscanf(fp_data, "%s\n", payload_noisy);
                            payload_len = strlen(payload_noisy);

                            n_bp = 0;
                            op_num = 0;
                            payload_polished_len = 0;
                            memset(payload_polished, 0, MAX_PAYLOAD_LEN);
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
                                    op_num = 0;
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
            fclose(fpr_sam);

            post_align_time = current_time_ms();
            align_used_time += difftime(post_align_time, pre_align_time);
        }
        fclose(fp_data);
    }


    (worker_hdlr->corrCtx + tid)->n_reads = n_reads;
    (worker_hdlr->corrCtx + tid)->n_read_base = tot_read_base;
    (worker_hdlr->corrCtx + tid)->n_payload_base = tot_payload_base;
    (worker_hdlr->corrCtx + tid)->n_reads_primer_valid = n_reads_primer_valid;
    (worker_hdlr->corrCtx + tid)->n_reads_valid_corr = n_reads_valid_corr;

    (worker_hdlr->corrCtx + tid)->n_reads_BWA = n_reads_BWA;
    (worker_hdlr->corrCtx + tid)->n_reads_valid_align = n_reads_valid_align;

    (worker_hdlr->corrCtx + tid)->primer_used_time = primer_used_time;
    (worker_hdlr->corrCtx + tid)->corr_used_time = corr_used_time;
    (worker_hdlr->corrCtx + tid)->align_used_time = align_used_time;

    for (int j = 0; j < n_component; ++j)
        free(*(crossCorrVal + j));
    free(crossCorrVal);
}

void ConsensusWorker(void *args)
{
    ConsensuContext *cons_ctx = (ConsensuContext *)args;
    pthread_mutex_lock(&cons_ctx->mutex);
    int tid = cons_ctx->tid;
    ++cons_ctx->tid;
    pthread_mutex_unlock(&cons_ctx->mutex);

    unsigned bits_offset = cons_ctx->process_bits * tid;
    unsigned process_bits = (cons_ctx->n_jobs == tid + 1 && cons_ctx->n_jobs != 1) ? cons_ctx->data_size - bits_offset
                                                                                   : cons_ctx->process_bits;

    if (process_bits <= 0)
        return;

    int cof0, cof1;
    int confidence;

    CorrelationContext *corrCtx = cons_ctx->corrCtx;

    double *lratio_codewords = cons_ctx->lratio_codewords + bits_offset;
    char *consensus_codewords = cons_ctx->consensus_codewords + bits_offset;

    unsigned avail_base_num = 0;
    unsigned n_era_consensus = 0;

    for (unsigned sb = 0; sb < process_bits; ++sb)
    {

        cof0 = cof1 = 0;
        for (int i = 0; i < cons_ctx->n_tid_file; ++i)
        {
            cof0 += (int)*(*((corrCtx + i)->unique_pos + 0) + bits_offset + sb);
            cof1 += (int)*(*((corrCtx + i)->unique_pos + 1) + bits_offset + sb);
        }

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
    int tid = decoder_ctx->tid_ldpc;
    ++decoder_ctx->tid_ldpc;
    pthread_mutex_unlock(&decoder_ctx->mutex);

    unsigned code_offset = decoder_ctx->n_ldpc_per_tid * tid;
    int process_ldpc_num = (decoder_ctx->n_jobs_ldpc == tid + 1 && decoder_ctx->n_jobs_ldpc != 1) ? decoder_ctx->tot_ldpc - code_offset
                                                                                                  : decoder_ctx->n_ldpc_per_tid;
    if (process_ldpc_num <= 0)
        return;

    char *codeword_status = decoder_ctx->codeword_status + code_offset;

    double **lratio_codewords = decoder_ctx->lratio_codewords + code_offset;

    char **bit_status = decoder_ctx->bit_status + code_offset;
    char **rec_codeword_errfree = decoder_ctx->rec_codeword_errfree + code_offset;

    char **product_codeword = decoder_ctx->product_codeword + code_offset;

    LDPCDecode *ldpc_decode = new LDPCDecode();
    if (ldpc_decode->decode_init(ldpc_pchk_file))
    {
        fprintf(stderr, "Failed to initialize LDPC decoder!\n");
        return;
    }

    int hard_bit;

    int ldpc_chk;
    int n_corrupted_ldpc = 0;
    for (int blk = 0; blk < process_ldpc_num; ++blk)
    {
        if (*(codeword_status + blk) == 0)
        {
            // Replace with error-free bits
            for (int ni = 0; ni < LDPC_K; ++ni)
            {
                if (*(*(bit_status + blk) + ni))
                {
                    hard_bit = *(*(rec_codeword_errfree + blk) + ni);
                    *(*(lratio_codewords + blk) + LDPC_M + ni) = exp(-6.3 * (1 - 2 * hard_bit));
                }
            }

            //---- LDPC decoding
            ldpc_chk = ldpc_decode->decode_LDPC_soft(*(product_codeword + blk), *(lratio_codewords + blk), maxIteration_ldpc);

            if (ldpc_chk)
                ++n_corrupted_ldpc;
            else
            {
                *(codeword_status + blk) = 1;
                for (int ni = 0; ni < LDPC_K; ++ni)
                {
                    *(*(bit_status + blk) + ni) = 1;
                    *(*(rec_codeword_errfree + blk) + ni) = *(*(product_codeword + blk) + ni);
                }
            }
        }
        else
        {
            for (int ni = 0; ni < LDPC_K; ++ni)
                *(*(product_codeword + blk) + ni) = *(*(rec_codeword_errfree + blk) + ni);
        }
    }
    *(decoder_ctx->n_corrupted_ldpc + tid) = n_corrupted_ldpc;

    fprintf(stderr, "Finished decoding of %d LDPC codes! Corrupted LDPC: %d\n", process_ldpc_num, n_corrupted_ldpc);

    ldpc_decode->mem_free_dec();
    delete ldpc_decode;
    ldpc_decode = NULL;
}


void RSDecodingWorker(void *args)
{
    DecoderContext *decoder_ctx = (DecoderContext *)args;

    pthread_mutex_lock(&decoder_ctx->mutex);
    int tid = decoder_ctx->tid_rs;
    ++decoder_ctx->tid_rs;
    pthread_mutex_unlock(&decoder_ctx->mutex);

    unsigned code_offset = decoder_ctx->n_rs_per_tid * tid;
    int process_rs_num = (decoder_ctx->n_jobs_rs == tid + 1 && decoder_ctx->n_jobs_rs != 1) ? decoder_ctx->tot_rs - code_offset
                                                                                            : decoder_ctx->n_rs_per_tid;
    if (process_rs_num <= 0)
        return;

    char **bit_status = decoder_ctx->bit_status;
    char **rec_codeword_errfree = decoder_ctx->rec_codeword_errfree;

    char **product_codeword = decoder_ctx->product_codeword;

    int zero_padding = RS_K - RS_K_SHORTEN;
    int rs_m = RS_N - RS_K;
    int num_GF;
    int rs_cw[RS_N] = {0};

    RSCodec *rs_decoder = new RSCodec(RS_ORDER, RS_N, RS_K);

    char success_flag = 0;
    char cur_bit;

    int n_corrupted_pre_rs = 0;
    int n_corrupted_post_rs = 0;
    for (int col = 0; col < process_rs_num; ++col)
    {
        success_flag = 0;
        memset(rs_cw, 0, zero_padding * sizeof(int));

        for (int row = 0; row < RS_N_SHORTEN; ++row)
        {
            num_GF = 0;
            for (int j = 0; j < RS_ORDER; ++j)
                num_GF += *(*(product_codeword + row) + (code_offset + col) * RS_ORDER + j) * (1 << (RS_ORDER - 1 - j));
            *(rs_cw + zero_padding + row) = num_GF;
        }

        if (rs_decoder->code_check(rs_cw))
            n_corrupted_pre_rs++;

        rs_decoder->decode_rs(rs_cw, 0, NULL);

        if (rs_decoder->code_check())
            n_corrupted_post_rs++;
        else
            success_flag = 1;

        for (int row = 0; row < RS_K_SHORTEN; ++row)
        {
            for (int j = 0; j < RS_ORDER; ++j)
            {
                cur_bit = (char)((*(rs_cw + zero_padding + row) >> j) % 2);
                *(*(product_codeword + row) + (code_offset + col) * RS_ORDER + RS_ORDER - 1 - j) = cur_bit;

                if (success_flag && *(*(bit_status + row) + (code_offset + col) * RS_ORDER + RS_ORDER - 1 - j) == 0)
                {
                    *(*(bit_status + row) + (code_offset + col) * RS_ORDER + RS_ORDER - 1 - j) = 1;
                    *(*(rec_codeword_errfree + row) + (code_offset + col) * RS_ORDER + RS_ORDER - 1 - j) = cur_bit;
                }
            }
        }

        if (success_flag)
        {
            for (int row = 0; row < rs_m; ++row)
            {
                for (int j = 0; j < RS_ORDER; ++j)
                {
                    if (*(*(bit_status + RS_K_SHORTEN + row) + (code_offset + col) * RS_ORDER + RS_ORDER - 1 - j) == 0)
                    {
                        cur_bit = (char)((*(rs_cw + RS_K + row) >> j) % 2);
                        *(*(bit_status + RS_K_SHORTEN + row) + (code_offset + col) * RS_ORDER + RS_ORDER - 1 - j) = 1;
                        *(*(rec_codeword_errfree + RS_K_SHORTEN + row) + (code_offset + col) * RS_ORDER + RS_ORDER - 1 - j) = cur_bit;
                    }
                }
            }
        }
    }
    *(decoder_ctx->n_corrupted_rs + tid) = n_corrupted_post_rs;

    fprintf(stderr, "Finished decoding of %d RS codes! Pre-corrupted: %d, post-corrupted: %3d\n",
            process_rs_num, n_corrupted_pre_rs, n_corrupted_post_rs);

    delete rs_decoder;
    rs_decoder = NULL;
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

    status = pclose(pipe);

    if (lines == set_lines * 4)
        return 1;
    else
        return 0;
}