/************************************
 * Copyright (C) 2025 Tianjin University (TJU)
 *
 * @brief: Error correction for oligo pool data storage
 *
 *************************************/
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cmath>
#include <cassert>

#include "pthread.h"
#include "ldpc_decode.h"
#include "utils.h"

typedef struct
{
    int tid;
    int code_num;
    int code_offset;

    char *mSeq1;
    char *mSeq2;

    unsigned *permutation;
    double *lratio_codewords;

    char **rec_payload;

} DecoderContext;

#define IMG1_BYTE 100161
#define IMG2_BYTE 95558

#define BLOCKCODE_NUM 29

#define LDPC_N 64800
#define LDPC_M 10800
#define LDPC_K 54000

#define THREAD_MAX 64

#define MSEQ1 65535
#define MSEQ2 31

const char *ldpc_pchk_file = "config/ldpc.pchk";

char bit2ascii(const char *inBit)
{
    char symbol = 0;
    for (int i = 0; i < 8; i++)
    {
        symbol += *(inBit + i) * (1 << (7 - i));
    }
    return symbol;
}
void *LDPCDecodingWorker(void *args);

int main(int argc, char *argv[])
{
    int n_thread = 1;
    char result_dir[500] = {0};
    char consensus_file[500] = {0};

    if (argc != 4)
    {
        fprintf(stderr, "Usage: "
                        "%s [Result dir] [consensus file] [thread]\n",
                argv[0]);
        return 1;
    }
    sscanf(argv[1], "%s", result_dir);
    sscanf(argv[2], "%s", consensus_file);
    sscanf(argv[3], "%d", &n_thread);

    time_t start_time, end_time;
    char timebuff[100];

    printf("-----Step 3: Error correction and data recovery----\n\n");

    time(&start_time);
    strftime(timebuff, 100, "%Y-%m-%d %H:%M:%S ", localtime(&start_time));
    printf("Start running: %s\n\n", timebuff);

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

    //----- Load consensus codeword
    double *lratio_codewords = (double *)alloc_init(BLOCKCODE_NUM * LDPC_N, sizeof(double));
    FILE *fpr_cons = NULL;
    if ((fpr_cons = fopen(consensus_file, "rb")) == NULL)
    {
        fprintf(stderr, "Could not open file '%s'!\n", consensus_file);
        return 1;
    }
    size_t rt = fread(lratio_codewords, sizeof(double), BLOCKCODE_NUM * LDPC_N, fpr_cons);
    assert(rt == BLOCKCODE_NUM * LDPC_N);
    fclose(fpr_cons);

    unsigned *randVal = (unsigned *)alloc_init(65535, sizeof(unsigned));
    unsigned *permutation = (unsigned *)alloc_init(LDPC_N, sizeof(unsigned));
    PRNG(randVal, 16, 65535, 64);
    for (unsigned i = 0, nn = 0; i < 65535; ++i)
    {
        if (*(randVal + i) <= LDPC_N)
        {
            *(permutation + nn) = *(randVal + i) - 1;
            ++nn;
        }
    }

    char **rec_payload = (char **)malloc(BLOCKCODE_NUM * sizeof(char *));
    for (int i = 0; i < BLOCKCODE_NUM; ++i)
        *(rec_payload + i) = (char *)alloc_init(LDPC_K, sizeof(char));

    int n_thread_on = 0;
    pthread_t tid[THREAD_MAX];
    DecoderContext decoder_ctx[THREAD_MAX];
    int divider = ceil(1.0 * BLOCKCODE_NUM / n_thread);
    n_thread_on = n_thread * divider > BLOCKCODE_NUM ? ceil(1.0 * BLOCKCODE_NUM / divider) : n_thread;

    for (int ni = 0; ni < n_thread_on; ++ni)
    {
        (decoder_ctx + ni)->tid = ni;
        (decoder_ctx + ni)->code_num = (ni == n_thread_on - 1) ? BLOCKCODE_NUM - divider * ni : divider;
        (decoder_ctx + ni)->code_offset = divider * ni;

        (decoder_ctx + ni)->mSeq1 = mSeq1;
        (decoder_ctx + ni)->mSeq2 = mSeq2;
        (decoder_ctx + ni)->lratio_codewords = lratio_codewords + LDPC_N * divider * ni;
        (decoder_ctx + ni)->permutation = permutation;

        (decoder_ctx + ni)->rec_payload = rec_payload + divider * ni;

        if (pthread_create(tid + ni, NULL, &LDPCDecodingWorker, (void *)(decoder_ctx + ni)) != 0)
        {
            fprintf(stderr, "Failed to create thread-%d!\n", ni);
            exit(1);
        }
    }
    for (int ni = 0; ni < n_thread_on; ++ni)
    {
        if (pthread_join(*(tid + ni), NULL) != 0)
        {
            fprintf(stderr, "Thread-%d failed!\n", ni);
        }
    }

    //---- Data recovery ----//
    size_t tot_bytes = BLOCKCODE_NUM * LDPC_K / 8;
    char *pbuff_ptr = NULL;
    char *rec_data_bytes = (char *)alloc_init(tot_bytes, sizeof(char));

    int n_byte = 0;
    for (int blk = 0; blk < BLOCKCODE_NUM; ++blk)
    {
        for (int ni = 0; ni < LDPC_K / 8; ++ni)
        {
            pbuff_ptr = *(rec_payload + blk) + ni * 8;
            *(rec_data_bytes + n_byte) = bit2ascii(pbuff_ptr);
            n_byte++;
        }
    }

    FILE *fpw = NULL;
    char img1_filename[500] = {0};
    char img2_filename[500] = {0};

    sprintf(img1_filename, "%s/img1.jpg", result_dir);
    sprintf(img2_filename, "%s/img2.jpg", result_dir);

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

    free(lratio_codewords);
    free(randVal);
    free(permutation);
    for (int i = 0; i < BLOCKCODE_NUM; ++i)
        free(*(rec_payload + i));
    free(rec_payload);

    time(&end_time);
    strftime(timebuff, 100, "%Y-%m-%d %H:%M:%S ", localtime(&end_time));
    printf("\nEnd running: %s\n", timebuff);
    printf("Elapsed time: %.2lf sec\n\n", difftime(end_time, start_time));
    return 0;
}

void *LDPCDecodingWorker(void *args)
{
    DecoderContext *decoder_ctx = (DecoderContext *)args;

    int code_offset = decoder_ctx->code_offset;
    char *mSeq1 = decoder_ctx->mSeq1;
    char *mSeq2 = decoder_ctx->mSeq2;

    double *lratio_codewords = decoder_ctx->lratio_codewords;
    unsigned *permutation = decoder_ctx->permutation;
    char **rec_payload = decoder_ctx->rec_payload;

    double *lratio = (double *)alloc_init(LDPC_N, sizeof(double));

    LDPCDecode *ldpc_decode = new LDPCDecode();
    if (ldpc_decode->decode_init(ldpc_pchk_file))
    {
        fprintf(stderr, "Failed to initialize LDPC decoder!\n");
        return NULL;
    }

    int ldpc_chk;
    for (int blk = 0; blk < decoder_ctx->code_num; ++blk)
    {
        //// de-permuation
        for (int nb = 0; nb < LDPC_N; ++nb)
        {
            *(lratio + *(permutation + nb)) = *(lratio_codewords + LDPC_N * blk + nb);
        }

        //---- LDPC(64800, 54000) decoding
        ldpc_chk = ldpc_decode->decode_LDPC_soft(*(rec_payload + blk), lratio, 200);

        //---- descramble
        for (int nb = 0; nb < LDPC_K; ++nb)
        {
            *(*(rec_payload + blk) + nb) = (*(*(rec_payload + blk) + nb) +
                                            *(mSeq1 + ((unsigned)LDPC_K * (code_offset + blk) + nb) % MSEQ1) +
                                            *(mSeq2 + ((unsigned)LDPC_K * (code_offset + blk) + nb) % MSEQ2)) %
                                           2;
        }
    }

    free(lratio);
    delete ldpc_decode;
}
