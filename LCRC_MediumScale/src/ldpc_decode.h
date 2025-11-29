#ifndef _LDPCDECODE_H
#define _LDPCDECODE_H

#include "mod2sparse.h"
#include "mod2dense.h"

#define ROW 30 * 360     // 10800
#define COLUMN 180 * 360 // 64800
#define CHECK_DEGREE 20

#define Prob_ALGORITHM

class LDPCDecode
{
public:
    explicit LDPCDecode();

    int decode_init(const char *pchk_file);
    void mem_free_dec();

    int decode_LDPC_soft(char *decMsg, double *lratio_cw, int maxIteration = 200);

    int code_check(char *dblk, char *cblk);

private:
    /* VARIABLES DECLARED .  These global variables are set to
       representations of the parity check and generator matrices by read_pchk
       and read_gen. */
    mod2sparse *H; /* Parity check matrix */

    int M; /* Number of rows in parity check matrix */
    int N; /* Number of columns in parity check matrix */

    char type; /* Type of generator matrix representation (s/d/m) */
    int *cols; /* Ordering of columns in generator matrix */

    mod2sparse *L, *U; /* Sparse LU decomposition, if type=='s' */
    int *rows;         /* Ordering of rows in generator matrix (type 's') */

    mod2dense *G; /* Dense or mixed representation of generator matrix,
             if type=='d' or type=='m' */

    /* READ PARITY CHECK MATRIX.  Sets the H, M, and N global variables.  If an
        error is encountered, a message is displayed on standard error, and the
        program is terminated. */
    int read_pchk(char *);
    int read_gen(char *, int, int);

    /* DECODING METHOD, ITS PARAMETERS, AND OTHER VARIABLES. */
    void initprp(mod2sparse *, double *, char *, double *);
    void iterprp(mod2sparse *, double *, char *, double *);
};

#endif