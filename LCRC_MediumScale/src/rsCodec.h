#ifndef _RSCODEC_H
#define _RSCODEC_H
#include <cstdlib>

class RSCodec
{
private:

	int m;				 
	int n, k;			  
	int check_s, t, d_min;

	int init_zero;
	int *g, *p, *index_of, *alpha_to;

	int *sdata, *b;

	int *rdata;

	void read_p();
	void generate_gf();
	void gen_poly();

public:
	explicit RSCodec(int _order, int _cwLen, int _infoLen, int _initZero = 0);
	~RSCodec();
	void encode_rs(int *parity, const int *msg_nb);
	void decode_rs(int *rvd_cw, const int numera = 0, const int *era = NULL);

	int code_check();
	int code_check(int *codeword);
};

#endif // !_RS_CODEC_
