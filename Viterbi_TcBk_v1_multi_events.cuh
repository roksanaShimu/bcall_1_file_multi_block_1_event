/*
 * Viterbi_TcBk_v1_multi_events.cuh
 *
 *  Created on: Jun 20, 2017
 *      Author: roksana
 */

#ifndef VITERBI_TCBK_V1_MULTI_EVENTS_CUH_
#define VITERBI_TCBK_V1_MULTI_EVENTS_CUH_






#include<iostream>
#include<cuda.h>
#include <ctime>
#include <math.h>
#include <cuda_runtime.h>
#include <vector>
#include <thrust/device_vector.h>

#define __BCALL_I__
#include "bcall-i.h"
#include "bcall.h"
#include "Kmer.cuh"
#include "Predefine_Values.cuh"
#include "Device_Memory_Allocation.cuh"

using namespace std;



#define Q 18

/* Conversion form float to fixed point */
#define QX(F) ((int32_t)((F) * (1 << (Q))))

/* Convert from fixed point to float point */
#define QXF(Q) ((float)(((float)(Q)) / (1.0 * (1 << (Q)))))

/* Multiply */
#define QXMUL(X1, X2) (int32_t)(((int64_t)(X1) * (int64_t)(X2)) >> (Q))

/* MinION event definition */
struct fix_event {
	int32_t mean;
	int32_t stdv;
	int32_t stdv_inv;
	int32_t log_stdv_1p5;
};

/* Pore model state */
struct fix_pm_state {
	int32_t level_mean;
	int32_t sd_mean;
	int32_t sd_mean_inv;
	int32_t level_stdv_inv_2;
	int32_t sd_lambda_p5;
	int32_t log_level_stdv_2pi;
	int32_t log_sd_lambda_p5;
};

/* Unroll the inner loop */
#define UNROLL_INNER 1

/* Transitions map */
struct fix_map {
	int32_t log_pr[MAX_TRANS];
	struct fix_pm_state pm;
};

/* Length of the trace back buffer */
#define TRACEBACK_LEN 128
/* To disable continuous traceback set TRACEBACK_LEN to N_EVENTS_MAX */
//#define TRACEBACK_LEN N_EVENTS_MAX

/* Period (in events) in which the fixpt_traceback() function is called */
#define TRACEBACK_CNT 32
/* To disable continuous traceback set TRACEBACK_CNT to N_EVENTS_MAX */
//#define TRACEBACK_CNT N_EVENTS_MAX

struct bcall_tcbk {
	/* Transitions map */
	struct fix_map map[N_STATES];
	/* previous state in the MLSS */
	uint16_t beta[TRACEBACK_LEN][N_STATES];
	/* Store the maximum probability for the most probable path */
	int64_t sum_log_p_max;
	/* Sequence of events */
	struct fix_event event[N_EVENTS_MAX];
	unsigned int n_events;
	uint16_t state_seq[N_EVENTS_MAX];
	/* Resulting sequence of bases (output) */
	char base_seq[N_EVENTS_MAX];
	unsigned int base_cnt;
};

/* ---------------------------------------------------------------------------
 * Viterbi algorithm
 * ---------------------------------------------------------------------------
 */


__constant__ int32_t d_ev[4];

static inline int32_t ln_pemission(struct fix_pm_state * pm,
							struct fix_event * e,
							int32_t p_max,  int32_t ln_pt)
{
	int32_t pm_mu = pm->sd_mean;
	int32_t pm_r_mu = pm->sd_mean_inv;
	int32_t pm_lmbd1 = pm->sd_lambda_p5;
	int32_t pm_mean = pm->level_mean;
	int32_t pm_r_stdv1 = pm->level_stdv_inv_2;
	int32_t pm_l_lmbd1 = pm->log_sd_lambda_p5;
	int32_t pm_l_stdv1 = pm->log_level_stdv_2pi;

	int32_t e_stdv = e->stdv;
	int32_t e_r_stdv = e->stdv_inv;
	int32_t e_mean = e->mean;
	int32_t e_l_stdv1 = e->log_stdv_1p5;

	/* First stage */
	int32_t x1;
	int32_t x2;
	int32_t x3;
	int32_t x4;
	int32_t x5;
	int32_t x6;
	int32_t x7;
	int32_t x8;
	/* Second stage */
	int32_t y1;
	int32_t y2;
	int32_t y3;
	int32_t y4;
	int32_t y5;
	/* Third stage */
	int32_t z1;
	int32_t z2;
	int32_t z3;
	int32_t z4;
	/* Forth stage */
	int32_t w1;
	int32_t w2;
	/* Fifth stage */
	int32_t v;

	/* First stage */
	x1 = ln_pt + pm_l_lmbd1;
	x2 = e_l_stdv1 + p_max;
	x3 = pm_l_stdv1;
	x5 = pm_r_stdv1;
	x8 = QXMUL(pm_lmbd1, e_r_stdv);
	x4 = e_mean - pm_mean;
	x6 = e_stdv - pm_mu;
	x7 = pm_r_mu;

	/* Second stage */
	y1 = x1 - x2;
	y2 = x3;
	y3 = QXMUL(x4, x5);
	y4 = QXMUL(x6, x7);
	y5 = x8;

	/* Third stage */
	z1 = y1 - y2;
	z2 = QXMUL(y3, y3);
	z3 = QXMUL(y4, y4);
	z4 = y5;

	/* Forth stage */
	w1 = z1 - z2;
	w2 = QXMUL(z3, z4);

	/* Fifth stage */
	v = w1 - w2;

	return v;
}


static inline void fixpt_traceback(struct bcall_tcbk * fix,
								   int to, unsigned int j_max)
{
	unsigned int from;
	unsigned int i;

	if (to < TRACEBACK_LEN) {
		from = 0;
	} else {
		from = to - TRACEBACK_LEN;
	}

	//DBG(DBG_INFO, "from=%d to=%d", from, to);
	//cout<<"from="<<from<<"   to="<<to<<endl;

	for (i = to - 1; i > from; --i) {
		fix->state_seq[i] = j_max;
		j_max = fix->beta[i % TRACEBACK_LEN][j_max];
	}

	fix->state_seq[i] = j_max;
}


__device__ static inline unsigned int dev_kmer_step(unsigned int state, unsigned int base) {
	return (state >> 2) + (base << (KMER_SIZE * 2 - 2));
}

__device__ static inline unsigned int dev_kmer_skip(unsigned int state, unsigned int seq) {
	return (state >> 4) + (seq << (KMER_SIZE * 2 - 4));
}


__device__ int32_t dev_ln_ptransition(int32_t * devPtr_log_pr, unsigned j, int32_t *alpha, uint16_t *beta){

	unsigned int j1;
	int32_t log_p;
	int32_t y;
	uint16_t temp_beta;



	//start index of devPtr_log_pr
	unsigned group_no=j/warp_size;
	unsigned start_index=group_no*warp_size*21 +j%warp_size;

	//if(j<10)printf("inside dev_ln_ptransition: j=%d   devPtr_log_pr[start_index]=%d \n", j,devPtr_log_pr[start_index] );

	/* Stay */
	log_p = devPtr_log_pr[start_index] + alpha[j];
	temp_beta = j;

	/* Step */
	j1 = dev_kmer_step(j, 0);
	y = devPtr_log_pr[start_index+1*warp_size]+ alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_step(j, 1);
	y = devPtr_log_pr[start_index+2*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_step(j, 2);
	y = devPtr_log_pr[start_index+3*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_step(j, 3);
	y = devPtr_log_pr[start_index+4*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}


	/* Skip */
	j1 = dev_kmer_skip(j, 0);
	y = devPtr_log_pr[start_index+5*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 1);
	y = devPtr_log_pr[start_index+6*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 2);
	y = devPtr_log_pr[start_index+7*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 3);
	y = devPtr_log_pr[start_index+8*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 4);
	y = devPtr_log_pr[start_index+9*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 5);
	y = devPtr_log_pr[start_index+10*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 6);
	y = devPtr_log_pr[start_index+11*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 7);
	y = devPtr_log_pr[start_index+12*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 8);
	y = devPtr_log_pr[start_index+13*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 9);
	y = devPtr_log_pr[start_index+14*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 10);
	y = devPtr_log_pr[start_index+15*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 11);
	y = devPtr_log_pr[start_index+16*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 12);
	y = devPtr_log_pr[start_index+17*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 13);
	y = devPtr_log_pr[start_index+18*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 14);
	y = devPtr_log_pr[start_index+19*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}
	j1 = dev_kmer_skip(j, 15);
	y = devPtr_log_pr[start_index+20*warp_size] + alpha[j1];
	if (y > log_p) {
		log_p = y;
		temp_beta = j1;
	}


	/*if(j<0 || j>=N_STATES){
		printf("j= %d is the out of range\n", j) ;
	}else{
		printf("j is in the range\n");
	}*/

	beta[j]=temp_beta;
	return log_p;

}


__device__ static inline int32_t dev_ln_pemission(int32_t *devPtr_pm, int32_t p_max,  int32_t ln_pt, unsigned j)
{


	//start index for devPtr_pm
	unsigned group_no=j/warp_size;
	unsigned start_index=group_no*warp_size*7 +j%warp_size;


	/* First stage */
	int32_t x1 = ln_pt + devPtr_pm[start_index+ 5*warp_size];
	int32_t x2 = d_ev[3] + p_max;
	int32_t x8 = QXMUL(devPtr_pm[start_index+ 2*warp_size], d_ev[1]);
	int32_t x4 = d_ev[2] - devPtr_pm[start_index+ 3*warp_size];
	int32_t x6 = d_ev[0] - devPtr_pm[start_index];


	/* Second stage */
	 x1 = x1 - x2;

	x2 = QXMUL(x4, devPtr_pm[start_index+ 4*warp_size]);
	x4 = QXMUL(x6, devPtr_pm[start_index+ 1*warp_size]);


	/* Third stage */
	x1 = x1 - devPtr_pm[start_index+ 6*warp_size];
	x2 = QXMUL(x2, x2);
	x4 = QXMUL(x4, x4);


	/* Forth stage */
	x1 = x1 - x2;
	x2 = QXMUL(x4, x8);

	/* Fifth stage */
	x1 = x1 - x2;

	return x1;
	//return 0;
}



__global__ void kernel(int32_t *devPtr_pm, int32_t *devPtr_log_pr, int32_t *alpha, uint16_t *beta, int32_t p_max, int32_t *d_amax, unsigned int *d_jmax, int32_t *d_mutex, int ith_event){

	unsigned id = threadIdx.x + blockIdx.x * threads_per_block;   //id is j from the serial code

	unsigned group_block=id/N_STATES;
	unsigned start_block=group_block* number_of_block_per_group;
	unsigned end_block= start_block+number_of_block_per_group -1;


	__syncthreads();

	
	//TD DO: try to move alpha to shared memory so that it can save time accessing global memory
	int32_t ln_pt=dev_ln_ptransition(devPtr_log_pr, id, alpha, beta);

	int32_t ln_pe = dev_ln_pemission(devPtr_pm, p_max,  ln_pt, id);

	alpha[id] = ln_pe;
			
	if(ith_event==1) if(threadIdx.x==0) printf("id=%d and ln_pt=%d and ln_pe=%d\n",id,ln_pt, ln_pe);
		

	__syncthreads();



	//find max ln_pe and the corresponding id
	__shared__ int cache[threads_per_block];
	cache[threadIdx.x]=ln_pe;
	__syncthreads();


	// reduction
	unsigned int i = blockDim.x/2;
	while(i != 0){
		if(threadIdx.x < i){
			cache[threadIdx.x] = max(cache[threadIdx.x], cache[threadIdx.x + i]);
		}

		__syncthreads();
		i /= 2;
	}
	__syncthreads();


	
	if(blockIdx.x >= start_block && blockIdx.x <= end_block){

		if(threadIdx.x == 0){
			while(atomicCAS(d_mutex,0,1) != 0);  //lock
			*d_amax = max(*d_amax, cache[0]);
			atomicExch(d_mutex, 0);  //unlock
		}
	}


	//---------------------------------------------------------------


	__syncthreads();
	
	if(ln_pe== *d_amax){
		
		*d_jmax=id;
		
	}
	__syncthreads();
	


}




void preprocess_for_parallel(struct bcall_tcbk * fix, int32_t alpha[N_STATES], int32_t p_max, int i, int32_t * devPtr_pm, int32_t * devPtr_log_pr, int32_t * dev_alpha,  uint16_t *dev_beta, uint16_t *h_beta, int32_t *h_amax, int32_t *d_amax, unsigned int *h_jmax, unsigned int *d_jmax, int32_t *d_mutex ){

	struct fix_event * e;

	int j, k;


	/* XXX: Read from memory */
	e = &fix->event[i];

	//event
	int32_t event_host[4];

	event_host[0]=e->stdv;
	event_host[1]=e->stdv_inv;
	event_host[2]=e->mean;
	event_host[3]=e->log_stdv_1p5;


	/* adjust the maximum probability for this round */
	*h_amax = INT32_MIN; *h_jmax=N_STATES;


	// set up timing variables
	/*float gpu_elapsed_time;
	cudaEvent_t gpu_start, gpu_stop;
	cudaEventCreate(&gpu_start);
	cudaEventCreate(&gpu_stop);

	// copy from host to device
	cudaEventRecord(gpu_start, 0);*/

/*
	const int num_streams = 2;

	cudaStream_t streams[num_streams];
	float *data[num_streams];

*/

	cudaMemset(d_mutex, 0, sizeof(int32_t));


	//copy event_host to constant memory
	cudaMemcpyToSymbol(d_ev, event_host, 4*sizeof(int32_t));

	cudaMemcpy(d_amax, h_amax, sizeof(int32_t), cudaMemcpyHostToDevice);
	cudaMemcpy(d_jmax, h_jmax, sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_alpha,alpha, sizeof(int32_t)*N_STATES, cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();
	kernel<<< number_of_blocks, threads_per_block >>>(devPtr_pm, devPtr_log_pr, dev_alpha, dev_beta, p_max, d_amax, d_jmax, d_mutex, i);
	cudaDeviceSynchronize();

	cudaMemcpy(alpha,dev_alpha, sizeof(int32_t)*N_STATES, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_beta,dev_beta, sizeof(uint16_t)*N_STATES, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_amax,d_amax, sizeof(int32_t), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_jmax,d_jmax, sizeof(int32_t), cudaMemcpyDeviceToHost);


	/*cudaEventRecord(gpu_stop, 0);
	cudaEventSynchronize(gpu_stop);
	cudaEventElapsedTime(&gpu_elapsed_time, gpu_start, gpu_stop);
	cudaEventDestroy(gpu_start);
	cudaEventDestroy(gpu_stop);
	*/



	for (j = 0; j < N_STATES; ++j) {
		fix->beta[i % TRACEBACK_LEN][j]=h_beta[j];
	}



}



struct bcall_tcbk * tcbk_blk_alloc(void)
{
	struct bcall_tcbk * fix;

	fix = (struct bcall_tcbk *)malloc(sizeof(struct bcall_tcbk));

	return fix;
}



/* Apply drift correction */
void tcbk_events_init(struct bcall_tcbk * fix,
					  const struct event_entry entry[],
					  unsigned int len, double drift)
{
	unsigned int i;

	assert(len < N_EVENTS_MAX);

	for (i = 0; i < len; ++i) {
		double mean;
		double stdv;

		mean = entry[i].mean - drift * entry[i].start;
		stdv = entry[i].stdv;

		fix->event[i].mean = QX(mean);
		fix->event[i].stdv = QX(stdv);
		fix->event[i].stdv_inv = QX(1.0/stdv);
		fix->event[i].log_stdv_1p5 = QX(1.5 * log(stdv));
	}

	//DBG(DBG_INFO, "entries=%d", i);
	cout<<" entries= "<< i<<endl;
	fix->n_events = i;
}

int tcbk_strand_load(struct bcall_tcbk * fix, struct sequence * seq,
					 unsigned int st_no, unsigned int pm_no)
{
	struct pm_state pm_state[N_STATES];
	unsigned int i;

	if (st_no >= seq->st_cnt)
		return -1;

	if (pm_no >= seq->st[st_no].pm_cnt)
		return -1;

	//DBG(DBG_INFO, "\"%s\"", seq->st[st_no].pm[pm_no].name);
	cout<<seq->st[st_no].pm[pm_no].name<<endl;

	pore_model_init(pm_state, seq->st[st_no].pm[pm_no].pm_entry, N_STATES);
	pore_model_scale(pm_state, seq->st[st_no].pm[pm_no].pm_param);

	for (i = 0; i < N_STATES; ++i) {
		fix->map[i].pm.level_mean = QX(pm_state[i].level_mean);
		fix->map[i].pm.sd_mean = QX(pm_state[i].sd_mean);
		fix->map[i].pm.sd_mean_inv = QX(1.0/pm_state[i].sd_mean);
		fix->map[i].pm.level_stdv_inv_2 = QX(1.0/ (sqrt(2.0) * pm_state[i].level_stdv));
		fix->map[i].pm.sd_lambda_p5 = QX(pm_state[i].sd_lambda / 2);
		fix->map[i].pm.log_level_stdv_2pi = QX(pm_state[i].log_level_stdv +
											   LOG_2PI);
		fix->map[i].pm.log_sd_lambda_p5 = QX(pm_state[i].log_sd_lambda / 2);
	}

	tcbk_events_init(fix, seq->st[st_no].ev, seq->st[st_no].ev_cnt,
					 seq->st[st_no].pm[pm_no].pm_param->drift);

	return seq->st[st_no].ev_cnt;
}

static void transition_add(struct bcall_tcbk * fix, unsigned int n,
						   unsigned int from, unsigned int to,
						   double p_stay, double p_step,
						   double p_skip_1)
{
	double p = get_trans_prob(from, to, p_stay, p_step, p_skip_1);
	fix->map[to].log_pr[n] = QX(log(p));
}
void tcbk_compute_transitions(struct bcall_tcbk * fix, struct st_params * param)
{
	double p_stay = param->p_stay;
	double p_skip = param->p_skip;
	double p_step = 1.0 - p_stay - p_skip;
	double p_skip_1 = p_skip / (p_skip + 1.0);
	unsigned int i;

	for (i = 0; i < N_STATES; ++i) {
		unsigned int n = 0;
		unsigned int j1;
		unsigned int k;

		/* Stay */
		transition_add(fix, n++, i, i, p_stay, p_step, p_skip_1);
		/* Step */
		for (k = 0; k < 4; ++k) {
			j1 = kmer_step(i, k);
			transition_add(fix, n++, j1, i, p_stay, p_step, p_skip_1);
		}
		/* Skip */
		for (k = 0; k < 16; ++k) {
			j1 = kmer_skip(i, k);
			transition_add(fix, n++, j1, i, p_stay, p_step, p_skip_1);
		}
	}
}

void tcbk_sequence_fill(struct bcall_tcbk * fix)
{
	/* Pr[ MLSS producing e_1 ... e_i, with S_i == j ] */
	int32_t alpha[N_STATES];
	struct fix_event * e;
	unsigned int j_max;
	int32_t a_max;
	unsigned int i;
	unsigned int j, k;
	/* Store the maximum probability in each round */
	int32_t p_max = INT32_MIN;

	//DBG(DBG_INFO, "begin");
	cout<<"begin"<<endl;
	fix->sum_log_p_max = 0;

	/* XXX: Read from memory */
	e = &fix->event[0];
	/* adjust the maximum probability for this round */
	a_max = INT32_MIN;
	j_max = N_STATES;
	for (j = 0; j < N_STATES; ++j) {
		struct fix_pm_state * pm = &fix->map[j].pm;
		int32_t ln_pe;

		ln_pe = ln_pemission(pm, e, p_max,  0);

		if (ln_pe > a_max) {
			a_max = ln_pe;
			j_max = j;
		}

		alpha[j] = ln_pe;
		fix->beta[0][j] = N_STATES;
	}


	/* Save p_max for next round */
	p_max = a_max;
	/* Accumulate max probability */
	fix->sum_log_p_max += a_max;

	//--------------------------------------------------------------------------
	//allocate event independent variables
	//--------------------------------------------------------------------------
	//pm
	int32_t pm_rearranged[7*N_STATES]; size_t pm_rearranged_size= sizeof(int32_t) * 7*N_STATES;

	for (j = 0; j < N_STATES; ++j) {
		struct fix_pm_state * pm = &fix->map[j].pm;

		unsigned group_no=j/warp_size;
		unsigned start_index=group_no*warp_size*7 +j%warp_size;

		pm_rearranged[start_index]=pm->sd_mean;
		pm_rearranged[start_index+1*warp_size]=pm->sd_mean_inv;
		pm_rearranged[start_index+2*warp_size]=pm->sd_lambda_p5;
		pm_rearranged[start_index+3*warp_size]=pm->level_mean;
		pm_rearranged[start_index+4*warp_size]=pm->level_stdv_inv_2;
		pm_rearranged[start_index+5*warp_size]=pm->log_sd_lambda_p5;
		pm_rearranged[start_index+6*warp_size]=pm->log_level_stdv_2pi;

	}

	//use page lock memory for pm
	int32_t * devPtr_pm=assign_page_locked_memory_int32_t(pm_rearranged, pm_rearranged_size);


	//log_pr
	int32_t log_pr_rearranged[21*N_STATES]; size_t log_pr_rearranged_size= sizeof(int32_t) *21* N_STATES;
	for (j = 0; j < N_STATES; ++j) {

		unsigned group_no=j/warp_size;
		unsigned start_index=group_no*warp_size*21 +j%warp_size;

		for(k=0;k<21;k++){
			log_pr_rearranged[start_index+ k*warp_size]=fix->map[j].log_pr[k];
		}
	}
	/*cout<<"log_pr from CPU"<<endl;
	for (j = 0; j < 10; ++j) {
		unsigned group_no=j/warp_size;
		unsigned start_index=group_no*warp_size*21 +j%warp_size;

		cout<<"  " << fix->map[j].log_pr[0]<<"   "<<log_pr_rearranged[start_index]<<endl;

	}
	cout<<endl<<endl;*/

	//use page-locked memory
	int32_t * devPtr_log_pr=assign_page_locked_memory_int32_t(log_pr_rearranged, log_pr_rearranged_size);



	//alpha
	int32_t *dev_alpha;
	cudaMalloc((void**)&dev_alpha, sizeof(int32_t)*N_STATES) ; // device
	//keep alpha in global memory


	//beta needs a write only memory : so global memory
	uint16_t * dev_beta;
	cudaMalloc((void**)&dev_beta, sizeof(uint16_t)*N_STATES) ; // device
	uint16_t * h_beta;
	h_beta=(uint16_t*)malloc(sizeof(uint16_t)*N_STATES);



	// a_max and j_max;
	int32_t *h_amax;
	h_amax = (int32_t*)malloc(sizeof(int32_t)); //host
	int32_t *d_amax;
	cudaMalloc((void**)&d_amax, sizeof(int32_t));//device

	unsigned int *h_jmax;
	h_jmax = (unsigned int*)malloc(sizeof(unsigned int)); //host
	unsigned int *d_jmax;
	cudaMalloc((void**)&d_jmax, sizeof(unsigned int));//device

	int32_t *d_mutex;
	cudaMalloc((void**)&d_mutex, sizeof(int32_t));
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------
	//--------------------------------------------------------------------------

	for (i = 1; i < fix->n_events; ++i) {

		if ((i % TRACEBACK_CNT) == 0)
			fixpt_traceback(fix, i, j_max); // need to update beta




		preprocess_for_parallel(fix, alpha, p_max, i, devPtr_pm, devPtr_log_pr, dev_alpha, dev_beta, h_beta, h_amax, d_amax, h_jmax, d_jmax, d_mutex );

		a_max=*h_amax; j_max=*h_jmax;



		/* Save p_max for next round */
		p_max = a_max;
		/* Accumulate max probability */
		fix->sum_log_p_max += a_max;



	}

	fixpt_traceback(fix, i, j_max); // need to update beta

	//DBG(DBG_INFO, "end");
	cout<<"end"<<endl;

}


int64_t tcbk_fill_state_seq(struct bcall_tcbk * fix)
{
	int64_t ret;

	ret = fix->sum_log_p_max;
#if 0
	free(fix);
#endif
	return ret;
}

char * tcbk_fill_base_seq(struct bcall_tcbk * fix)
{
	int base_cnt;

	base_cnt = encode_base_seq(fix->base_seq, fix->state_seq, fix->n_events);
	fix->base_cnt = base_cnt;

	return fix->base_seq;
}

void tcbk_blk_free(struct bcall_tcbk * fix)
{
	free(fix);
}



/*
__global__ void kernel(int32_t *devPtr_pm, int32_t *devPtr_log_pr, int32_t *alpha, uint16_t *beta, int32_t p_max, int32_t *d_amax, unsigned int *d_jmax, int32_t *d_mutex, int ith_event){

	unsigned id = threadIdx.x + blockIdx.x * threads_per_block;   //id is j from the serial code

	unsigned group_block=id/N_STATES;
	unsigned start_block=group_block* number_of_block_per_group;
	unsigned end_block= start_block+number_of_block_per_group -1;


	__syncthreads();

	//const int number_of_iteration=ceil((float)N_STATES/threads_per_block);
	/*if (threadIdx.x==0){
		//printf("block id= %d and number_of_iteration=%d\n",blockIdx.x, number_of_iteration);
		printf("block id= %d\n",blockIdx.x);
	}*/

/*
	int32_t ln_pe_list[number_of_iteration];
	__syncthreads();
	unsigned number_of_elements_in_ln_pe_list=0;
	for (int k=0;k<number_of_iteration;k++){
		unsigned temp_id=id+k*threads_per_block;
		if( temp_id < N_STATES){
			//TD DO: try to move alpha to shared memory so that it can save time accessing global memory
			int32_t ln_pt=dev_ln_ptransition(devPtr_log_pr, temp_id, alpha, beta);


			int32_t ln_pe = dev_ln_pemission(devPtr_pm, p_max,  ln_pt, temp_id);


			alpha[temp_id] = ln_pe;
			ln_pe_list[k]=ln_pe; number_of_elements_in_ln_pe_list++;

			if(ith_event==1)
				if(threadIdx.x==0)
					printf("k=%d and ln_pt=%d \n",k,ln_pt);
		
	

		}
	}
	if(ith_event==1){
		if(threadIdx.x==0){
			for (int k=0;k<number_of_iteration;k++){
				printf("k=%d and ln_pe=%d \n",k,ln_pe_list[k]);

			}
		}
	}

	int32_t max_ln_pe_in_thread=ln_pe_list[0];
	unsigned max_temp_id_in_thread=id;
	for(int k=1; k<number_of_elements_in_ln_pe_list; k++ ){
		if (max_ln_pe_in_thread<ln_pe_list[k]){
			 max_ln_pe_in_thread=ln_pe_list[k];
			 max_temp_id_in_thread=id+k*threads_per_block;
		}
	}


	if(ith_event==1){
		if(threadIdx.x==0){
			printf("max_ln_pe_in_thread=%d and max_temp_id_in_thread=%d \n",max_ln_pe_in_thread,max_temp_id_in_thread);
		}
	}
	__syncthreads();



	//find max ln_pe and the corresponding id
	__shared__ int cache[threads_per_block];
	cache[threadIdx.x]=max_ln_pe_in_thread;
	__syncthreads();


	// reduction
	unsigned int i = blockDim.x/2;
	while(i != 0){
		if(threadIdx.x < i){
			cache[threadIdx.x] = max(cache[threadIdx.x], cache[threadIdx.x + i]);
		}

		__syncthreads();
		i /= 2;
	}
	__syncthreads();



	/*if(cache[0]==max_ln_pe_in_thread){
		//while(atomicCAS(d_mutex,0,1) != 0);  //lock
		*d_amax=cache[0];
		*d_jmax=  max_temp_id_in_thread;
		//atomicExch(d_mutex, 0);  //unlock
	}

	__syncthreads();*/
	
/*	if(blockIdx.x >= start_block && blockIdx.x <= end_block){

		if(threadIdx.x == 0){
			while(atomicCAS(d_mutex,0,1) != 0);  //lock
			*d_amax = max(*d_amax, cache[0]);
			atomicExch(d_mutex, 0);  //unlock
		}
	}


	//---------------------------------------------------------------


	__syncthreads();
	
	if(ln_pe_list[0]== *d_amax){
		
		*d_jmax=id;
		
	}
	__syncthreads();
	


}*/




#endif /* VITERBI_TCBK_V1_MULTI_EVENTS_CUH_ */
