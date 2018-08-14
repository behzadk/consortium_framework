#include "stdio.h"
#include <sstream>

#define NSPECIES 7
#define NPARAM 12
#define NREACT 2

#define D tex2D(param_tex, 0, tid)
#define KB_mccI tex2D(param_tex, 1, tid)
#define K_CV514 tex2D(param_tex, 2, tid)
#define K_MG tex2D(param_tex, 3, tid)
#define g_CV514 tex2D(param_tex, 4, tid)
#define g_MG tex2D(param_tex, 5, tid)
#define kA_AHL_1 tex2D(param_tex, 6, tid)
#define kBmax_mccI tex2D(param_tex, 7, tid)
#define mu_max_CV514 tex2D(param_tex, 8, tid)
#define mu_max_MG tex2D(param_tex, 9, tid)
#define nB_mccI tex2D(param_tex, 10, tid)
#define omega_mccI tex2D(param_tex, 11, tid)

struct myFex{
__device__ void operator()(int *neq, double *t, double *y, double *ydot){

int tid = blockDim.x * blockIdx.x + threadIdx.x;

// Order is: A_AHL_1, B_mccI, N_CV514, N_MG, N_XX, S_glu, S_leu, 

ydot[0] = + kA_AHL_1 * N_CV514+ kA_AHL_1 * y[3] - D * A_AHL_1;
ydot[1] =  + ( kBmax_mccI  * ( powf( A_AHL1 , nB_mccI ) / ( powf( KB__mccI, nB_mccI ) + powf( A_AHL1 , nB_mccI ) ) )) * y[2] + ( kBmax_mccI  * ( powf( A_AHL1 , nB_mccI ) / ( powf( KB__mccI, nB_mccI ) + powf( A_AHL1 , nB_mccI ) ) )) * y[3] + ( kBmax_mccI ) * y[4] - D * B_mccI;
ydot[2] = (( mu_max_CV514 * y[6] / K_CV514 + S_leu)  + omega_mccI * B_mccI- D) * N_CV514;
ydot[3] = (( mu_max_MG * y[5] / K_MG + S_glu) - D) * N_MG;
ydot[4] = (( mu_max_XX * y[5] / K_XX + S_glu) - D) * N_XX;
ydot[5] = D * (S0_glu - S_glu)- ( mu_MG * y[3] / g_MG )- ( mu_XX * y[4] / g_XX );
ydot[6] = D * (S0_leu - S_leu)- ( mu_CV514 * y[2] / g_CV514 );
}
};struct myJex{
    __device__ void operator()(int *neq, double *t, double *y, int ml, int mu, double *pd, int nrowpd){
return;
} 
 };