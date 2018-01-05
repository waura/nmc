#ifndef _NMC_UCU_FUNC_H_
#define _NMC_UCU_FUNC_H_

__host__ __device__ static __inline__ float binomial_coefficientf(int n, int m)
{
#ifdef __DEVICE_EMULATION__
	if(n < 0) printf("error: binomial_coefficientf\n");
	if(m < 0) printf("error: binomial_coefficientf\n");
#endif
	if(n == m) return 1.0;
	if(m == 0) return 1.0;

	int i;
	float nn=1.0, nm=1.0, mm=1.0;
	for(i=1; i<=n; i++){
		nn *= i;
	}
	for(i=1; i<=(n-m); i++){
		nm *= i;
	}
	for(i=1; i<=m; i++){
		mm *= i;
	}

	return nn /(nm*mm);
}

#endif //_NUM_UCU_FUNC_H_