#include "nmc/def.h"
#include "nmc/LU.h"
#include "nmc/CSM.h"
#include "nmc/BiCGSTAB.h"
#include "nmc/Vector2D.h"
#include "nmc/FMCSM2D.h"
#include "nmc/FMCSM2D_Helmholtz_Kernel.h"
#include "nmc/FMCSM2D_Laplace_Kernel.h"


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

class CSMKernelOp {
public:
	nmc::Complexd operator()(const nmc::Vector2Dd& spt, const nmc::Vector2Dd& cpt)
	{
		double r = sqrt((spt[0]-cpt[0])*(spt[0]-cpt[0]) + (spt[1]-cpt[1])*(spt[1]-cpt[1]));
		nmc::Complexd ret(log(r), 0.0);
		return ret;
	}
};

class CSMBoundaryOp {
public:
	nmc::Complexd operator()(const nmc::Vector2Dd& cpt)
	{
		nmc::Complexd ret(cpt[0]*cpt[0] - cpt[1]*cpt[1]);
		return ret;
	}
};

bool isInArea(double rho, double x, double y){
	if(x*x + y*y > rho*rho){
		return true;
	}
	return false;
}

int main(int argc, char* arg[]){
	int i,j;
	//int n=atoi(arg[1]);
	//double rho = atof(arg[2]);
	//double q = atof(arg[3]);
	//char* cmd = arg[4];
	int n=1000;
	double rho =1.0;
	double q = 0.2;
	char* cmd = "dif";

	//char* maxN_eq_file;
	//if(argc > 5){
	//	maxN_eq_file = arg[5];
	//}

	nmc::Vector< nmc::Vector2Dd > collocation_points(n);
	nmc::Vector< nmc::Vector2Dd > source_points(n);

	//電荷点・拘束点配置
	for(i=0; i<n; i++){
		nmc::Vector2Dd t1;
		nmc::Vector2Dd t2;
		t1.SetVal(0, rho * cos(i*(2.0*NMC_PI)/n));
		t1.SetVal(1, rho * sin(i*(2.0*NMC_PI)/n));
		collocation_points.SetVal(i, t1);
		t2.SetVal(0, q * t1.GetVal(0));
		t2.SetVal(1, q * t1.GetVal(1));
		source_points.SetVal(i, t2);
	}

	CSMBoundaryOp csm_boundaryOp;
	CSMKernelOp csm_kernelOp;
	nmc::Matrixcmd A,B;
	nmc::Vectorcmd f_s;

	nmc::CSM2Dcmd::CreateMatrix<CSMBoundaryOp, CSMKernelOp>(
		csm_boundaryOp,
		csm_kernelOp,
		source_points,
		collocation_points,
		A,
		f_s);

	//A.dump(); printf("\n");
	//f_s.dump(); printf("\n");

	nmc::Vectori P(A.row());
	nmc::Vectorcmd Q(A.row());
	B = A;

	nmc::LUcmd::Decomp(A, P);
	nmc::LUcmd::Solve(A, P, f_s, Q);
	std::cout << Q << std::endl;

	//double conv_ratio = 1.0e-16;
	//unsigned int iteration = n;
	//for(i=0; i<n; i++){
	//	Q[i] = 0.0;
	//}
	//Vectorcmd tmp(n);
	//tmp = A * Q;
	//for(i=0; i<n; i++){
	//	f_s[i] -= tmp[i];
	//}
	//BiCGSTAB<nmc::Complexd>::Solve(conv_ratio, iteration, A, f_s, Q);

	//f_s.dump(); printf("\n");

	/*
	double conv_ratio = 1.0e-10;
	unsigned int iteration = n;
	for(i=0; i<n; i++){
		Q[i] = 0.0;
	}
	*/

	nmc::Vectorcmd Ax(n);
	//Ax = B * Q;
	//Ax.dump(); printf("\n");
	//nmc::FMCSM2D_Helmholtz_Kernel kernel;
	//kernel.SetWaveNumber(1.0);
	nmc::FMCSM2D_Laplace_Kernel kernel;
	nmc::FMCSM2DTree<nmc::Complexd>* tree = nmc::FMCSM2D<nmc::Complexd>::CreateTree(
		3, source_points, collocation_points);
	kernel.Init( tree );
	////nmc::FMCSM2D<nmc::Complexd>::PlotTreeBox( tree );
	////nmc::FMCSM2D<nmc::Complexd>::PlotSourcePoint( tree );
	////nmc::FMCSM2D<nmc::Complexd>::PlotBundaryPoint( tree );
	nmc::FMCSM2D<nmc::Complexd>::Evaluate(tree, kernel, Q, Ax);
	//nmc::FMCSM2D<nmc::Complexd>::Solve(
	//	tree,
	//	kernel,
	//	FMCSM_SOLVER_BiCGSTAB,
	//	conv_ratio,
	//	iteration,
	//	f_s,
	//	Q);
	nmc::FMCSM2D<nmc::Complexd>::DestoryTree( tree );
	//f_s.dump(); printf("\n");
	//Q.dump(); printf("\n");
	//printf("%d \n", iteration);
	std::cout << Ax << std::endl;

	if(strcmp("dif", cmd) == 0){
		double maxDif = 0.0;
		nmc::Complexd tmp;
		nmc::Vector2Dd pt;
		for(i=0; i<1000; i++){
			pt[0] = rho * cos(i*(2.0*NMC_PI)/1000.0);
			pt[1] = rho * sin(i*(2.0*NMC_PI)/1000.0);
			
			tmp = nmc::CSM2Dcmd::u_N<CSMKernelOp>(csm_kernelOp, Q, source_points, pt);
			tmp -= csm_boundaryOp(pt);
			double dif = tmp.Abs(); //dif = |u_N - f|

			if(dif > maxDif){
				maxDif = dif;
			}	
		}
		//
		printf("%f %f %d %.15f \n", rho, q, n, log10(maxDif));
	}
	else if(strcmp("dif-ext2", cmd) == 0){
	  /*
	        srand(time(NULL));
		bool is_in_area;
		double maxDif = 0.0;
		nmc::Complexd tmp;
		for(i=0; i<1000; i++){
			double x, y;
			is_in_area=false;
			while(!is_in_area){
				x = (rand()%210 - 100)/10.0;
				y = (rand()%210 - 100)/10.0;
				if(isInArea(rho, x, y)){
					
					double r = sqrt(x*x + y*y);
					double theta = atan2(y, x);
					tmp = u_N(k, Q, source_points, x, y);
					tmp -= u(rho, m, k, r, theta);
					double dif = tmp.Abs();

					if(dif > maxDif){
						maxDif = dif;
					}
					is_in_area = true;
				}
			}
		}
		//
		printf("%f %f %d %.15f \n", rho, q, n, log10(maxDif));
	  */
	}
	else if(strcmp("cond", cmd) == 0){
		nmc::Complexd cm = nmc::LUcmd::Cond(A, P, A.l1_norm());
		printf("%f %f %d %.15f \n", rho, q, n, log10(static_cast<double>(cm.Re())));
	}
	else if(strcmp("plot", cmd) == 0){
		double x, y;
		printf("# x  y  z \n");
		nmc::Vector2Dd pt;
		for(x=-20.0; x<=20.0; x+=0.12){
			for(y=-20.0; y<=20.0; y+=0.12){
				if(isInArea(rho, x, y)){
				  pt[0] = x;
				  pt[1] = y;
				  nmc::Complexd cm = nmc::CSM2Dcmd::u_N<CSMKernelOp>(csm_kernelOp, Q, source_points, pt);
				  printf("%f %f %.13f \n", (double)x, (double)y, cm.Re());
				}
				else{
					printf("%f %f -4000000000 \n",(double)x, (double)y);
				}	
			}
			printf("\n");
		}
	}

	return 0;
}
