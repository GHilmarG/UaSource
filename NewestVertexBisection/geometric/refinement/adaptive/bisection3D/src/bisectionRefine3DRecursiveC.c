#include "mex.h"
#include "structs.c"
#include "buildMesh.c"
#include "divideSimplex.c"
#include "refine.c"
#include "convertToMatlab.c"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int i;
	//### Element- und Koordinaten-Informationen
    double* elements = mxGetPr(prhs[0]);
	int nE = mxGetM(prhs[0]);
	int dim = mxGetN(prhs[0])-1;
    double* elementgeneration = mxGetPr(prhs[1]);
    double* coordinates= mxGetPr(prhs[2]);
    int nC = (int) mxGetM(prhs[2]);
    double* elem2neigh = mxGetPr(prhs[3]);
    //### Rand-Information
    int nB = (nrhs-5)/2;
    int* nBE = (int*) mxMalloc(nB*sizeof(int));
    double** elem2bdNumber = (double**) mxMalloc(nB*sizeof(double*));
    double** boundaries = (double**) mxMalloc(nB*sizeof(double*));
    //### Markierte Elemente
    double* marked = mxGetPr(prhs[nrhs-1]);
    int nM = mxGetM(prhs[nrhs-1])*mxGetN(prhs[nrhs-1]);    
    //### Netz
    mesh* M;
    //### Rand aufbauen
    for (i=0;i<nB;i++){
    	boundaries[i] = mxGetPr(prhs[4+i]);
    	nBE[i] = mxGetM(prhs[4+i]);
    	elem2bdNumber[i] = mxGetPr(prhs[4+nB+i]);
    }
  	//### Aufbauen, Verfeinern, Zurueckgeben
   	M = buildMesh(dim,nE,elements,elementgeneration,nC,coordinates,elem2neigh,nB,nBE,boundaries,elem2bdNumber,nM,marked);  
	refine(M);
	convertToMatlab(nlhs, plhs, M);
}
