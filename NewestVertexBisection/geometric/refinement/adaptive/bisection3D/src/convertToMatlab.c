//### convertToMatlab
void convertToMatlab(int nlhs, mxArray *plhs[], mesh* M){
	int k,i,j;
	double* elements;
	double* elementgeneration;
	double* coordinates;
	double* boundary;
	simplex* aktsimplex;
	simplex* aktboundary;
	node* aktnode;
	//### elements
	if(nlhs>=1){
	    plhs[0] = mxCreateDoubleMatrix(M->nE,M->dim+1,mxREAL);
	    elements = mxGetPr(plhs[0]);
	    // 
		i = 0;
	   	for(j=0;j<M->dim+1;j++){
	    	aktsimplex = M->firstElement;
		    while(aktsimplex!=NULL){
		    	elements[i] = aktsimplex->node[j]->num;
		    	aktsimplex = aktsimplex->next;
		    	i++;
		    }
		}
	}
	//### elementgeneration
	if(nlhs>=2){
	    plhs[1] = mxCreateDoubleMatrix(M->nE,1,mxREAL);
	    elementgeneration = mxGetPr(plhs[1]);
	    aktsimplex = M->firstElement;
	    i=0;
	    while(aktsimplex!=NULL){
	    	elementgeneration[i] = aktsimplex->generation;
	    	aktsimplex = aktsimplex->next;
	    	i++;
	    }
	}
	//### coordinates
	if(nlhs>=3){
   		plhs[2] = mxCreateDoubleMatrix(M->nC,3,mxREAL);
    	coordinates = mxGetPr(plhs[2]);
		i=0;
	   	for(j=0;j<3;j++){
			aktnode = M->firstNode;
		    while(aktnode!=NULL){
		    	coordinates[i] = aktnode->co[j];
		    	aktnode = aktnode->next;
		    	i++;
		    }
		}
	}
    //### boundaries
   	for (k=0;k<M->nB;k++){
   		if(nlhs>=4+k){
   			plhs[3+k] = mxCreateDoubleMatrix(M->nBE[k],M->dim,mxREAL);
			boundary = mxGetPr(plhs[3+k]);
			i = 0;
		   	for(j=0;j<M->dim;j++){
		    	aktboundary = M->firstBdElement[k];
			    while(aktboundary!=NULL){
			    	boundary[i] = aktboundary->node[j]->num;
			    	aktboundary = aktboundary->next;
			    	i++;
			    }
			}
   		}
   	}
}