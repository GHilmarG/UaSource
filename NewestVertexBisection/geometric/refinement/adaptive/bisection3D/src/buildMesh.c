//### buildMesh:
mesh* buildMesh(int dim, int nE, double* elements,double* elementgeneration, int nC,
				double* coordinates,double* elem2neigh,int nB,int*nBE,
				double** boundariesMAT,double** elem2bdNumber,int nM, double* marked){
	int i,j,k;
	mesh* M = (mesh*) mxMalloc(sizeof(mesh));
	simplex* simplices;
	simplex** boundaries;
	simplex* undefBd = genUndefBd();//### Dummy-Rand, um haeufige NULL-Abfragen zu vermeiden.
	node* nodes;
	
	M->dim = dim;
	M->nE = nE;
	M->nB = nB;
	M->nBE = (int*) mxMalloc(nB*sizeof(int));
    for (i=0;i<nB;i++){
    	M->nBE[i] = nBE[i];
    }
    M->nC = nC;
    //### Knoten anlegen
    nodes = (node*) mxMalloc(nC*sizeof(node));
	for (i=0;i<nC;i++){
		(&nodes[i])->co[0] = coordinates[i];
		(&nodes[i])->co[1] = coordinates[i+nC];
		(&nodes[i])->co[2] = coordinates[i+2*nC];
	}
	for (i=0;i<nC-1;i++){
		(&nodes[i])->num = i+1;
		(&nodes[i])->next = &(nodes[i+1]);
	}
	nodes[nC-1].num = nC;
	nodes[nC-1].next = NULL;
	M->firstNode = &(nodes[0]);
	M->lastNode = &(nodes[nC-1]);
    //### Simplizes anlegen
    simplices = (simplex*) mxMalloc(nE*sizeof(simplex));

    for (i=0;i<nE;i++){
    	simplices[i].boundaryNum = 0;
    	simplices[i].generation = elementgeneration[i];
    	simplices[i].dim = dim;
    	simplices[i].child0 = NULL;
    	simplices[i].childn = NULL;
    	simplices[i].parent = NULL;
    	simplices[i].neighbor = (simplex**) mxMalloc((dim+1)*sizeof(simplex*));
    	simplices[i].node = (node**) mxMalloc((dim+1)*sizeof(node*));
      	for (j=0;j<dim+1;j++){
      		int nachb = (int)elem2neigh[j*nE+i];
      		if (nachb!=0){
      			simplices[i].neighbor[j] = &(simplices[nachb-1]);
      		}else{
      			simplices[i].neighbor[j] = undefBd;
      		}
      		simplices[i].node[j] = &(nodes[(int)elements[j*nE+i]-1]);
    	}
    	
    }
   	simplices[0].before = NULL;
   	if (nE>1){
   		simplices[0].next = &(simplices[1]);
   	}
   	for (i=1;i<nE-1;i++){
   		simplices[i].next = &(simplices[i+1]);
   		simplices[i].before = &(simplices[i-1]);
   	}
   	if (nE>1){
   		simplices[nE-1].before = &(simplices[nE-2]);
   	}
   	simplices[nE-1].next = NULL;

    M->firstElement = &(simplices[0]);
	//### Raender anlegen
	if (nB!=0){
		boundaries = (simplex**) mxMalloc(nB*sizeof(simplex*));	    	
		M->firstBdElement = (simplex**) mxMalloc(nB*sizeof(simplex*));
	
		for (k=0;k<nB;k++){
			boundaries[k] = (simplex*) mxMalloc(nBE[k]*sizeof(simplex));
		    for (i=0;i<nBE[k];i++){
		    	boundaries[k][i].boundaryNum = k+1;
		    	boundaries[k][i].generation = 0;
		    	boundaries[k][i].dim = dim-1;
		    	boundaries[k][i].child0 = NULL;
		    	boundaries[k][i].childn = NULL;
		    	boundaries[k][i].parent = NULL;
		    	boundaries[k][i].neighbor = NULL;
		    	boundaries[k][i].node = (node**) mxMalloc(dim*sizeof(node*));
		      	for (j=0;j<dim;j++){
		      		boundaries[k][i].node[j] = &(nodes[(int)boundariesMAT[k][j*nBE[k]+i]-1]);
		    	}	    	
		    }
		    if (nBE[k]!=0){
			   	boundaries[k][0].before = NULL;
			   	if (nBE[k]>1){
			   		boundaries[k][0].next = &(boundaries[k][1]);
			   	}
			   	for (i=1;i<nBE[k]-1;i++){
			   		boundaries[k][i].next = &(boundaries[k][i+1]);
			   		boundaries[k][i].before = &(boundaries[k][i-1]);
			   	}
			   	if (nBE[k]>1){
			   		boundaries[k][nBE[k]-1].before = &(boundaries[k][nBE[k]-2]);
			   	}
		   		boundaries[k][nBE[k]-1].next = NULL;
		   	}
	
		    M->firstBdElement[k] = boundaries[k];
			//### Nachbarschaftsrelationen auf boundaries
			for(i=0;i<nE;i++){
				for(j=0;j<dim+1;j++){
					int randflaeche = (int) elem2bdNumber[k][i+j*nE];
					if (randflaeche!=0){
						simplices[i].neighbor[j] = &(boundaries[k][randflaeche-1]);
					}
				}
			}
		}
	}
	M->todoSize = nM;
	M->todo = (simplex**) mxMalloc(nM*sizeof(simplex*));

	for (i=0;i<nM;i++){
		M->todo[i] = &(simplices[(int)marked[i]-1]);
	}
	
	return M;
}
