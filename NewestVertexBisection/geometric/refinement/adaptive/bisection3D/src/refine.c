//### Refine:
void refine(mesh* M){
	simplex* simpl;
	int i,k;
	for (i=0;i<M->todoSize;i++){
		divideSimplex(M->todo[i],NULL,M);
	}
	//### Neues firstElement finden
	simpl = M->firstElement;
	while (simpl->child0!=NULL){
		simpl = simpl->child0;
	}
	M->firstElement = simpl;
	//### neue firstBdElements finden
	for (k=0;k<M->nB;k++){
		if (M->nBE[k]!=0){
			simpl=M->firstBdElement[k];
			while (simpl->child0!=NULL){
				simpl = simpl->child0;
			}
			M->firstBdElement[k] = simpl;
		}
	}
}