//### divideSimplex:
//### S wird geteilt, Aufruf durch T
void divideSimplex (simplex* S, simplex* T, mesh* M){
	int n = S->dim;
	int i,p;
	int gamma;
	node* newnode;
	if (S->child0 != NULL || isundefBd(S)){//### Kein Teilen nötig
		return;
	}
	
	if (n == M->dim){//### Randflaechen duerfen keine Teilung initiieren
		for (i=1;i<n;i++){
			if (S->neighbor[i]->dim == M->dim){//### Randflaechen benoetigen keine Mehrfachteilung 
				if (S->neighbor[i]->generation < S->generation){
					divideSimplex(S->neighbor[i],NULL,M);
				}
			}
		}
	}
	if (T==NULL){//### Erster Tetraeder des Edgepatches
		newnode = generateMidpoint(S->node[0],S->node[n]);
		newnode->num = M->lastNode->num+1;
		newnode->next = NULL;
		M->lastNode->next = newnode;
		M->lastNode = newnode;
		M->nC = M->nC+1;
	}else{
		newnode = T->child0->node[1];
	}
	S->child0 = newSimplex(S,S->generation+1);
	S->childn = newSimplex(S,S->generation+1);
	if (S->boundaryNum == 0){
		(M->nE)++;
	}else{
		(M->nBE[S->boundaryNum-1])++;
	}
	//### Verkettete Liste aktualisieren
	if (S->before != NULL){
		S->before->next = S->child0;
	}
	S->child0->before = S->before;
	S->child0->next = S->childn;
	S->childn->before = S->child0;
	S->childn->next = S->next;	
	if(S->next!=NULL){
		S->next->before = S->childn;
	}
	//### Knoten anpassen
	S->child0->node[0] = S->node[0];
	S->childn->node[0] = S->node[n];

	S->child0->node[1] = newnode;
	S->childn->node[1] = newnode;
	
	gamma = (S->generation)%n;
	for (p=1;p<n;p++){
		int pn=p+1;
		int p0=p+1;
		if (p>gamma){
			pn = n-p+gamma+1;
		}
		S->childn->node[pn] = S->node[p];
		S->child0->node[p0] = S->node[p];
		if (n == M->dim){
			if (S->neighbor[p]->child0 != NULL){
				if (S->neighbor[p]->node[0]==S->node[0]){
					S->child0->neighbor[p0] = S->neighbor[p]->child0;
					S->childn->neighbor[pn] = S->neighbor[p]->childn;
				}else{
					S->child0->neighbor[p0] = S->neighbor[p]->childn;
					S->childn->neighbor[pn] = S->neighbor[p]->child0;
				}
			}else{
				S->childn->neighbor[pn] = S->neighbor[p];
				S->child0->neighbor[p0] = S->neighbor[p];
			}
		}
	}
	if (n == M->dim){//### Falls Rand, bleibt neighbor unangetastet
		S->child0->neighbor[0] = S->childn;
		S->childn->neighbor[0] = S->child0;
		
		if (S->neighbor[n]->child0 != NULL){
			if (hasSameRefEdge(S->child0,S->neighbor[n]->child0)){// S->child0 passt mit S->neighbor[n]->child0 zusammen
				S->child0->neighbor[1] = S->neighbor[n]->child0;
			}else{// S->child0 passt mit S->neighbor[n]->childn zusammen
				S->child0->neighbor[1] = S->neighbor[n]->childn;
			}
		}else{
			S->child0->neighbor[1] = S->neighbor[n];
		}
		if (S->neighbor[0]->child0 != NULL){
			if (hasSameRefEdge(S->childn,S->neighbor[0]->childn)){// S->childn passt mit S->neighbor[0]->childn zusammen
				S->childn->neighbor[1] = S->neighbor[0]->childn;				
			}else{//S passt mit S->neighbor->child0 zusammen
				S->childn->neighbor[1] = S->neighbor[0]->child0;
			}
		}else{
			S->childn->neighbor[1] = S->neighbor[0];
		}
		for (i=0;i<M->dim+1;i++){
			if (S->child0->neighbor[i]->generation==S->child0->generation){
				S->child0->neighbor[i]->neighbor[i] = S->child0;
			}
			if (S->childn->neighbor[i]->generation==S->childn->generation){
				S->childn->neighbor[i]->neighbor[i] = S->childn;
			}
		}
		for (i=1;i<n;i++){
			divideSimplex(S->neighbor[i],S,M);
		}
	}else{//### Hier ist S Randflaeche und wurde von T aufgerufen
		// Die Kinder von T muessen auf soeben erstellte Randflaechen verweisen.
		for (i=0;i<M->dim+1;i++){
			if (T->child0->neighbor[i]==S){
				if (T->child0->node[0]==S->node[0])
					T->child0->neighbor[i] = S->child0;
				else
					T->child0->neighbor[i] = S->childn;
			}
			if (T->childn->neighbor[i]==S){
				if (T->childn->node[0]==S->node[0])
					T->childn->neighbor[i] = S->child0;
				else
					T->childn->neighbor[i] = S->childn;
			}
		}
	}
}
