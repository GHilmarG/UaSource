typedef struct node{
    int num;
    double co[3];
    struct node* next;
} node;

typedef struct simplex{
    int boundaryNum;/*boundarNum = 
   -1: Rand, aber nicht uebergeben worden, 
    0: kein Rand sondern Tetraeder, 
    1: boundary[0],
    2: boundary[1], ...*/
    int generation;
    int dim;
    node** node;
    struct simplex** neighbor;/* neighbor bleibt bei Randflaechen leer*/
    struct simplex* child0;
    struct simplex* childn;
    struct simplex* parent;
	struct simplex* next;
	struct simplex* before;
} simplex;

typedef struct mesh{
	int dim;
	int nE;
	int nB;
	int* nBE;
	int nC;
	node* firstNode;
	node* lastNode;
	simplex* firstElement;
	simplex** firstBdElement; /* firstBdElement[0] pointer auf ersten boundary[0]*/
	simplex** todo;
	int todoSize;
} mesh;

node* newNode(double cx, double cy, double cz){
    node* newnode = (node*)mxMalloc(sizeof(node));
    newnode->co[0] = cx;
    newnode->co[1] = cy;
    newnode->co[2] = cz;
    return newnode;
}

node* generateMidpoint(node* a, node* b){
    return newNode((a->co[0]+b->co[0])/2,(a->co[1]+b->co[1])/2,(a->co[2]+b->co[2])/2);
}

simplex* newSimplex(simplex* parent, int generation){
	int i;
    simplex* newsimplex = (simplex*)mxMalloc(sizeof(simplex));
    newsimplex->dim = 0;
    newsimplex->boundaryNum = 0;
    newsimplex->child0 = NULL;
    newsimplex->childn = NULL;
    newsimplex->neighbor = NULL;
    newsimplex->parent = NULL;
	newsimplex->next = NULL;
	newsimplex->before = NULL;
	newsimplex->node = NULL;
    if (parent!=NULL){
	    newsimplex->node=(node**)mxMalloc((parent->dim+1)*sizeof(node*));
		for (i=0;i<parent->dim+1;i++){
			newsimplex->node[i]=NULL;
		}
	    newsimplex->parent=parent;
	    newsimplex->dim = parent->dim;
	    newsimplex->boundaryNum = parent->boundaryNum;
	    if (parent->neighbor!=NULL){/* D.h. kein Randelement*/
	    	newsimplex->neighbor = (simplex**) mxMalloc((newsimplex->dim+1)*sizeof(simplex*));
			for (i=0;i<parent->dim+1;i++){
				newsimplex->neighbor[i]=NULL;
			}
	    }
	}
	newsimplex->generation = generation;
    return newsimplex;
}
simplex* genUndefBd(){
    simplex* newsimplex = (simplex*) mxMalloc(sizeof(simplex));
    newsimplex->dim = 0;
	newsimplex->boundaryNum = -1;
	newsimplex->generation = -1;
    newsimplex->child0 = NULL;
    newsimplex->childn = NULL;
    newsimplex->parent = NULL;
	newsimplex->neighbor = NULL;
	newsimplex->next = NULL;
	newsimplex->before = NULL;
	newsimplex->node=NULL;
    return newsimplex;
}

int isundefBd(simplex* S){
	if (S!=NULL){
		return (S->boundaryNum==-1);
	}else{
		return 1;
	}
}
int hasSameRefEdge(simplex* S, simplex* T){
	if ((S->node[0]==T->node[0]) && (S->node[S->dim]==T->node[S->dim]  )){
		return 1;
	}
	if ((S->node[0]==T->node[S->dim]) && (S->node[S->dim]==T->node[0] )){
		return 1;
	}
	return 0;
}
