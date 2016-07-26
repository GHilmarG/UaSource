/* 
   Copyright (c) 2012 by Marcin Krotkiewski, University of Oslo
   See ../License.txt for License Agreement.
*/

#include "mesh.h"
#include "mexparams.h"

#ifdef MATLAB_MEX_FILE

t_mesh mex2mesh(const mxArray *mesh_struct){

  mxArray *field;
  t_mesh mesh = {0};
  int i;

  if(!mxIsStruct(mesh_struct)){
    USERERROR("mesh_struct is not a structure", MUTILS_INVALID_MESH);
  }

  /* ELEMS */
  mesh.n_elems = 0;
  mesh.n_elem_nodes = 0;
  field = mxGetField(mesh_struct, 0, "ELEMS");
  mesh.elems = mex_get_matrix_Int(field, &mesh.n_elem_nodes, &mesh.n_elems, 
				  "MESH.ELEMS", "number of element nodes", 
				  "number of elements", 0);

  /* NODES */
  mesh.n_nodes = 0;
  mesh.n_dim = 2;
  field = mxGetField(mesh_struct, 0, "NODES");
  mesh.nodes = mex_get_matrix_Double(field, &mesh.n_dim, &mesh.n_nodes, 
				     "MESH.NODES", "2", 
				     "number of nodes", 0);

  /* NEIGHBORS */
  mesh.n_neighbors = 0;
  field = mxGetField(mesh_struct, 0, "NEIGHBORS");
  mesh.neighbors = mex_get_matrix_Int(field, &mesh.n_neighbors, &mesh.n_elems,
				      "MESH.NEIGHBORS", "number of element neighbors", 
				      "number of elements", 1);


  /* data validation */
  for(i=0; i<mesh.n_elems*mesh.n_elem_nodes; i++){
    if(mesh.elems[i]>mesh.n_nodes){
      USERERROR("Illegal mesh structure: ELEMS access non-existant node IDs.", MUTILS_INVALID_MESH);
    }
  }
  
  if(mesh.neighbors){
    for(i=0; i<mesh.n_elems*mesh.n_neighbors; i++){
      if(mesh.neighbors[i]>mesh.n_elems){
	USERERROR("Illegal mesh structure: NEIGHBORS access non-existant elemsent IDs.", MUTILS_INVALID_MESH);
      }    
    }
  }
  
  return mesh;
}

#endif /* MATLAB_MEX_FILE */
