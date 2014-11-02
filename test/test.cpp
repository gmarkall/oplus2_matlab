#include <iostream>

//
// OP header file
//

#include "op_seq.h"

#define NN       2

int cell = 0;

// C-style kernel computing the center of mass of an element

extern "C" {

void kernel( const op_data<float,4,2> x, op_data<float,2> c) {
  c(0) = (x(0,0) + x(1,0) + x(2,0) + x(3,0))*0.25f;
  c(1) = (x(0,1) + x(1,1) + x(2,1) + x(3,1))*0.25f;
  std::cout << "cell " << cell++ << ": (" << x(0,0) << "," << x(0,1) << ") (" << x(1,0) << "," << x(1,1) << ") (" << x(2,0) << "," << x(2,1) << ") (" << x(3,0) << "," << x(3,1) << ")" << std::endl;
}

}

int main(int argc, char **argv) {

  float x_d[2*NN*NN];
  float c_d[(NN-1)*(NN-1)];
  int m[4*(NN-1)*(NN-1)];

  int nnode = NN*NN;
  int ncell = (NN-1)*(NN-1);
  float dx    = 1.0f / ((float) NN-1);

  // Coordinates
  for (int i = 0; i < NN; ++i) {
    for (int j = 0; j < NN; ++j) {
      x_d[2*(j+i*NN)] = i*dx;
      x_d[2*(j+i*NN)+1] = j*dx;
    }
  }

  // Mapping
  for (int i = 0; i < NN-1; ++i) {
    for (int j = 0; j < NN-1; ++j) {
      m[4*(j+i*(NN-1))]   = j  +i*NN;
      m[4*(j+i*(NN-1))+1] = j+1+i*NN;
      m[4*(j+i*(NN-1))+2] = j  +(i+1)*NN;
      m[4*(j+i*(NN-1))+3] = j+1+(i+1)*NN;
    }
  }

  op_set nodes, cells;
  op_map ppcell;
  op_dat x, c;

  // OP initialisation

  op_init(argc,argv,5);

  // declare sets, pointers, and datasets

  op_decl_set(nnode, nodes,"nodes");
  op_decl_set(ncell, cells,"cells");

  op_decl_map(cells,nodes,4,m, ppcell,"ppcell");

  op_decl_dat(cells,1,"float", c_d, c, "c" );
  op_decl_dat(nodes,2,"float", x_d, x, "x" );

  op_par_loop(kernel,"kernel", cells,
              x, OP_ALL,ppcell, 2,"float", OP_READ,
              c,OP_NONE, OP_ID, 1,"float", OP_WRITE);

  for (int i = 0; i < NN-1; ++i) {
    for (int j = 0; j < NN-1; ++j) {
      std::cout << "(" << c_d[2*(j+i*(NN-1))] << "," << c_d[2*(j+i*(NN-1))+1] << ") ";
    }
    std::cout << std::endl;
  }

  return 0;
}
