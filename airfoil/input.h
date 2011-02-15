void input(int &nnode, int &ncell, int &nedge, float * &x, float * &q,
           int * &cell, int * &edge, int * &ecell, int * &boun) {

  float p,r,u,e, xt,yt;

// set global constants

  gam = 1.4f;
  gm1 = gam - 1.0f;
  cfl = 0.9f;
  eps = 0.05f;

// read in data from grid file

  printf("reading in grid data\n");

  FILE *fp;
  fp = fopen("grid.dat","r");

  fscanf(fp,"%d %d %d \n",&nnode, &ncell, &nedge);
  ncell = ncell+1;   // far-field dummy cell

  x     = (float *) malloc(2*nnode*sizeof(float));
  q     = (float *) malloc(4*ncell*sizeof(float));
  cell  = (int *) malloc(4*ncell*sizeof(int));
  edge  = (int *) malloc(2*nedge*sizeof(int));
  ecell = (int *) malloc(2*nedge*sizeof(int));
  boun  = (int *) malloc(  nedge*sizeof(int));

  for (int n=0; n<nnode; n++)
    fscanf(fp,"%f %f \n",&x[2*n], &x[2*n+1]);

  for (int n=0; n<ncell-1; n++) {
    fscanf(fp,"%d %d %d %d \n",&cell[4*n], &cell[4*n+1], &cell[4*n+2], &cell[4*n+3]);
    for(int m=0; m<4; m++) cell[4*n+m]--;             // FORTRAN->C offset
  }
  for (int n=0; n<4; n++) cell[4*(ncell-1)+n] = 0;    // far-field dummy cell

  for (int n=0; n<nedge; n++) {
    fscanf(fp,"%d %d %d %d %d \n",&edge[2*n],&edge[2*n+1],&ecell[2*n],&ecell[2*n+1],&boun[n]);
    for(int m=0; m<2; m++) edge[2*n+m]--;             // FORTRAN->C offset
    for(int m=0; m<2; m++) ecell[2*n+m]--;            // FORTRAN->C offset

    if (boun[n]==1) {
      ecell[2*n] = ecell[2*n+1];  // set interior dummy cell to point to one outside
    }
    else if (boun[n]==2) {
      ecell[2*n+1] = ncell-1;     // set exterior dummy cell to point to far-field cell
      boun[n] = 0;                // then treat like a regular edge
    }
  }

  fclose(fp);

// read in data from flow file, initialising if necessary

  fp = fopen("flow.dat","r");
  fscanf(fp,"%f %f %f %f \n",&p,&r,&mach,&alpha);
  alpha = alpha*atan(1.0f)/45.0f;
  p = 1.0f;
  r = 1.0f;
  u = sqrt(gam*p/r)*mach;
  e = p/(r*gm1) + 0.5f*u*u;

  for (int n=0; n<ncell; n++) {
    q[4*n  ] = r;
    q[4*n+1] = r*u;
    q[4*n+2] = 0.0f;
    q[4*n+3] = r*e;
  }

  fclose(fp);

// rotate grid to specified angle of attack

  for (int n=0; n<nnode; n++) {
    xt = x[2*n  ];
    yt = x[2*n+1];
    x[2*n  ] = cos(alpha)*xt + sin(alpha)*yt;
    x[2*n+1] = cos(alpha)*yt - sin(alpha)*xt;
  }
}
