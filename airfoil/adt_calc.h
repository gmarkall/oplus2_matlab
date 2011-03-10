inline void adt_calc(op_data<float,2,4> x,op_data<float,4> q,op_data<float> adt){
  float dx,dy, ri,u,v,c;

  ri =  1.0f/q(0);
  u  =   ri*q(1);
  v  =   ri*q(2);
  c  = sqrt(gam*gm1*(ri*q(3)-0.5f*(u*u+v*v)));

  dx = x(0,1) - x(0,0);
  dy = x(1,1) - x(1,0);
  (*adt)  = fabs(u*dy-v*dx) + c*sqrt(dx*dx+dy*dy);

  dx = x(0,2) - x(0,1);
  dy = x(1,2) - x(1,1);
  (*adt) += fabs(u*dy-v*dx) + c*sqrt(dx*dx+dy*dy);

  dx = x(0,3) - x(0,2);
  dy = x(1,3) - x(1,2);
  (*adt) += fabs(u*dy-v*dx) + c*sqrt(dx*dx+dy*dy);

  dx = x(0,0) - x(0,3);
  dy = x(1,0) - x(1,3);
  (*adt) += fabs(u*dy-v*dx) + c*sqrt(dx*dx+dy*dy);

  (*adt) = (*adt) / cfl;
}


