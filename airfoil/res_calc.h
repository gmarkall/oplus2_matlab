inline void res_calc(op_data<float,2,2> x,  op_data<float,4,2> q,
                     op_data<float,2> adt,op_data<float,4,2> res) {
  float dx,dy,mu, ri, p1,vol1, p2,vol2, f;

  dx = x(0,0) - x(0,1);
  dy = x(1,0) - x(1,1);

  ri   = 1.0f/q(0,0);
  p1   = gm1*(q(3,0)-0.5f*ri*(q(1,0)*q(1,0)+q(2,0)*q(2,0)));
  vol1 =  ri*(q(1,0)*dy - q(2,0)*dx);

  ri   = 1.0f/q(0,1);
  p2   = gm1*(q(3,1)-0.5f*ri*(q(1,1)*q(1,1)+q(2,1)*q(2,1)));
  vol2 =  ri*(q(1,1)*dy - q(2,1)*dx);

  mu = 0.5f*(adt(0)+adt(1))*eps;

  f = 0.5f*(vol1* q(0,0)         + vol2* q(0,1)        ) + mu*(q(0,0)-q(0,1));
  res(0,0) += f;
  res(0,1) -= f;
  f = 0.5f*(vol1* q(1,0) + p1*dy + vol2* q(1,1) + p2*dy) + mu*(q(1,0)-q(1,1));
  res(1,0) += f;
  res(1,1) -= f;
  f = 0.5f*(vol1* q(2,0) - p1*dx + vol2* q(2,1) - p2*dx) + mu*(q(2,0)-q(2,1));
  res(2,0) += f;
  res(2,1) -= f;
  f = 0.5f*(vol1*(q(3,0)+p1)     + vol2*(q(3,1)+p2)    ) + mu*(q(3,0)-q(3,1));
  res(3,0) += f;
  res(3,1) -= f;
}
