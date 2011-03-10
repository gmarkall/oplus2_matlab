inline void bres_calc(op_data<float,2,2> x, op_data<float,4> q,
                      op_data<float> adt, op_data<float,4> res,op_data<int> bound) {
  float dx,dy,mu, ri, p1,vol1, p2,vol2, f;

  dx = x(0,0) - x(0,1);
  dy = x(1,0) - x(1,1);

  ri = 1.0f/q(0);
  p1 = gm1*(q(3)-0.5f*ri*(q(1)*q(1)+q(2)*q(2)));

  if ((*bound)==1) {
    res(1) += + p1*dy;
    res(2) += - p1*dx;
  }
  else {
    vol1 =  ri*(q(1)*dy - q(2)*dx);

    ri   = 1.0f/qinf[0];
    p2   = gm1*(qinf[3]-0.5f*ri*(qinf[1]*qinf[1]+qinf[2]*qinf[2]));
    vol2 =  ri*(qinf[1]*dy - qinf[2]*dx);

    mu = (*adt)*eps;

    f = 0.5f*(vol1* q(0)         + vol2* qinf[0]        ) + mu*(q(0)-qinf[0]);
    res(0) += f;
    f = 0.5f*(vol1* q(1) + p1*dy + vol2* qinf[1] + p2*dy) + mu*(q(1)-qinf[1]);
    res(1) += f;
    f = 0.5f*(vol1* q(2) - p1*dx + vol2* qinf[2] - p2*dx) + mu*(q(2)-qinf[2]);
    res(2) += f;
    f = 0.5f*(vol1*(q(3)+p1)     + vol2*(qinf[3]+p2)    ) + mu*(q(3)-qinf[3]);
    res(3) += f;
  }
}
