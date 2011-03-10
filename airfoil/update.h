inline void update(op_data<float,4> qold, op_data<float,4> q, op_data<float,4> res, op_data<float> adt, op_data<float> rms){
  float del, adti;

  adti = 1.0f/(*adt);

  for (int n=0; n<4; n++) {
    del    = adti*res(n);
    q(n)   = qold(n) - del;
    res(n) = 0.0f;
    *rms  += del*del;
  }
}
