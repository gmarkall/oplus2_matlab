void update(float *r, float *du, float *u){
  *u += *du + alpha * (*r);
  *du = 0.0f;
}
