void update(float *r, float *du, float *u){
  *u += *du + *r;
  *du = 0.0f;
}
