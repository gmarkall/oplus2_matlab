inline void update(const op_data<float> r, op_data<float> du, op_data<float> u, op_data<float> u_sum, op_data<float> u_max){
  *u += *du + alpha * (*r);
  *du = 0.0f;
  *u_sum += (*u)*(*u);
  *u_max = MAX(*u_max,*u);
}
