inline void res(const op_data<float> A, const op_data<float> u, op_data<float> du, const op_data<float> beta){
  *du += (*beta)*(*A)*(*u);
}
