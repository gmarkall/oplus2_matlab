inline void save_soln(op_data<float,4> q, op_data<float,4> qold){
  for (int n=0; n<4; n++) qold(n) = q(n);
}
