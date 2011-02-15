
//
// atomic addition functions from Jon Cohen at NVIDIA
//

static __inline__ __device__ float atomicAdd(float *addr, float val)
{
  float old=*addr, assumed;

  do {
    assumed = old;
    old = int_as_float(__iAtomicCAS((int*)addr,
                       float_as_int(assumed),
                       float_as_int(val+assumed)));
  } while( assumed!=old );

  return old;
}


static __inline__ __device__ double atomicAdd(double *addr, double val)
{
  double old=*addr, assumed;

  do {
    assumed = old;
    old = __longlong_as_double( atomicCAS((unsigned long long int*)addr,
                                __double_as_longlong(assumed),
                                __double_as_longlong(val+assumed)));
  } while( assumed!=old );

  return old;
}


//
// my own min/max atomic functions, based on Jon Cohen's addition
//

static __inline__ __device__ float atomicMin(float *addr, float val)
{
  float old=*addr, assumed=CUDART_NAN_F;

  while (old>val && assumed!=old) {
    assumed = old;
    old = int_as_float(__iAtomicCAS((int*)addr,
                       float_as_int(assumed),
                       float_as_int(val)));
  }

  return old;
}


static __inline__ __device__ double atomicMin(double *addr, double val)
{
  double old=*addr, assumed=CUDART_NAN;

  while (old>val && assumed!=old) {
    assumed = old;
    old = __longlong_as_double( atomicCAS((unsigned long long int*)addr,
                                __double_as_longlong(assumed),
                                __double_as_longlong(val)));
  }

  return old;
}

static __inline__ __device__ float atomicMax(float *addr, float val)
{
  float old=*addr, assumed=CUDART_NAN_F;

  while (old<val && assumed!=old) {
    assumed = old;
    old = int_as_float(__iAtomicCAS((int*)addr,
                       float_as_int(assumed),
                       float_as_int(val)));
  }

  return old;
}


static __inline__ __device__ double atomicMax(double *addr, double val)
{
  double old=*addr, assumed=CUDART_NAN;

  while (old<val && assumed!=old) {
    assumed = old;
    old = __longlong_as_double( atomicCAS((unsigned long long int*)addr,
                                __double_as_longlong(assumed),
                                __double_as_longlong(val)));
  }

  return old;
}


//
// OP2 atomic operations for global reductions
//

template < class T >
__inline__ __device__ void op_atomic_inc(T *a, T *b, int dim) {

  for (int d=0; d<dim; d++) atomicAdd(&a[d],b[d]);
}

template < class T >
__inline__ __device__ void op_atomic_min(T *a, T *b, int dim) {

  for (int d=0; d<dim; d++) atomicMin(&a[d],b[d]);
}

template < class T >
__inline__ __device__ void op_atomic_max(T *a, T *b, int dim) {

  for (int d=0; d<dim; d++) atomicMax(&a[d],b[d]);
}

