/*
 * Naive implementations of the matrix algorithms.
 * Note that the result matrix r is not initisialised to 0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "smat.h"
#include "parallel.h"

#define LOOP_UNROLLING
#define UNROLLING_SIZE 1
#define NUM_THREADS 2
#define BLOCK_SIZE 400

inline int min(int a, int b) {
  return a < b ? a: b;
}

void init_thread_args(struct work_struct *job, int id, const struct smat* a, 
		      const struct smat* b, struct smat* r, int start, int end)
{
  job->id = id;
  job->amat = a;
  job->bmat = b;
  job->result = r;
  job->start_row = start;
  job->end_row = end;
  job->size = a->rows;
}


void cleanup(struct work_struct* jobs[]) {
  int i;
  for (i = 0; i < NUM_THREADS; i++)
    free(jobs[i]);
}

void smat_mult_t(void* args) {
  struct work_struct *s = args;

  double** adata = s->amat->data;
  double** bdata = s->bmat->data;
  double** rdata = s->result->data;

  int size = s->size;

  double t;
  int i, j, k, jj, kk;

  for (i = s->start_row; i < s->end_row; i++)
    memset(rdata[i], 0, size * sizeof(double));

  
  for (jj = 0; jj < size; jj += BLOCK_SIZE)
    for (kk = 0; kk < size; kk += BLOCK_SIZE)
      for (i = s->start_row; i < s->end_row; i++) {
	double *a = adata[i];
	int jend = min(jj+BLOCK_SIZE, size);
	int kend = min(kk+BLOCK_SIZE, size);
	for (j = jj; j < jend; j++) {
	  double *b = bdata[j];
	  double *r = rdata[i];

	  t = 0;
	  int offset = (kend - kk) & (UNROLLING_SIZE - 1);    
	  for (k = kk; k < kk + offset; k++)
	    t = t + a[k] * b[k];

	  for (k = kk  + offset; k < kend; k+=UNROLLING_SIZE) {
	    t = t + a[k] * b[k];

#ifdef LOOP_UNROLLING_2
	    t += (a[k+1] * b[k+1]);
#endif

#ifdef LOOP_UNROLLING_4
	    t += (a[k+1] * b[k+1]);
	    t += (a[k+2] * b[k+2]);
	    t += (a[k+3] * b[k+3]);
#endif
	  }
	  r[j] = r[j] + t;
	}
      }

  if (s->id)
    job_exit();
}


/* Matrix multiplication, r = a x b */
void smat_mult(const struct smat *a, const struct smat *b, struct smat *r)
{
  struct work_struct *jobs[NUM_THREADS];
  int i = 0;

  int size = a->rows;

  // transpose
  struct smat transpose;
  transpose.rows = size;
  transpose.cols = size;
  transpose.data = calloc_2d_double(size, size);
  
  int j;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      transpose.data[j][i] = b->data[i][j];
	
  int nrows = a->rows / NUM_THREADS;

  for (i = 0; i < NUM_THREADS; i++) {
    jobs[i] = malloc(sizeof(struct work_struct));
    if (jobs[i] == NULL) {
      perror(NULL);
      exit(1);
    }
    int start, end;
    start = i * nrows;
    end = (i==NUM_THREADS-1) ? a->rows : start + nrows;

    init_thread_args(jobs[i], i, a, &transpose, r, start, end);
  }

  job_init();

  for (i = 1; i < NUM_THREADS; i++)
    job_create(smat_mult_t, jobs[i], 0);

  smat_mult_t(jobs[0]);

  for (i = 1; i < NUM_THREADS; i++)
    job_join(jobs[i]);

  cleanup(jobs);

  return;
}


void smat_vect_t(void* args) {
  struct work_struct *s = args;

  int i, j;
  double** adata = s->amat->data;
  double* vdata = s->bmat->data[0];
  double* rdata = s->result->data[0];

  int size = s->size;

  for (i = s->start_row; i < s->end_row; i++) {
    double t = 0;
    double *a = adata[i];

    int offset = s->size & (UNROLLING_SIZE - 1);

    for (j = 0; j < offset; j++)
      t += (a[j] * vdata[j]);

    for (j = offset; j < size; j+=UNROLLING_SIZE) {

      t += (a[j] * vdata[j]);
#ifdef LOOP_UNROLLING_2
      t += (a[j+1] * vdata[j+1]);
#endif

#ifdef LOOP_UNROLLING_4
      t += (a[j+1] * vdata[j+1]);
      t += (a[j+2] * vdata[j+2]);
      t += (a[j+3] * vdata[j+3]);
#endif

    }

    rdata[i] = t;
  }

  if (s->id)
    job_exit();
}




/*
 * Matrix-vector multiplication, r = a x v
 * NOTE: a is the matrix, v is the vector. This means v->cols is always 1.
 * The result is always a vector (r->cols == 1).
 */
void smat_vect(const struct smat *a, const struct smat *v, struct smat *r)
{
  struct work_struct *jobs[NUM_THREADS];
  int i = 0;
	
  int nrows = a->rows / NUM_THREADS;

  for (i = 0; i < NUM_THREADS; i++) {
    jobs[i] = malloc(sizeof(struct work_struct));
    if (jobs[i] == NULL) {
      perror(NULL);
      exit(1);
    }
    int start, end;
    start = i * nrows;
    end = (i==NUM_THREADS-1) ? a->rows : start + nrows;

    init_thread_args(jobs[i], i, a, v, r, start, end);
  }

  job_init();

  for (i = 1; i < NUM_THREADS; i++)
    job_create(smat_vect_t, jobs[i], 0);

  smat_vect_t(jobs[0]);

  for (i = 1; i < NUM_THREADS; i++)
    job_join(jobs[i]);

  cleanup(jobs);

  return;
}

void smat_add_t(void* args) {
  struct work_struct *s = args;

  int i, j;
  double** adata = s->amat->data;
  double** bdata = s->bmat->data;
  double** rdata = s->result->data;

  for (i = s->start_row; i < s->end_row; i++) {
    double *r = rdata[i];
    double *a = adata[i];
    double *b = bdata[i];

    int offset = s->size & (UNROLLING_SIZE - 1);

    for (j = 0; j < offset; j++)
      r[j] = a[j] + b[j];

    for (j = offset; j < s->size; j+=UNROLLING_SIZE) {

      r[j] = a[j] + b[j];
#ifdef LOOP_UNROLLING_2
      r[j+1] = a[j+1] + b[j+1];
#endif

#ifdef LOOP_UNROLLING_4
      r[j+1] = a[j+1] + b[j+1];
      r[j+2] = a[j+2] + b[j+2];
      r[j+3] = a[j+3] + b[j+3];
#endif
    }
      
  }

  if (s->id)
    job_exit();
}


/* Matrix addition, i.e. r = a + b */
void smat_add(const struct smat *a, const struct smat *b, struct smat *r)
{
  struct work_struct *jobs[NUM_THREADS];
  int i = 0;
	
  int nrows = a->rows / NUM_THREADS;

  for (i = 0; i < NUM_THREADS; i++) {
    jobs[i] = malloc(sizeof(struct work_struct));
    if (jobs[i] == NULL) {
      perror(NULL);
      exit(1);
    }
    int start, end;
    start = i * nrows;
    end = (i==NUM_THREADS-1) ? a->rows : start + nrows;

    init_thread_args(jobs[i], i, a, b, r, start, end);
  }

  job_init();

  for (i = 1; i < NUM_THREADS; i++)
    job_create(smat_add_t, jobs[i], 0);

  smat_add_t(jobs[0]);

  for (i = 1; i < NUM_THREADS; i++)
    job_join(jobs[i]);


  cleanup(jobs);

  return;
}


void smat_scale_t(void* args) {
  struct work_struct *s = args;
  
  int i, j;
  double** adata = s->amat->data;
  double factor = s->bmat->data[0][0];
  double** rdata = s->result->data;

  for (i = s->start_row; i < s->end_row; i++) {
    double *r = rdata[i];
    double *a = adata[i];

    int offset = s->size & (UNROLLING_SIZE - 1);    

    // process [0, offset)
    for (j = 0; j < offset; j++)
      r[j] = a[j] * factor;

    // process UNROLLING_SIZE elements at a time
    for (j = offset; j < s->size; j += UNROLLING_SIZE) {      
      r[j] = a[j] * factor;
#ifdef LOOP_UNROLLING_2
      r[j+1] = a[j+1] * factor;
#endif

#ifdef LOOP_UNROLLING_4
      r[j+1] = a[j+1] * factor;
      r[j+2] = a[j+2] * factor;
      r[j+3] = a[j+3] * factor;
#endif
    }
  }

  if (s->id)
    job_exit();
}


/* Scale matrix a by constant alpha */
void smat_scale(const struct smat *a, const struct smat *alpha, struct smat *r)
{
  struct work_struct *jobs[NUM_THREADS];

  int i;
  int nrows = a->rows / NUM_THREADS;

  for (i = 0; i < NUM_THREADS; i++) {
    jobs[i] = malloc(sizeof(struct work_struct));
    if (jobs[i] == NULL) {
      perror(NULL);
      exit(1);
    }
    int start, end;
    start = i * nrows;
    end = (i==NUM_THREADS-1) ? a->rows : start + nrows;

    init_thread_args(jobs[i], i, a, alpha, r, start, end);
  }

  job_init();

  for (i = 1; i < NUM_THREADS; i++)
    job_create(smat_scale_t, jobs[i], 0);

  smat_scale_t(jobs[0]);

  for (i = 1; i < NUM_THREADS; i++)
    job_join(jobs[i]);


  cleanup(jobs);

  return;
}
