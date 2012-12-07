/*
 * Naive implementations of the matrix algorithms.
 * Note that the result matrix r is not initisialised to 0.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "smat.h"
#include "parallel.h"

#define NUM_THREADS 2
#define BLOCK_SIZE 32

inline int min(int a, int b) {
  return a < b ? a: b;
}

/* Matrix multiplication, r = a x b */
void smat_mult(const struct smat *a, const struct smat *b, struct smat *r)
{
  int jj, kk, i, j, k;
  double** adata = a->data;
  double** bdata = b->data;
  double** rdata = r->data;

  int size = a->rows;
  
  for (i = 0; i < size; i++)
    memset(rdata[i], 0, size * sizeof(double));

  
  for (jj = 0; jj < size; jj += BLOCK_SIZE)
    for (kk = 0; kk < size; kk += BLOCK_SIZE)
      for (i = 0; i < size; i++)
	for (j = jj; j < min(jj+BLOCK_SIZE,size); j++) {
	  double t = 0;
	  for (k = kk; k < min(kk+BLOCK_SIZE,size); k++)
	    t = t + adata[i][k] * bdata[k][j];
	  rdata[i][j] = rdata[i][j] + t;
	}
}


void smat_vect_t(void* args) {
  struct work_struct *s = args;

  int i, j;
  double** adata = s->matrix;
  double* vdata = s->matrixb[0];
  double* rdata = s->result[0];

  memset(rdata, 0, sizeof(double) * size);

  for (i = 0; i < s->end_row; i++) {
    double t = 0;
    for (j = 0; j < size; j++)
      t += (adata[i][j] * vdata[j]);
    rdata[i] += t;
  }


  int i, j;
  const struct smat* a = s->matrix;
  const struct smat* b = s->matrixb;
  struct smat* r = s->result;

  for (i = s->start_row; i < s->end_row; i++)
    for (j = 0; j < a->cols; j++)
      r->data[i][j] = a->data[i][j] + b->data[i][j];

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
  
      
}

void smat_add_t(void* args) {
  struct work_struct *s = args;

  int i, j;
  const struct smat* a = s->matrix;
  const struct smat* b = s->matrixb;
  struct smat* r = s->result;

  for (i = s->start_row; i < s->end_row; i++)
    for (j = 0; j < a->cols; j++)
      r->data[i][j] = a->data[i][j] + b->data[i][j];

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
    jobs[i]->id = i;
    jobs[i]->matrix = a;
    jobs[i]->matrixb = b;
    jobs[i]->result = r;
    jobs[i]->start_row = i * nrows;
    jobs[i]->end_row = (i+1) * nrows;
    if (i == NUM_THREADS - 1)
      jobs[i]->end_row = a->rows;
  }

  job_init();

  for (i = 1; i < NUM_THREADS; i++)
    job_create(smat_add_t, jobs[i], 0);

  smat_add_t(jobs[0]);

  for (i = 1; i < NUM_THREADS; i++)
    job_join(jobs[i]);

  return;
}




void smat_scale_t(void* args) {
  struct work_struct *s = args;

  int i, j;
  double factor = s->factor;
  const struct smat *a = s->matrix;
  const struct smat *r = s->result;

  for (i = s->start_row; i < s->end_row; i++)
    for (j = 0; j < a->cols; j++)
      r->data[i][j] = a->data[i][j] * factor;

  if (s->id)
    job_exit();
}


/* Scale matrix a by constant alpha */
void smat_scale(const struct smat *a, const struct smat *alpha, struct smat *r)
{
  double factor = alpha->data[0][0];

  struct work_struct *jobs[NUM_THREADS];
  int i = 0;
	
  int nrows = a->rows / NUM_THREADS;

  for (i = 0; i < NUM_THREADS; i++) {
    jobs[i] = malloc(sizeof(struct work_struct));
    if (jobs[i] == NULL) {
      perror(NULL);
      exit(1);
    }
    jobs[i]->id = i;
    jobs[i]->matrix = a;
    jobs[i]->result = r;
    jobs[i]->factor = factor;
    jobs[i]->start_row = i * nrows;
    jobs[i]->end_row = (i+1) * nrows;
    if (i == NUM_THREADS - 1)
      jobs[i]->end_row = a->rows;
  }

  job_init();

  for (i = 1; i < NUM_THREADS; i++)
    job_create(smat_scale_t, jobs[i], 0);

  smat_scale_t(jobs[0]);

  for (i = 1; i < NUM_THREADS; i++)
    job_join(jobs[i]);

  return;

}
