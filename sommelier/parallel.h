#ifndef _PARALLEL_H_
#define _PARALLEL_H_

#include "smat.h"

#ifdef SESC
#include "sescapi.h"
#else
#include <pthread.h>
#endif

// a thread structure shared by
// all four operations
struct work_struct {
  const struct smat* amat; // input matrix A
  const struct smat* bmat; // input matrix B
  struct smat* result; // output matrix R
  int size; // size of the matrix
  int start_row; // starting row of a thread
  int end_row; // end row of a thread
  int id; // thread id
#ifndef SESC
  pthread_t pthread;
#endif /* !SESC */
};

#ifdef SESC

static inline void job_create(void (*func), void *arg, long flags)
{
  sesc_spawn(func, arg, flags);
}

static inline int job_join(struct work_struct *work)
{
  sesc_wait();
  return 0;
}

static inline void job_init(void)
{
  sesc_init();
}

static inline void job_exit(void)
{
  sesc_exit(0);
}

#else

static inline void job_create(void (*func), void *arg, long flags)
{
  struct work_struct *work = arg;

  pthread_create(&work->pthread, NULL, func, arg);
}

static inline int job_join(struct work_struct *work)
{
  return pthread_join(work->pthread, NULL);
}

static inline void job_init(void)
{
}

static inline void job_exit(void)
{
  pthread_exit(NULL);
}

#endif /* SESC */

#endif /* _PARALLEL_H_ */
