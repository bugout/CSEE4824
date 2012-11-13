/*
 * job_api_example.c
 * Simple and stupid example to illustrate how to use the tiny encapsulation
 * layer "job_api" around both pthread/sesc_api. Using this layer allows
 * us to write a single code file, that can then be compiled natively
 * and cross-compiled for SESC.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/*
 * From your implementation, you should include "job_api.h", and modify
 * the work_struct in that file. Leave the example one untouched.
 */
#include "job_api_example.h"

/*
 * NOTE: To achieve real parallelism on SESC, the number of strings here should
 * be equal to the number of cores in your SESC config file.
 */
#define NUM_STRINGS	2

static char *strings[] = {
	"Live by the foma that make you brave and kind and healthy and happy.",
	"One of my most productive days was throwing away 1000 lines of code.",
};

static void do_work(void *args)
{
	struct work_struct *s = args;
	size_t len = strlen(s->string);

	printf("Thread %d: the string '%s' is %Zu characters long.\n",
		s->id, s->string, len);

	/*
	 * The parent thread has id == 0. Since it is not a child thread,
	 * we don't call job_exit() for it.
	 */
	if (s->id)
		job_exit();
}

int main(int argc, char *argv[])
{
	struct work_struct *jobs[NUM_STRINGS];
	int i;

	for (i = 0; i < NUM_STRINGS; i++) {
		jobs[i] = malloc(sizeof(struct work_struct));
		if (jobs[i] == NULL) {
			perror(NULL);
			exit(1);
		}

		jobs[i]->string = strings[i];
		jobs[i]->id = i;
	}

	/* Initialise the job API */
	job_init();

	/*
	 * Spawn child threads.
	 * We leave the 'flags' argument empty--we could use it though
	 * to attach the spawned thread to a particular core (although this
	 * would only work on SESC). See sesacpi.h in case you really need
	 * to do this--it only makes sense if you are designing a multi-core
	 * CPU with cores of different sizes.
	 */
	for (i = 1; i < NUM_STRINGS; i++)
		job_create(do_work, jobs[i], 0);

	/*
	 * thread 0 is the parent thread. It doesn't need to be spawned, but
	 * it also has some work to do--do it now.
	 */
	do_work(jobs[0]);

	/* wait untill all the spawned processes have finished */
	for (i = 1; i < NUM_STRINGS; i++)
		job_join(jobs[i]);

	return 0;
}
