typedef struct {
        double maxDiff;
        int counter;
        pthread_mutex_t *m;
        pthread_cond_t *c;
	int iters;
	double tol;
} shmem;

void heat2dSolve(long rank, long count, shmem *, int start_row, int end_row, int N, double eps, int print, double **u);
