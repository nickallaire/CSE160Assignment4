/*
Name: Nicholas Allaire
Email: nallaire@ucsd.edu
PID: A10639753
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "heat2d_solver.h"
/* Head2D Solver 
 *
 * Dirichlet boundary conditions 
 *
 */

/* heat2dSolve 
 * 	M - number of rows (input)
 * 	N - number of cols (output)
 * 	u - temperature distribution(input/output)
 *	eps - tolerance
 *	print - print iteration information (boolean)
 *
 * 	returns
 * 	    - number of iterations
 * 	    - u contains the final temperature distribution 
*/
void heat2dSolve(long rank, long count, shmem *shared, int start_row, int end_row, int N, double eps, int print, double **u) {

	int iterations = 0;
	int iterations_print = 1;
	int i,j;
	double diff = 2.0 * eps;
	long my_rank = (long) rank;
	long thread_count = (long) count;
	double *rowPrev; /* copy of the previous row in u */
	double *rowCurr; /* copy of the current row in u */
	double *rowTmp;
	double *ghost1; /* copy of overlapping row in threads */
	double *ghost2; /* copy of overlapping row in threads */
			
	shared -> maxDiff = diff;

	rowPrev = calloc(N, sizeof(double));
	rowCurr = calloc(N, sizeof(double));
	ghost1 = calloc(N, sizeof(double));
	ghost2 = calloc(N, sizeof(double));

	if (print && my_rank == 0) 
		printf( "\n Iteration  Change\n" );

	while ( eps <= shared -> maxDiff ) {
		/*
			Initialize copy of "current" row 
		*/
		memcpy(rowCurr, u[start_row - 1], N * sizeof(double));
		
		/*
 			Set ghost row boundary values
		*/
		ghost1[0] = u[1][0];
                ghost2[0] = u[1][0];
                ghost1[N - 1] = u[1][N - 1];
                ghost2[N - 1] = u[1][N - 1];
		
		/*
 			Copy over ghost rows
		*/
		if (my_rank == 0) {
			memcpy(ghost1, u[end_row + 1], N * sizeof(double));
		} else if (my_rank == (thread_count - 1)) {
			memcpy(ghost1, u[start_row - 1], N * sizeof(double));
		} else {
			memcpy(ghost1, u[start_row - 1], N * sizeof(double));
			memcpy(ghost2, u[end_row + 1], N * sizeof(double));
		}
	
		/*
			Wait until all threads catch up
		*/
		pthread_mutex_lock(shared -> m);
		(shared -> counter)++;
		if (shared -> counter == thread_count) {
			shared -> counter = 0;
			shared -> maxDiff = 0.0;
			pthread_cond_broadcast(shared -> c);
		} else {
			while (pthread_cond_wait(shared -> c, shared -> m) != 0);
		}
		pthread_mutex_unlock(shared -> m);
		/*
			Determine the new estimate of the solution at the interior points.
			The new solution W is the average of north, south, east and west 
			neighbors.
		*/
		diff = 0.0;
		for (i = start_row; i <= end_row; i++) {
			/*
 				swap rowPrev and rowCurr pointers. Save the current row
			*/
			rowTmp = rowPrev; rowPrev=rowCurr; rowCurr=rowTmp;
			memcpy(rowCurr, u[i], N * sizeof(double));
	
			for (j = 1; j < N - 1; j++) {
				/*
					Calculate new cell: if a boundary row, use ghost row value, else 
					calculate with 4 adjacent neighbors
				*/
				pthread_mutex_lock(shared -> m);
				if (my_rank == 0) {
					if (i == end_row) {
						u[i][j] = (rowPrev[j] + ghost1[j] + rowCurr[j - 1] + rowCurr[j + 1]) / 4.0;
					} else {
                                        	u[i][j] = (rowPrev[j] + u[i+1][j] + rowCurr[j - 1] + rowCurr[j + 1]) / 4.0;
					}
				} else if (my_rank == (thread_count -1)) {
					if (i == start_row) {
						u[i][j] = (ghost1[j] + u[i + 1][j] + rowCurr[j - 1] + rowCurr[j + 1]) / 4.0;
					} else {
						u[i][j] = (rowPrev[j] + u[i + 1][j] + rowCurr[j - 1] + rowCurr[j + 1]) / 4.0;
					}
				} else {
					if (i == start_row) {
                                                u[i][j] = (ghost1[j] + u[i + 1][j] + rowCurr[j - 1] + rowCurr[j + 1]) / 4.0;
					} else if (i == end_row) {
						u[i][j] = (rowPrev[j] + ghost2[j] + rowCurr[j - 1] + rowCurr[j + 1]) / 4.0;
					} else {
                                                u[i][j] = (rowPrev[j] + u[i + 1][j] + rowCurr[j - 1] + rowCurr[j + 1]) / 4.0;
					}
				}

				double delta = fabs(rowCurr[j] - u[i][j]);
				if (diff < delta) {
                                	diff = delta;
                                }
				pthread_mutex_unlock(shared -> m);
			}
		}
		
		iterations++;
		if (diff > shared -> maxDiff) {
			shared -> maxDiff = diff;
		}

		/*
 			Wait here until all threads catch up
		*/
		pthread_mutex_lock(shared -> m);
                (shared -> counter)++;
                if (shared -> counter == thread_count) {
                        shared -> counter = 0;
                        pthread_cond_broadcast(shared -> c);
                } else {
                        while (pthread_cond_wait(shared -> c, shared -> m) != 0);
                }
                pthread_mutex_unlock(shared -> m);

		/*
 			Print iterations message
		*/
		if (my_rank == 0 && print && iterations == iterations_print) {
			printf ( "  %8d  %f\n", iterations, shared -> maxDiff );
			iterations_print *= 2;
		}
	} 
	/* memory cleanup */
	free(rowCurr);
	free(rowPrev);
	free(ghost1);
	free(ghost2);

	if (my_rank == 0) {
		shared -> tol = shared -> maxDiff;
		shared -> iters = iterations;
	}
	//return iterations;
}
