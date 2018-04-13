#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
int heat2dSolve(int M, int N, double eps, int print, double **u, double *tol)
{

	int iterations = 0;
	int iterations_print = 1;
	int i,j;
	double diff = 2.0 * eps;
	double *rowPrev; /* copy of the previous row in u */
	double *rowCurr; /* copy of the current row in in */
	double *rowTmp;

	rowPrev = calloc(N, sizeof(double));
	rowCurr = calloc(N, sizeof(double));
	if (print) 
		printf( "\n Iteration  Change\n" );

	while ( eps <= diff )
	{
		/*
			Initialize copy of "current" row 
		*/
		memcpy(rowCurr, u[0], N*sizeof(double));
		/*
		Determine the new estimate of the solution at the interior points.
		The new solution W is the average of north, south, east and west 
		neighbors.  */
		diff = 0.0;
		for ( i = 1; i < M - 1; i++ )
		{
			/* swap rowPrev and rowCurr pointers. Save the current row */
			rowTmp = rowPrev; rowPrev=rowCurr; rowCurr=rowTmp;
			memcpy(rowCurr, u[i], N*sizeof(double));
	
			for ( j = 1; j < N - 1; j++ )
			{
				u[i][j] = (rowPrev[j] + u[i+1][j] + 
					rowCurr[j-1] + rowCurr[j+1] ) / 4.0;
	
				double delta = fabs(rowCurr[j] - u[i][j]);
				if ( diff < delta ) 
				{
					diff = delta; 
				}
			}
		}
		iterations++;
		if ( print && iterations == iterations_print )
		{
			printf ( "  %8d  %f\n", iterations, diff );
			iterations_print *= 2;
		}
	} 
	/* memory cleanup */
	free(rowCurr);
	free(rowPrev);
	*tol = diff;
	return iterations;
}
