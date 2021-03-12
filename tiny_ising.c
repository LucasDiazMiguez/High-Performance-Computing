/*
 * Tiny Ising model.
 * Loosely based on  "q-state Potts model metastability
 * study using optimized GPU-based Monte Carlo algorithms",
 * Ezequiel E. Ferrero, Juan Pablo De Francesco, Nicolás Wolovick,
 * Sergio A. Cannas
 * http://arxiv.org/abs/1101.0876
 *
 * Debugging: Ezequiel Ferrero
 */

#include "ising.h"
#include "params.h"

#include <stdlib.h> // rand()
#include <math.h> // expf()
#include <stdio.h> // printf()
#include <time.h> // time()
#include <assert.h>
#include <limits.h> // UINT_MAX
#include <omp.h> // omp_get_wtime()


// Internal definitions and functions
// out vector size, it is +1 since we reach TEMP_
#define NPOINTS (1+(int)((TEMP_MAX-TEMP_MIN)/DELTA_TEMP))
#define N (L*L) // system size
#define SEED (time(NULL)) // random seed


// temperature, E, E^2, E^4, M, M^2, M^4
struct statpoint {
	double t;
	double e; double e2; double e4;
	double m; double m2; double m4;
};

// The grid: global array
static int grid[L][L] = {{0}};


static
void
cycle(int grid[L][L],
      const double min, const double max,
      const double step, const unsigned int calc_step,
      struct statpoint stats[]) {

	assert((0<step && min<=max) || (step<0 && max<=min));
	int modifier = (0<step)?1:-1;

	unsigned int index = 0;
	for (double temp=min; modifier*temp<=modifier*max; temp+=step) {
		printf("Temp: %f\n", temp);

		// equilibrium phase
		for (unsigned int j=0; j<TRAN; ++j) {
			update(temp, grid);
		}

		// measurement phase
		unsigned int measurements = 0;
		double e=0.0, e2=0.0, e4=0.0, m=0.0, m2=0.0, m4=0.0;
		for (unsigned int j=0; j<TMAX; ++j) {
			update(temp, grid);
			if (j%calc_step==0) {
				double energy = 0.0, mag = 0.0;
				int M_max = 0;
				energy = calculate(grid, &M_max);
				mag = abs(M_max)/(float)N;
				e  += energy;
				e2 += energy*energy;
				e4 += energy*energy*energy*energy;
				m  += mag;
				m2 += mag*mag;
				m4 += mag*mag*mag*mag;
				++measurements;
			}
		}
		assert(index < NPOINTS);
		stats[index].t = temp;
		stats[index].e += e/measurements;
		stats[index].e2 += e2/measurements;
		stats[index].e4 += e4/measurements;
		stats[index].m += m/measurements;
		stats[index].m2 += m2/measurements;
		stats[index].m4 += m4/measurements;
		++index;
	}
}


static void init() {
	for (unsigned int i=0; i<L; ++i) {
		for (unsigned int j=0; j<L; ++j) {
			grid[i][j] = 1;
		}
	}
}


int
main(void)
{
	// the stats
	struct statpoint stat[NPOINTS];
	for (unsigned int i=0; i<NPOINTS; ++i) {
			stat[i].t = 0.0;
			stat[i].e = stat[i].e2 = stat[i].e4 = 0.0;
			stat[i].m = stat[i].m2 = stat[i].m4 = 0.0;
	}

	// parameter checking
	assert(TEMP_MIN<=TEMP_MAX);
	assert(0<DELTA_T && DELTA_T<TMAX); // at least one calculate()
	assert(TMAX%DELTA_T==0); // take equidistant calculate()
	assert((L*L/2)*4L<UINT_MAX); // max energy, that is all spins are the same, fits into a ulong

	// print header
	printf("# L: %i\n", L);
	printf("# Minimum Temperature: %f\n", TEMP_MIN);
	printf("# Maximum Temperature: %f\n", TEMP_MAX);
	printf("# Temperature Step: %.12f\n", DELTA_TEMP);
	printf("# Equilibration Time: %i\n", TRAN);
	printf("# Measurement Time: %i\n", TMAX);
	printf("# Data Acquiring Step: %i\n", DELTA_T);
	printf("# Number of Points: %i\n", NPOINTS);

	// configure RNG
	srand(SEED);

	// start timer
	double start = omp_get_wtime();

	// clear the grid
	init();

	// temperature increasing cycle
	cycle(grid, TEMP_MIN, TEMP_MAX, DELTA_TEMP, DELTA_T, stat);

	// stop timer
	double elapsed = omp_get_wtime()-start;
	printf("# Total Simulation Time (sec): %lf\n", elapsed);

	printf("# Temp\tE\tE^2\tE^4\tM\tM^2\tM^4\n");
	for (unsigned int i=0; i<NPOINTS; ++i) {
		printf ("%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",
			stat[i].t,
			stat[i].e/((double)N),
			stat[i].e2/((double)N*N),
			stat[i].e4/((double)N*N*N*N),
			stat[i].m,
			stat[i].m2,
			stat[i].m4);
	}

	return 0;
}
