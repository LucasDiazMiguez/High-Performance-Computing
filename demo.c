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

#include "colormap.h"
#include "gl2d.h"
#include "ising.h"
#include "params.h"
#include "wtime.h"

#include <assert.h>
#include <limits.h> // UINT_MAX
#include <stdio.h> // printf()
#include <stdlib.h> // rand()
#include <string.h>
#include <time.h> // time()

#define MAXFPS 60
#define N (L * L) // system size
#define SEED (time(NULL)) // random seed


/**
 * GL output
 */
static void draw(gl2d_t gl2d, float t_now, float t_min, float t_max, int grid[L][L])
{
    static double last_frame = 0.0;

    double current_time = wtime();
    if (current_time - last_frame < 1.0 / MAXFPS) {
        return;
    }
    last_frame = current_time;

    float row[L * 3];
    float color[3];
    colormap_rgbf(COLORMAP_VIRIDIS, t_now, t_min, t_max, &color[0], &color[1], &color[2]);
    for (int i = 0; i < L; ++i) {
        memset(row, 0, sizeof(row));
        for (int j = 0; j < L; ++j) {
            if (grid[i][j] > 0) {
                row[j * 3] = color[0];
                row[j * 3 + 1] = color[1];
                row[j * 3 + 2] = color[2];
            }
        }
        gl2d_draw_rgbf(gl2d, 0, i, L, 1, row);
    }
    gl2d_display(gl2d);
}


static void cycle(gl2d_t gl2d, const double min, const double max, const double step, int grid[L][L])
{
    assert((0 < step && min <= max) || (step < 0 && max <= min));
    int modifier = (0 < step) ? 1 : -1;

    for (double temp = min; modifier * temp <= modifier * max; temp += step) {
        printf("Temp: %f\n", temp);
        for (unsigned int j = 0; j < TRAN + TMAX; ++j) {
            update(temp, grid);
            draw(gl2d, temp, min < max ? min : max, min < max ? max : min, grid);
        }
    }
}


static void init(int grid[L][L])
{
    for (unsigned int i = 0; i < L; ++i) {
        for (unsigned int j = 0; j < L; ++j) {
            grid[i][j] = (rand() / (float)RAND_MAX) < 0.5f ? -1 : 1;
        }
    }
}


int main(void)
{
    // parameter checking
    assert(TEMP_MIN <= TEMP_MAX);
    assert(0 < DELTA_T && DELTA_T < TMAX); // at least one calculate()
    assert(TMAX % DELTA_T == 0); // take equidistant calculate()
    assert((L * L / 2) * 4L < UINT_MAX); // max energy, that is all spins are the same, fits into a ulong

    // print header
    printf("# L: %i\n", L);
    printf("# Minimum Temperature: %f\n", TEMP_MIN);
    printf("# Maximum Temperature: %f\n", TEMP_MAX);
    printf("# Temperature Step: %.12f\n", DELTA_TEMP);
    printf("# Equilibration Time: %i\n", TRAN);
    printf("# Measurement Time: %i\n", TMAX);

    // configure RNG
    srand(SEED);

    gl2d_t gl2d = gl2d_init("tiny_ising", L, L);

    // start timer
    double start = wtime();

    // clear the grid
	int grid[L][L] = { { 0 } };
    init(grid);

    // temperature increasing cycle
    cycle(gl2d, TEMP_MIN, TEMP_MAX, DELTA_TEMP, grid);

    // stop timer
    double elapsed = wtime() - start;
    printf("# Total Simulation Time (sec): %lf\n", elapsed);

    gl2d_destroy(gl2d);

    return 0;
}
