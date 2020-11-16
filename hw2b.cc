#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#define PNG_NO_SETJMP
#define R_M 1 //row-major
#define C_M 2 //column-major
#define DEBUG_MODE 0
#define IS_TIMING 1
#define SHOW_LOAD 1
#define CORES_PER_CPU 1
#define CHUNKSIZE 1000

#include <sched.h>
#include <assert.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>

#include <iostream>
#include <iomanip>

/* HYBRID VERSION 1
 *
 * LET CPUs divide work load evenly
 *
 * */

int mode;
int num_cpu;
int SIZE;
int iters;
double left;
double right;
double lower;
double upper;
int width;
int height;
int* image;
int* T_image;

unsigned long long* Load_count;
unsigned long long* T_Load_count;
struct timespec begin,end;
double elapsed;

void write_png(const char* filename, int iters, int width, int height, const int* buffer) {
	FILE* fp = fopen(filename, "wb");
	assert(fp);
	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	assert(png_ptr);
	png_infop info_ptr = png_create_info_struct(png_ptr);
	assert(info_ptr);
	png_init_io(png_ptr, fp);
	png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
		     PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	png_set_filter(png_ptr, 0, PNG_NO_FILTERS);
	png_write_info(png_ptr, info_ptr);
	png_set_compression_level(png_ptr, 1);
	size_t row_size = 3 * width * sizeof(png_byte);
	png_bytep row = (png_bytep)malloc(row_size);
	for (int y = 0; y < height; ++y) {
		memset(row, 0, row_size);
		for (int x = 0; x < width; ++x) {
			int p = buffer[(height - 1 - y) * width + x];
			png_bytep color = row + x * 3;
			if (p != iters) {
				if (p & 16) {
					color[0] = 240;
					color[1] = color[2] = p % 16 * 16;
				} else {
					color[0] = p % 16 * 16;
				}
			}
		}
		png_write_row(png_ptr, row);
	}
	free(row);
	png_write_end(png_ptr, NULL);
	png_destroy_write_struct(&png_ptr, &info_ptr);
	fclose(fp);
}

int main(int argc, char** argv) {
	/* detect how many CPUs are available */
	cpu_set_t cpu_set;
	sched_getaffinity(0, sizeof(cpu_set), &cpu_set);
	num_cpu = CPU_COUNT(&cpu_set);
	if(DEBUG_MODE) std::cout << num_cpu << " cpus available" << std::endl;
	if(IS_TIMING) clock_gettime(CLOCK_MONOTONIC, &begin);

	/* argument parsing */
	assert(argc == 9);
	const char* filename = argv[1];
	iters = strtol(argv[2], 0, 10);
	left = strtod(argv[3], 0);
	right = strtod(argv[4], 0);
	lower = strtod(argv[5], 0);
	upper = strtod(argv[6], 0);
	width = strtol(argv[7], 0, 10);
	height = strtol(argv[8], 0, 10);
	unsigned long long total_Load = 0;

	mode = (height>width)?C_M:R_M;

	/* allocate memory for image */
	image = (int*)malloc(width * height * sizeof(int));
	assert(image);
	T_image = (int*)malloc(width * height * sizeof(int));
	assert(T_image);
	memset(image, 0, width * height *sizeof(int));
	memset(T_image, 0, width * height *sizeof(int));


	MPI_Init(&argc, &argv);
	int mpi_rank, mpi_ranks, omp_threads, omp_thread;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_ranks);
	Load_count = (unsigned long long*)malloc(mpi_ranks * num_cpu * sizeof(unsigned long long));
	T_Load_count = (unsigned long long*)malloc(mpi_ranks * num_cpu * sizeof(unsigned long long));
	for(int i = 0;i<mpi_ranks*num_cpu;i++){
		Load_count[i] = 0;
		T_Load_count[i] = 0;
	}
	if(DEBUG_MODE) printf("MPI rank %d\n",mpi_rank);

	/* mandelbrot set */

	omp_set_num_threads(num_cpu);
	if(DEBUG_MODE) printf("\n");
	#pragma omp parallel private(omp_thread)
	{
	omp_thread = omp_get_thread_num();
//	if(DEBUG_MODE) printf("thread_num %d of Rank %d\n",omp_thread, mpi_rank);
	#pragma omp for schedule(dynamic,1)
	for (int j = mpi_rank; j < height; j+=mpi_ranks) {
		if(DEBUG_MODE) printf("Row %d done by thread_num %d of Rank %d, id = %d\n", j, omp_thread, mpi_rank, mpi_rank*num_cpu + omp_thread);
		double y0 = j * ((upper - lower) / height) + lower;
		for (int i = 0; i < width; ++i) {
			double x0 = i * ((right - left) / width) + left;

			int repeats = 0;
			double x = 0;
			double y = 0;
			double length_squared = 0;
			while (repeats < iters && length_squared < 4) {
				double temp = x * x - y * y + x0;
				y = 2 * x * y + y0;
				x = temp;
				length_squared = x * x + y * y;
				++repeats;
				++Load_count[mpi_rank*num_cpu+omp_thread];
			}
			image[j * width + i] = repeats;
//			if(DEBUG_MODE) printf("image(%d,%d) :%d\n",j,i,repeats);
		}
	}
	}

	MPI_Reduce(image, T_image, width * height, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(Load_count, T_Load_count, mpi_ranks * num_cpu, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
	/* draw and cleanup */
	if(mpi_rank==0){
		write_png(filename, iters, width, height, T_image);
	}
	free(image);
	free(T_image);
	if(mpi_rank==0){
		if(IS_TIMING){
			clock_gettime(CLOCK_MONOTONIC, &end);
			elapsed = end.tv_sec - begin.tv_sec;
			elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
			std::cout << "Time used: " << elapsed << std::endl;
		}
		if(SHOW_LOAD){
			for(int i = 0;i<mpi_ranks * num_cpu;i++){
				std::cout << "Load[" <<i<<"]: "<< T_Load_count[i] << std::endl;
				total_Load += T_Load_count[i];
			}
			std::cout << "Total work load: " << total_Load << std::endl;
			for(int i = 0;i<mpi_ranks * num_cpu;i++){
				std::cout << "Thread " << i;
				std::cout << "runs " << std::fixed << std::setw(2) << std::setprecision(2) << (double(T_Load_count[i])/double(total_Load))*100 << "% Work Loads" << std::endl;
			}
		}
	}
}
