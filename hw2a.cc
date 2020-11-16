#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#define PNG_NO_SETJMP
#define R_M 1 //row-major
#define C_M 2 //column-major
#define DEBUG_MODE 1
#define SHOW_LOAD 1
#define IS_TIMING 1
#define CORES_PER_CPU 1

#include <sched.h>
#include <pthread.h>
#include <assert.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <iomanip>

/* Pthread version1
 * tries to distribute pixels in min(row,col)-major
 * then do static static assigning
 *
 * Maybe just do loop +size distribute
 *
 * */


typedef struct thread_Info{
	int id;
}T_I;

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
unsigned long long* Load_count;
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

void* Thread_mandelbrot(void* arg){
	T_I* data = (T_I*) arg;
	int id = data->id;
	if(DEBUG_MODE) std::cout << "Thread " << id << "ran by CPU " << sched_getcpu() << std::endl;
/*	if(DEBUG_MODE) printf("Thread id: %d\n",id);
	if(DEBUG_MODE){
		for(int i=0;i<1000;){
			i++;//DO NOTHING
		}
	}
	if(DEBUG_MODE) printf("Mode: %d\n",mode);
*/
	unsigned long long pixel_count;
//	if(mode == R_M){
		for (int j = id; j < height; j+=SIZE) {
			double y0 = j * ((upper - lower) / height) + lower;
			for (int i = 0; i < width; ++i) {
				double x0 = i * ((right - left) / width) + left;
				pixel_count++;

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
					++Load_count[id];
				}
				image[j * width + i] = repeats;
			}
		}
//	}
/*	else{
		for (int i = id; i < width; i+=num_cpu) {
			double x0 = i * ((right - left) / width) + left;
			for (int j = 0; j < height; ++j) {
				double y0 = j * ((upper - lower) / height) + lower;
				pixel_count++;

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
					++Load_count;
				}
				image[j * width + i] = repeats;
			}
		}
	}
*/
	if(DEBUG_MODE) std::cout << Load_count[id] << " work load in thread" << id << std::endl;
//	if(DEBUG_MODE) printf("%llu pixels in cpu %d\n", pixel_count, id);
	pthread_exit(NULL);
}
/*
void master_mandelbrot(){
	unsigned long long master_pixel_count = 0;
//	if(mode == R_M){
		for (int j = 0; j < height; j+=SIZE) {
			double y0 = j * ((upper - lower) / height) + lower;
			for (int i = 0; i < width; ++i) {
				double x0 = i * ((right - left) / width) + left;
				master_pixel_count++;

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
					master_Load_count[0];
				}
				image[j * width + i] = repeats;
			}
		}
//	}
	else{
		for (int i = 0; i < width; i+=num_cpu) {
			double x0 = i * ((right - left) / width) + left;
			for (int j = 0; j < height; ++j) {
				double y0 = j * ((upper - lower) / height) + lower;
				master_pixel_count++;

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
					++master_Load_count;
				}
				image[j * width + i] = repeats;
			}
		}
	}
	if(DEBUG_MODE) printf("%llu work load in master\n", master_Load_count);
//	if(DEBUG_MODE) printf("%llu pixels in master\n", master_pixel_count);
}
*/
int main(int argc, char** argv) {
	/* detect how many CPUs are available */
	cpu_set_t cpu_set;
	sched_getaffinity(0, sizeof(cpu_set), &cpu_set);
	num_cpu = CPU_COUNT(&cpu_set);
	if(DEBUG_MODE) std::cout << num_cpu << "cpus available" << std::endl;
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

	mode = (height>width)?C_M:R_M;

	/* allocate memory for image */
	image = (int*)malloc(width * height * sizeof(int));
	Load_count = (unsigned long long*)malloc(num_cpu * sizeof(unsigned long long));
	for(int i = 0;i<SIZE;i++){
		Load_count[i] = 0;
	}
	assert(image);

	/* mandelbrot set */
	/*

	for (int j = 0; j < height; ++j) {
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
			}
			image[j * width + i] = repeats;
		}
	}

	*/
		
	/* pthread mandelbrot set*/
	SIZE = (CORES_PER_CPU)*num_cpu;
	T_I DATA[SIZE];
	pthread_t threads[SIZE];
	for(int i = 0;i<SIZE;i++){
		DATA[i].id = i;
		pthread_create(&threads[i], NULL, Thread_mandelbrot, (void*) &DATA[i]);
	}

//	master_mandelbrot();
	unsigned long long total_Load = 0;

	for(int i =0;i<SIZE;i++){
		pthread_join(threads[i], NULL);
		total_Load+=Load_count[i];
	}

	if(SHOW_LOAD) std::cout << "total_Load: " << total_Load << std::endl;

	if(SHOW_LOAD){
	for(int i =0;i<SIZE;i++){
		std::cout << "Thread " << i;
		std::cout << "runs " << std::fixed << std::setw(2) << std::setprecision(2) << (double(Load_count[i])/double(total_Load))*100 << "% Work Loads" << std::endl;
	}
	}

	/* draw and cleanup */
	write_png(filename, iters, width, height, image);
	free(image);
	if(IS_TIMING){
		clock_gettime(CLOCK_MONOTONIC, &end);
		elapsed = end.tv_sec - begin.tv_sec;
		elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
		std::cout << "Time used: " << elapsed << std::endl;
	}
}
