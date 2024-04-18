#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>

#define NUM_THREADS 64

double total_time = 0;

typedef struct {
	double val;
	char pad[128];
} sum_structure;

void generate_random(double* input, size_t size)
{
	for (size_t i = 0; i < size;i++)
	{
		input[i] = rand() / (double)(RAND_MAX);
	}
}

double serial_sum(double* x, size_t size)
{
	sum_structure sum_val[NUM_THREADS];
	double final_sum = 0;
#pragma omp parallel shared(sum_val)
	{
		int id = omp_get_thread_num();
#pragma omp  for
		for (size_t i = 0; i < size;i++)
		{
			sum_val[id].val += x[i];
		}
	}

	for (int i = 0;i < NUM_THREADS;i++)
	{
		final_sum += sum_val[i].val;
	}

	return final_sum;
}



int main()
{
	size_t array_size = 100000;
	double* input = (double*)malloc(array_size * sizeof(double));
	int iterations = 1000;

	omp_set_num_threads(NUM_THREADS);

	generate_random(input, array_size);

	for (int i = 0; i < iterations; i++)
	{
		clock_t start = clock();

		serial_sum(input, array_size);

		clock_t finish = clock();
		total_time += finish - start;
	}

	double average_time = (total_time / iterations) / CLOCKS_PER_SEC;
	printf("\nAVERAGE TIME: %.9f\n\n", average_time);

	double squaredSum = 0.0;
	for (int i = 0; i < iterations; ++i) {

		clock_t start = clock();

		// Call the function to compute the sum
		serial_sum(input, array_size);

		clock_t finish = clock();

		double elapsedTime = ((double)(finish - start)) / CLOCKS_PER_SEC; // Convert to milliseconds
		squaredSum += (elapsedTime - average_time) * (elapsedTime - average_time);
	}
	double stdDeviation = sqrt(squaredSum / iterations);
	printf("\nSTANDARD DEVIATION: %.9f\n\n", stdDeviation);


	return 0;
}