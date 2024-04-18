#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>

double total_time = 0;


void generate_random(double* input, size_t size)
{
	for (size_t i = 0; i < size;i++)
	{

		input[i] = rand() / (double)(RAND_MAX);
	}
}

double serial_sum(double* x, size_t size)
{
	double sum_val = 0.0;
#pragma omp parallel for
	{
		for (size_t i = 0; i < size;i++)
		{
			sum_val += x[i];
		}
	}
	return sum_val;
}



int main()
{
	size_t array_size = 100000;
	double* input = (double*)malloc(array_size * sizeof(double));
	int iterations = 1000;

	omp_set_num_threads(32);

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