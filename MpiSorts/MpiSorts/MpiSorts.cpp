#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include<string.h>

int MPI_Is_sorted(int n, double* array, int* isSorted, int root, MPI_Comm comm);
int MPI_Bitonic_sort(int n, double* a, int root, MPI_Comm comm);
int MPI_Shell_sort(double* a, int n, int root, MPI_Comm comm);
int MPI_Sort_bucket(int n, double* a, double max, int root, MPI_Comm comm);
int MPI_Sort_oddEven(int n, double* array, int root, MPI_Comm comm);
int MPI_Exchange(int n, double* array, int rank1, int rank2, MPI_Comm comm);
int MPI_Direct_Sort(int n, double* array, int root, MPI_Comm comm);

void exchange(double *a, int i, int j);
void compare(double* a, int i, int j, int dir);
void merge_sort(int n, double* a);
void swap(double* a, double* b);
void sort(int n, double* a);
double* merge_array(int n, double* a, int m, double* b);
double* merge(int n, double* array, int m, double* b);


int main(int argc, char* argv[])
{
	int rank, size, i, n, * recvbuf, * bitonicRecvbuf;
	double* a;
	int isSorted;
	double m = 10.0;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	n = 1000000;

	a = (double*)calloc(n, sizeof(double));


	if (rank == 0)
	{
		srand(((unsigned)time(NULL) + rank));

		for (i = 0; i < n; i++)
		{
			a[i] = ((double)rand() / RAND_MAX) * m;
	
		}
	}
	//copy of the vector

	double* a_copy = (double*)calloc(n, sizeof(double));
	memcpy(a_copy, a, n * sizeof(double));

	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);
	//SHELL SORT

	printf("1.\n");
	start_time = MPI_Wtime();
	MPI_Shell_sort(a,n,0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		end_time = MPI_Wtime();

		MPI_Is_sorted(n / size, a_copy, &isSorted, 0, MPI_COMM_WORLD);
		if (isSorted)
		{
			printf("Shell sort: Array is sorted\n");
		}
		else
		{
			printf("Shell sort: Array is NOT sorted\n");
		}

		printf("Shell sort: Time taken = %lf seconds\n", end_time - start_time);

	}

	//BUCKET SORT
	printf("2.\n");
	memcpy(a_copy, a, n * sizeof(double));
	start_time = MPI_Wtime();
	MPI_Sort_bucket(n, a_copy, m, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		end_time = MPI_Wtime();
		MPI_Is_sorted(n / size, a_copy, &isSorted, 0, MPI_COMM_WORLD);
		if (isSorted)
		{
			printf("Bucket sort: Array is sorted\n");
		}
		else
		{
			printf("Bucket sort: Array is NOT sorted\n");
		}

		printf("Bucket sort: Time taken = %lf seconds\n", end_time - start_time);
	}

	//ODD EVEN SORT
	printf("3.\n");
	start_time = MPI_Wtime();
	memcpy(a_copy, a, n * sizeof(double));
	MPI_Sort_oddEven(n, a_copy, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		end_time = MPI_Wtime();
		MPI_Is_sorted(n / size, a_copy, &isSorted, 0, MPI_COMM_WORLD);
		if (isSorted)
		{
			printf("Odd even sort: Array is sorted\n");
		}
		else
		{
			printf("Odd even sort: Array is NOT sorted\n");
		}

		printf("Odd even sort: Time taken = %lf seconds\n", end_time - start_time);

	}

	//DIRECT SORT
	//memcpy(a_copy, a, n * sizeof(double));
	//printf("4.\n");
	//start_time = MPI_Wtime();
	//MPI_Direct_Sort(n, a, 0, MPI_COMM_WORLD);
	//if (rank == 0)
	//{
	//	end_time = MPI_Wtime();
	//	MPI_Is_sorted(n / size, a_copy, &isSorted, 0, MPI_COMM_WORLD);
	//	if (isSorted)
	//	{
	//		printf("Direct sort: Array is sorted\n");
	//	}
	//	else
	//	{
	//		printf("Direct sort: Array is NOT sorted\n");
	//	}

	//	printf("Direct sort: Time taken = %lf seconds\n", end_time - start_time);

	//}
	//free(a_copy);

	//BITONIC SORT
	printf("5.\n");
	start_time = MPI_Wtime();
	MPI_Bitonic_sort(n,a,0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		end_time = MPI_Wtime();
		MPI_Is_sorted(n / size, a, &isSorted, 0, MPI_COMM_WORLD);
		if (isSorted)
		{
			printf("Bitonic sort: Array is sorted\n");
		}
		else
		{
			printf("Bitonic sort: Array is NOT sorted\n");
		}

		printf("Bitonic sort: Time taken = %lf seconds\n", end_time - start_time);

		free(a);
	}


	MPI_Finalize();
	free(a_copy);
	return 0;
}
int MPI_Shell_sort(double* a, int n, int root, MPI_Comm comm)
{
	int rank, size;
	int gap_size = 8;
	int gap_array[8] = { 1, 4, 10, 23, 57, 132, 301, 701 };

	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	int local_n = n / size;
	double* local_array = (double*)malloc(local_n * sizeof(double));

	MPI_Scatter(a, local_n, MPI_DOUBLE, local_array, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);



	for (int i = gap_size - 1; i >= 0; i--) {
		int gap = gap_array[i];
		for (int j = gap; j < n; j++) {
			double temp = a[j];
			int k;
			for (k = j; k >= gap && a[k - gap] > temp; k -= gap) {
				a[k] = a[k - gap];
			}
			a[k] = temp;
		}
		MPI_Allreduce(MPI_IN_PLACE, a, n, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	}

	MPI_Gather(local_array, local_n, MPI_DOUBLE, a, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	return MPI_SUCCESS;
}

void exchange(double* a, int i, int j)
{
	int t = a[i];
	a[i] = a[j];
	a[j] = t;
}
void compare(double* a, int i, int j, int dir)
{
	if (dir == (a[i] > a[j]))
		exchange(a, i, j);
}
void bitonicMerge(double* a, int low, int cnt, int dir)
{
	if (cnt > 1)
	{
		int k = cnt / 2;
		for (int i = low; i < low + k; i++)
			compare(a, i, i + k, dir);
		bitonicMerge(a, low, k, dir);
		bitonicMerge(a, low + k, k, dir);
	}
}
void bitonicSort(double* a, int low, int cnt, int dir)
{
	if (cnt > 1)
	{
		int k = cnt / 2;
		bitonicMerge(a, low, k, dir);
		bitonicMerge(a, low + k, k, dir);
		bitonicSort(a, low, k, dir);
		bitonicSort(a, low + k, k, dir);
	}
}

int MPI_Bitonic_sort(int n, double* a, int root, MPI_Comm comm)
{
	int rank, size, local_n;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	local_n = n / size;

	double* local_a = (double*)malloc(local_n * sizeof(double));

	MPI_Scatter(a, local_n, MPI_DOUBLE, local_a, local_n, MPI_DOUBLE, root, MPI_COMM_WORLD);
	bitonicSort(local_a, 0, local_n, 1);
	MPI_Gather(local_a, local_n, MPI_DOUBLE, a, local_n, MPI_DOUBLE, root, MPI_COMM_WORLD);

	free(local_a);
	return MPI_SUCCESS;
}

int isSorted(double* arr, int n) {
	int i;
	for (i = 0; i < n - 1; i++) {
		if (arr[i] > arr[i + 1]) {
			return 0;
		}
	}
	return 1;
}

int MPI_Sort_bucket(int n, double* a, double max, int root, MPI_Comm comm)
{
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// allocate the extra memory / arrays needed
	double* bucket = (double*)calloc(n, sizeof(double));

	// Brodcast the array to the processor
	int error = MPI_Bcast(a, n, MPI_DOUBLE, root, comm);
	if (error != MPI_SUCCESS)
	{
		return error;
	}

	//scan a to collect the elements of bucket rank
	int count = 0;
	for (int i = 0; i < n; i++) {
		if (rank * max / size <= a[i] && a[i] < (rank + 1) * max / size) {
			bucket[count++] = a[i];
		}
	}

	//sort bucket
	merge_sort(count, bucket);

	//gather coutn to counts
	int* counts = (int*)calloc(size, sizeof(int));
	error = MPI_Gather(&count, 1, MPI_INT, counts, 1, MPI_INT, root, comm);
	if (error != MPI_SUCCESS)
	{
		return error;
	}

	int* displs = (int*)calloc(size, sizeof(int));
	displs[0] = 0;

	for (int i = 1; i < size; i++) {
		displs[i] = displs[i - 1] + counts[i - 1];
	}

	error = MPI_Gatherv(bucket, count, MPI_DOUBLE, a, counts, displs, MPI_DOUBLE, root, comm);
	if (error != MPI_SUCCESS)
	{
		return error;
	}

	return MPI_SUCCESS;
}

int MPI_Is_sorted(int n, double* array, int* isSorted, int root, MPI_Comm comm)
{
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double* first = (double*)calloc(size, sizeof(double));
	double* last = (double*)calloc(size, sizeof(double));
	MPI_Gather(&array[0], 1, MPI_DOUBLE, first, 1, MPI_DOUBLE, root, comm);
	MPI_Gather(&array[n - 1], 1, MPI_DOUBLE, last, 1, MPI_DOUBLE, root, comm);
	if (rank == root)
	{
		*isSorted = 1;
		for (int i = 0; i < size - 1; i++)
			if (last[i] > first[i + 1])
			{
				*isSorted = 0;
				break;
			}

	}

	MPI_Bcast(isSorted, 1, MPI_INT, root, comm);
	free(first);
	free(last);
	return MPI_SUCCESS;
}
int MPI_Sort_oddEven(int n, double* array, int root, MPI_Comm comm) {

	// get rank and size of comm
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//allocate space for numElements/numProcessors amount of doubles
	double* local_a = (double*)calloc(n / size, sizeof(double));

	//scatter a to local_a
	int error = MPI_Scatter(array, n / size, MPI_DOUBLE, local_a, n / size, MPI_DOUBLE, root, comm);
	if (error != MPI_SUCCESS)
	{
		free(local_a);
		return error;
	}

	//sort local_a using mergeSort
	merge_sort(n / size, local_a);

	//odd-even iterations
	for (int i = 0; i < size; i++)
	{
		if ((i + rank) % 2 == 0)
		{
			if (rank < size - 1)
				MPI_Exchange(n / size, local_a, rank, rank + 1, comm);
		}
		else
		{
			if (rank > 0)
				MPI_Exchange(n / size, local_a, rank - 1, rank, comm);
		}
	
		MPI_Barrier(comm);
	}

	//gather local_a
	error = MPI_Gather(local_a, n / size, MPI_DOUBLE, array, n / size, MPI_DOUBLE, root, comm);
	if (error != MPI_SUCCESS)
	{
		free(local_a);
		return error;
	}
	free(local_a);

	return MPI_SUCCESS;
}

int MPI_Exchange(int n, double* array, int rank1, int rank2, MPI_Comm comm) {
	int rank, size, result, i, tag1 = 0, tag2 = 1;
	double* b = (double*)calloc(n, sizeof(double));
	double* c = NULL;

	MPI_Status status;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	if (rank == rank1)
	{
		result = MPI_Send(&array[0], n, MPI_DOUBLE, rank2, tag1, comm);
		result = MPI_Recv(&b[0], n, MPI_DOUBLE, rank2, tag2, comm, &status);
		c = merge(n, array, n, b);
		for (i = 0; i < n; i++)
		{
			array[i] = c[i];
		}
	}
	else if (rank == rank2)
	{
		result = MPI_Recv(&b[0], n, MPI_DOUBLE, rank1, tag1, comm, &status);
		result = MPI_Send(&array[0], n, MPI_DOUBLE, rank1, tag2, comm);
		c = merge(n, array, n, b);
		for (i = 0; i < n; i++)
		{
			array[i] = c[i + n];
		}
	}
	free(b);
	if (c)
		free(c);
	return MPI_SUCCESS;
}
double* merge(int n, double* a, int m, double* b) {
	int i, j, k;
	double* c = (double*)calloc(n + m, sizeof(double));

	for (i = j = k = 0; (i < n) && (j < m); )
	{
		if (a[i] <= b[j])
		{
			c[k++] = a[i++];
		}
		else
		{
			c[k++] = b[j++];
		}
	}
	if (i == n)
	{
		for (; j < m; )
		{
			c[k++] = b[j++];
		}
	}
	else
	{
		for (; i < n; )
		{
			c[k++] = a[i++];
		}
	}
	return c;
}
int MPI_Direct_Sort(int n, double* array, int root, MPI_Comm comm)
{
	int rank, size;
	MPI_Comm_size(comm, &size);
	MPI_Comm_size(comm, &rank);
	double* local_a = (double*)calloc(n / size, sizeof(double));

	int err = MPI_Scatter(array, n / size, MPI_DOUBLE, local_a, n / size, MPI_DOUBLE, root, comm);
	if (err != MPI_SUCCESS)
	{
		free(local_a);
		return err;
	}
	sort(n / size, local_a);

	int err1 = MPI_Gather(local_a, n / size, MPI_DOUBLE, array, n / size, MPI_DOUBLE, root, comm);
	if (err1 != MPI_SUCCESS)
	{
		free(local_a);
		return err;
	}

	if (rank == root)
	{
		for (int i = 0; i < size; i++)
		{
			double* aux = merge_array((i + 1) * n / size, array, n / size, array + (i + 1) * n / size);
			for (int j = 0; j < (i + 2) * n / size; j++)
			{
				array[j] = aux[j];
			}
		}
	}
	free(local_a);

	return MPI_SUCCESS;
}
void merge_sort(int n, double* a) {

	double* c;
	int i;

	if (n <= 1) return;
	if (n == 2) {

		if (a[0] > a[1]) {
			swap(&a[0], &a[1]);
		}
		return;
	}

	merge_sort(n / 2, a);
	merge_sort(n - n / 2, a + n / 2);

	c = merge_array(n / 2, a, n - n / 2, a + n / 2);

	for (i = 0; i < n; i++) {
		a[i] = c[i];
	}
}
void swap(double* a, double* b) {
	double temp;

	temp = *a;
	*a = *b;
	*b = temp;
}
double* merge_array(int n, double* a, int m, double* b) {

	int i, j, k;
	double* c = (double*)calloc(n + m, sizeof(double));

	for (i = j = k = 0; (i < n) && (j < m);) {
		if (a[i] <= b[j]) {
			c[k++] = a[i++];
		}
		else {
			c[k++] = b[j++];
		}
	}

	if (i == n) {
		for (; j < m; )
		{
			c[k++] = b[j++];
		}
	}
	else {
		for (; i < n; )
		{
			c[k++] = a[i++];
		}
	}
	return c;
}
void sort(int n, double* a) {

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (a[i] > a[j])
			{
				swap(&a[i], &a[j]);
			}
		}
	}
}
