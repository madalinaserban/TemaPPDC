#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
void bitonicSort(int a[], int low, int cnt, int dir);
void bitonicMerge(int a[], int low, int cnt, int dir);
void exchange(int a[], int i, int j);
void compare(int a[], int i, int j, int dir);
int isSorted(int arr[], int n);
void MPI_Bitonic_sort(int a[], int n, int dir);
void MPI_Shell_sort(int arr[], int n);
void MPI_Bucket_sort(int arr[], int n);


int main(int argc, char* argv[])
{
	int rank, size, i, n, * arr, * recvbuf, * bitonicRecvbuf;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	n = 1000000;
	arr = (int*)malloc(n * sizeof(int));

	if (rank == 0)
	{
		for (i = 0; i < n; i++)
			arr[i] = rand();
	}

	double start_time, end_time;
	MPI_Barrier(MPI_COMM_WORLD);
	start_time = MPI_Wtime();

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int chunk_size = n / size;
	int remainder = n % size;
	int* sendcounts = (int*)malloc(size * sizeof(int));
	int* displs = (int*)malloc(size * sizeof(int));

	for (i = 0; i < size; i++)
	{
		sendcounts[i] = chunk_size;
		if (remainder > 0)
		{
			sendcounts[i]++;
			remainder--;
		}
		displs[i] = (i > 0) ? (displs[i - 1] + sendcounts[i - 1]) : 0;
	}

	recvbuf = (int*)malloc(sendcounts[rank] * sizeof(int));
	bitonicRecvbuf = (int*)malloc(sendcounts[rank] * sizeof(int));

	MPI_Scatterv(arr, sendcounts, displs, MPI_INT, recvbuf, sendcounts[rank], MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Shell_sort(recvbuf, sendcounts[rank]);

	MPI_Gatherv(recvbuf, sendcounts[rank], MPI_INT, arr, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		end_time = MPI_Wtime();

		if (isSorted(arr, n))
		{
			printf("Shell sort: Array is sorted\n");
		}
		else
		{
			printf("Shell sort: Array is NOT sorted\n");
		}

		printf("Shell sort: Time taken = %lf seconds\n", end_time - start_time);

		start_time = MPI_Wtime();
	}

	MPI_Bucket_sort(arr, n);

	if (rank == 0)
	{
		end_time = MPI_Wtime();

		if (isSorted(arr, n))
		{
			printf("Bucket sort: Array is sorted\n");
		}
		else
		{
			printf("Bucket sort: Array is NOT sorted\n");
		}

		printf("Bucket sort: Time taken = %lf seconds\n", end_time - start_time);

		start_time = MPI_Wtime();
	}


	MPI_Bitonic_sort(recvbuf, sendcounts[rank], 1);

	MPI_Gatherv(recvbuf, sendcounts[rank], MPI_INT, arr, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		end_time = MPI_Wtime();

		if (isSorted(arr, n))
		{
			printf("Bitonic sort: Array is sorted\n");
		}
		else
		{
			printf("Bitonic sort: Array is NOT sorted\n");
		}

		printf("Bitonic sort: Time taken = %lf seconds\n", end_time - start_time);

		free(arr);
	}

	free(recvbuf);
	free(sendcounts);
	free(displs);
	free(bitonicRecvbuf);

	MPI_Finalize();
	return 0;
}
void MPI_Shell_sort(int arr[], int n)
{
	int i, j, gap, temp;
	for (gap = n / 2; gap > 0; gap /= 2)
	{
		for (i = gap; i < n; i++)
		{
			for (j = i - gap; j >= 0 && arr[j] > arr[j + gap]; j -= gap)
			{
				temp = arr[j];
				arr[j] = arr[j + gap];
				arr[j + gap] = temp;
			}
		}
	}
}

void exchange(int a[], int i, int j)
{
	int t = a[i];
	a[i] = a[j];
	a[j] = t;
}
void compare(int a[], int i, int j, int dir)
{
	if (dir == (a[i] > a[j]))
		exchange(a, i, j);
}
void bitonicMerge(int a[], int low, int cnt, int dir)
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
void bitonicSort(int a[], int low, int cnt, int dir)
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


void MPI_Bitonic_sort(int a[], int n, int dir)
{
	bitonicSort(a, 0, n, dir);
}

int isSorted(int arr[], int n)
{
	for (int i = 1; i < n; i++)
	{
		if (arr[i] < arr[i - 1])
		{
			return 0;
		}
	}
	return 1;
}
void MPI_Bucket_sort(int arr[], int n)
{
	int i, j, count;
	int* bucket = (int*)malloc(n * sizeof(int));

	for (i = 0; i < n; i++)
	{
		bucket[i] = 0;
	}

	for (i = 0; i < n; i++)
	{
		(bucket[arr[i]])++;
	}

	for (i = 0, j = 0; i < n; i++)
	{
		for (count = bucket[i]; count > 0; count--)
		{
			arr[j++] = i;
		}
	}

	free(bucket);
}
