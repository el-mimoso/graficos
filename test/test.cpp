#include <stdio.h>
#include <omp.h>

int main(void)
{
    int t = 100;
    int matrix[t][t];

#pragma omp parallel
    {
        int my_rank = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
        #pragma omp for schedule(dynamic, 3)
        for (int i = 0; i < t; i++)
        {
            for (int j = 0; j < t; j++)
            {
                matrix[i][j] = omp_get_thread_num();
            }
        }
    }

    FILE *f = fopen("MatrixTest.txt", "w");
    for (int i = 0; i < t; i++)
    {
        for (int j = 0; j < t; j++)
        {
            fprintf(f,"%d ", matrix[i][j]);
        }
        fprintf(f, "\n");
    }

    // #pragma omp parallel
    //     {
    //         int my_rank = omp_get_thread_num();
    //         int num_threads = omp_get_num_threads();
    //         printf("Un saludo desde el hilo %d de %d \n", my_rank, num_threads);
    //     }
    return 0;
}