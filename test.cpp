#include <stdio.h>
#include <omp.h>
int main()
{
    int i = 256; // a shared variable
#pragma omp parallel
    {
        int x; // a variable local or private to each thread
        x = omp_get_thread_num();
        printf("x = %d, i = %d\n", x, i);
    }
}