#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>
typedef int pixel;
// #define max(a,b) (((a) > (b)) ? (a) : (b))
// #define min(a,b) (((a) < (b)) ? (a) : (b))
#define fastmin(a, b) (a < b ? a : b)
#define fastmax(a, b) (a > b ? a : b)
#define MAT_SIZE 4096
#define COPY(d, s) *(d) = *(s)
#define assign_sum_to_pixel_fast(a, b) *(a) = b;
int RIDX(int x, int y, int n)
{
    return x * n + y;
}
int max(int a, int b)
{
    return a > b ? a : b;
}
int min(int a, int b)
{
    return a > b ? b : a;
}
void initialize_pixel_sum(int *sum)
{
    *sum = 0;
}
void accumulate_sum(int *sum, int src)
{
    *sum += src;
}
void assign_sum_to_pixel(pixel *obj, int sum)
{
    *obj = sum;
}
static pixel avg(int dim, int i, int j, pixel *src)
{
    int sum = 0;
    initialize_pixel_sum(&sum);
    for (int jj = max(j - 1, 0); jj <= min(j + 1, dim - 1); jj++)
        for (int ii = max(i - 1, 0); ii <= min(i + 1, dim - 1); ii++)
            accumulate_sum(&sum, src[RIDX(ii, jj, dim)]);
    assign_sum_to_pixel(&src[RIDX(i, j, dim)], sum);
    return src[RIDX(i, j, dim)];
}
void show(pixel *src, int dim)
{
    printf("picture matrix:\n");
    for (int i = 0; i < dim; i++)
    {

        for (int j = 0; j < dim; j++)
        {
            printf("%d ", src[RIDX(i, j, dim)]);
        }
        printf("\n");
    }
    printf("end---------------\n");
}
void fill(pixel *src, int dim)
{
    srand(2);
    for (int i = 0; i < dim * dim; i++)
    {
        src[i] = rand() % 10;
    }
}
void naive_rotate(int dim, pixel *src, pixel *dst)
{
    int i, j;

    for (int i = 0; i < dim; i++)
    {

        for (int j = 0; j < dim; j++)
        {
            dst[RIDX(dim - 1 - j, i, dim)] = src[RIDX(i, j, dim)];
        }
    }
}

void rotate_1_1(int dim, pixel *src, pixel *dst)
{
    int i, j, ii, jj;
    for (ii = 0; ii < dim; ii += 4)
        for (jj = 0; jj < dim; jj += 4)
            for (i = ii; i < ii + 4; i++)
                for (j = jj; j < jj + 4; j++)
                    dst[RIDX(dim - 1 - j, i, dim)] = src[RIDX(i, j, dim)];
}
void rotate_1_2(int dim, pixel *src, pixel *dst)
{
    int i, j, ii, jj;
    for (ii = 0; ii < dim; ii += 32)
        for (jj = 0; jj < dim; jj += 32)
            for (i = ii; i < ii + 32; i += 4)
                for (j = jj; j < jj + 32; j += 4)
                {
                    dst[RIDX(dim - 1 - j, i, dim)] = src[RIDX(i, j, dim)];
                    dst[RIDX(dim - 1 - j, i + 1, dim)] = src[RIDX(i + 1, j, dim)];
                    dst[RIDX(dim - 1 - j, i + 2, dim)] = src[RIDX(i + 2, j, dim)];
                    dst[RIDX(dim - 1 - j, i + 3, dim)] = src[RIDX(i + 3, j, dim)];
                    dst[RIDX(dim - 1 - j - 1, i, dim)] = src[RIDX(i, j + 1, dim)];
                    dst[RIDX(dim - 1 - j - 1, i + 1, dim)] = src[RIDX(i + 1, j + 1, dim)];
                    dst[RIDX(dim - 1 - j - 1, i + 2, dim)] = src[RIDX(i + 2, j + 1, dim)];
                    dst[RIDX(dim - 1 - j - 1, i + 3, dim)] = src[RIDX(i + 3, j + 1, dim)];
                    dst[RIDX(dim - 1 - j - 2, i, dim)] = src[RIDX(i, j + 2, dim)];
                    dst[RIDX(dim - 1 - j - 2, i + 1, dim)] = src[RIDX(i + 1, j + 2, dim)];
                    dst[RIDX(dim - 1 - j - 2, i + 2, dim)] = src[RIDX(i + 2, j + 2, dim)];
                    dst[RIDX(dim - 1 - j - 2, i + 3, dim)] = src[RIDX(i + 3, j + 2, dim)];
                    dst[RIDX(dim - 1 - j - 3, i, dim)] = src[RIDX(i, j + 3, dim)];
                    dst[RIDX(dim - 1 - j - 3, i + 1, dim)] = src[RIDX(i + 1, j + 3, dim)];
                    dst[RIDX(dim - 1 - j - 3, i + 2, dim)] = src[RIDX(i + 2, j + 3, dim)];
                    dst[RIDX(dim - 1 - j - 3, i + 3, dim)] = src[RIDX(i + 3, j + 3, dim)];
                }
}
void rotate_1_3(int dim, pixel *src, pixel *dst)
{
    int i, j, ii, jj;
    for (ii = 0; ii < dim; ii += 4)
        for (jj = 4; jj < dim + 4; jj = (jj + 4) % dim)
        {
            for (i = ii; i < ii + 4; i++)
                for (j = jj; j < jj + 4; j++)
                    dst[RIDX(dim - 1 - j, i, dim)] = src[RIDX(i, j, dim)];
            if (jj == 0)
                break;
        }
}
void rotate_1_4(int dim, pixel *src, pixel *dst)
{
    int i, j;
    for (i = 0; i < dim; i += 32)
        for (j = dim - 1; j >= 0; j -= 1)
        {
            pixel *dptr = dst + RIDX(dim - 1 - j, i, dim);
            pixel *sptr = src + RIDX(i, j, dim);
            COPY(dptr, sptr);
            sptr += dim;
            COPY(dptr + 1, sptr);
            sptr += dim;
            COPY(dptr + 2, sptr);
            sptr += dim;
            COPY(dptr + 3, sptr);
            sptr += dim;
            COPY(dptr + 4, sptr);
            sptr += dim;
            COPY(dptr + 5, sptr);
            sptr += dim;
            COPY(dptr + 6, sptr);
            sptr += dim;
            COPY(dptr + 7, sptr);
            sptr += dim;
            COPY(dptr + 8, sptr);
            sptr += dim;
            COPY(dptr + 9, sptr);
            sptr += dim;
            COPY(dptr + 10, sptr);
            sptr += dim;
            COPY(dptr + 11, sptr);
            sptr += dim;
            COPY(dptr + 12, sptr);
            sptr += dim;
            COPY(dptr + 13, sptr);
            sptr += dim;
            COPY(dptr + 14, sptr);
            sptr += dim;
            COPY(dptr + 15, sptr);
            sptr += dim;
            COPY(dptr + 16, sptr);
            sptr += dim;
            COPY(dptr + 17, sptr);
            sptr += dim;
            COPY(dptr + 18, sptr);
            sptr += dim;
            COPY(dptr + 19, sptr);
            sptr += dim;
            COPY(dptr + 20, sptr);
            sptr += dim;
            COPY(dptr + 21, sptr);
            sptr += dim;
            COPY(dptr + 22, sptr);
            sptr += dim;
            COPY(dptr + 23, sptr);
            sptr += dim;
            COPY(dptr + 24, sptr);
            sptr += dim;
            COPY(dptr + 25, sptr);
            sptr += dim;
            COPY(dptr + 26, sptr);
            sptr += dim;
            COPY(dptr + 27, sptr);
            sptr += dim;
            COPY(dptr + 28, sptr);
            sptr += dim;
            COPY(dptr + 29, sptr);
            sptr += dim;
            COPY(dptr + 30, sptr);
            sptr += dim;
            COPY(dptr + 31, sptr);
        }
}
int main()
{
    int *a = (int *)malloc(MAT_SIZE * MAT_SIZE * sizeof(int));
    int *b = (int *)malloc(MAT_SIZE * MAT_SIZE * sizeof(int));
    fill(a, MAT_SIZE);
    //show(a, MAT_SIZE);
    long begin = clock();
    // naive_rotate(MAT_SIZE, a, b);
    // printf("无改进耗时: %lf\n", (double)(clock() - begin) / CLOCKS_PER_SEC);
    // begin = clock();
    // rotate_1_1(MAT_SIZE, a, b);
    // printf("分块耗时: %lf\n", (double)(clock() - begin) / CLOCKS_PER_SEC);
    // begin = clock();
    // rotate_1_2(MAT_SIZE, a, b);
    // printf("循环展开耗时: %lf\n", (double)(clock() - begin) / CLOCKS_PER_SEC);
    // begin = clock();
    // rotate_1_3(MAT_SIZE, a, b);
    // printf("不同巡回路线耗时: %lf\n", (double)(clock() - begin) / CLOCKS_PER_SEC);
    begin = clock();
    rotate_1_4(MAT_SIZE, a, b);
    printf("最后尝试耗时: %lf\n", (double)(clock() - begin) / CLOCKS_PER_SEC);
    //show(b, MAT_SIZE);
    return 0;
}