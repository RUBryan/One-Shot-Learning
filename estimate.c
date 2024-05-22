#include <stdio.h>
#include <stdlib.h>

double **alloc_matrix(int rows, int cols)
{
    double **m = malloc(rows * sizeof(double *));
    m[0] = malloc(rows * cols * sizeof(double));
    for (int i = 1; i < rows; i++)
        m[i] = &m[0][i * cols];
    return m;
}

void free_matrix(double **m, int rows, int cols)
{
    free(m[0]);
    free(m);
}

double **gaussJordan(int n, double **M)
{
    double **N = alloc_matrix(n,n);
    for (int i = 0; i < n; i++)
    {
        N[i] = &N[0][i * n];
        for (int j = 0; j < n; j++)
        {
            N[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    double factor;
    int i, j, k;

    for (k = 0; k < n; k++)
    {
        factor = M[k][k];

    
        for (j = 0; j < n; j++)
        {
            M[k][j] /= factor;
            N[k][j] /= factor;
        }

        for (i = 0; i < n; i++)
        {
            if (i != k)
            {
                factor = M[i][k];
                for (j = 0; j < n; j++)
                {
                    M[i][j] -= factor * M[k][j];
                    N[i][j] -= factor * N[k][j];
                }
            }
        }
    }
    for (k = n - 1; k >= 0; k--)
    {
        for (i = k - 1; i >= 0; i--)
        {
            factor = M[i][k];
            for (j = 0; j < n; j++)
            {
                M[i][j] -= factor * M[k][j];
                N[i][j] -= factor * N[k][j];
            }
        }
    }

    return N;
}

void matrixProduct(int rows1, int cols1, double **Amatrix, int cols2, double **Bmatrix, double **Cmatrix)
{
    for (int i = 0; i < rows1; i++)
    {
        for (int j = 0; j < cols2; j++)
        {
            Cmatrix[i][j] = 0;
            for (int k = 0; k < cols1; k++)
            {
                Cmatrix[i][j] += Amatrix[i][k] * Bmatrix[k][j];
            }
        }
    }
}

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        fprintf(stderr, "Usage: %s <input_file> <input_file2>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    FILE *f = fopen(argv[1], "r");
    if (f == NULL)
    {
        fprintf(stderr, "Unable to open %s\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    char train[5];
    int rows, cols, result;
    result = fscanf(f, "%5c %d %d", train, &cols, &rows);

    if (result != 3)
    {
        fprintf(stderr, "Failed to read matrix dimensions.\n");
        fclose(f);
        exit(EXIT_FAILURE);
    }

    cols = 1 + cols;

    // Create space for new data (malloc)
    double **A = alloc_matrix(rows, cols);

    double **lastNumbers = alloc_matrix(rows, 1);
    // Initialize the first column of the matrix with 1
    for (int i = 0; i < rows; i++)
    {
        A[i][0] = 1.0;
    }

    int value = 0;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 1; j <= cols; j++)
        {
            if (j == cols)
            {
                // Read and discard the last number in the row
                fscanf(f, "%lf", &lastNumbers[value][0]);
                value++;
            }
            else
            {
                // Read and store the other numbers in the matrix
                fscanf(f, "%lf", &A[i][j]);
            }
        }
    }


    int transpose_rows1 = cols;
    int transpose_colums1 = rows;

    // Create space for new data (malloc)
    double **transpose_matrix1 = (double **)malloc(transpose_rows1 * sizeof(double *));
    transpose_matrix1[0] = (double *)malloc(transpose_rows1 * transpose_colums1 * sizeof(double));
    for (int i = 0; i < transpose_rows1; i++)
    {
        transpose_matrix1[i] = &transpose_matrix1[0][i * transpose_colums1];
    }

    for (int i = 0; i < cols; i++)
        for (int j = 0; j < rows; j++)
            transpose_matrix1[i][j] = A[j][i];

    double **XTX = alloc_matrix(transpose_rows1, transpose_rows1);
    matrixProduct(transpose_rows1, transpose_colums1, transpose_matrix1, cols, A, XTX);
    free_matrix(A, rows, cols);

    double **inverseMatrix = gaussJordan(transpose_rows1, XTX);
    free_matrix(XTX, transpose_rows1, transpose_rows1);

    double **XTY = alloc_matrix(rows, 1);

    matrixProduct(transpose_rows1, transpose_colums1, transpose_matrix1, 1, lastNumbers, XTY);
    free_matrix(lastNumbers, rows, 1);
    free_matrix(transpose_matrix1, transpose_rows1, transpose_colums1);

    double **W = alloc_matrix(rows, 1);

    matrixProduct(transpose_rows1, transpose_rows1, inverseMatrix, 1, XTY, W);
    free_matrix(inverseMatrix, transpose_rows1, transpose_colums1);
    free_matrix(XTY, transpose_rows1, 1);

    fclose(f);

    FILE *z = fopen(argv[2], "r");
    if (z == NULL)
    {
        fprintf(stderr, "Unable to open %s\n", argv[2]);
        exit(EXIT_FAILURE);
    }

    char data[4];
    int datarows, datacols, result2;
    result2 = fscanf(z, "%4c %d %d", data, &datacols, &datarows);

    if (result2 != 3)
    {
        fprintf(stderr, "Failed to read matrix dimensions.\n");
        fclose(f);
        exit(EXIT_FAILURE);
    }

    datacols = datacols + 1;

    double **dataMatrix = alloc_matrix(datarows, datacols);

    // printf("break\n");

    for (int i = 0; i < datarows; i++)
    {
        dataMatrix[i][0] = 1.0;
    }
    for (int i = 0; i < datarows; i++)
    {
        for (int j = 1; j < datacols; j++)
        {
            // Read and store the other numbers in the matrix
            fscanf(z, "%lf", &dataMatrix[i][j]);
        }
    }

    double **Y = alloc_matrix(rows, 1);
    matrixProduct(datarows, datacols, dataMatrix, 1, W, Y);
    free_matrix(dataMatrix, datarows, datacols);
    free_matrix(W, datarows, 1);

    for (int i = 0; i < datarows; i++)
    {

        printf("%.0f", Y[i][0]);

        printf("\n");
    }

    free_matrix(Y, datarows, 1);
    fclose(z);

    return EXIT_SUCCESS;
}