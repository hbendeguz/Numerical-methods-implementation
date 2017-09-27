#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    double *matrix;
    double *vec;
    double *mult;
    int rows,cols;
    int i,j;
    
   /* if (argc != 3)
    {
        fprintf(stderr, "Usage: rowToCol <rows> <cols>  < infile > outfile\n");
        exit(-1);
    }
    */
    rows = atoi(argv[1]);
    cols = atoi(argv[2]);
    
    matrix = (double*)calloc(rows * cols, sizeof(double));
    FILE *f;
    f=fopen(argv[3], "r");
    for (j = 0; j < rows; j++)
    {
        for(i = 0; i < cols; i++)
        {
            fscanf(f, "%lf" ,&matrix[cols * j + i]);
        }
    }
    fclose(f);
    
   /* for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            printf("%f ", matrix[cols * i + j]);
        }
        printf("\n");
    }*/
    vec=(double*)calloc(cols, sizeof(double));
    mult=(double*)calloc(rows, sizeof(double));
    
    FILE *g;
    g=fopen(argv[4], "r");
    for (j=0; j < cols; j++){
        fscanf(g, "%lf", &vec[j]);
    }
   
    for(j=0;j < cols; j++){
    for(i = 0; i < rows; i++){
        mult[i] = mult[i]+vec[j] * matrix[cols * i + j];
    }
    }
    fclose(g);
    free(matrix);
    free(vec);
    for (i = 0; i < rows; i++){
        printf( "%f\n",mult[i]); /* Sorvektor az eredmény, de jobban áttekinthetőnek tűnik, ha oszlopvektorként van kiíratva. Ha sorvektort szeretnénk csak a \n kell kicserélni \t-re*/
    }
    free(mult);
    
    return 0;
    
}
