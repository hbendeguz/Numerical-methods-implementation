#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void kiir(double * mat, int N){
    for(int i=0; i<N; i++){
        for(int j=0;j<N; j++){
            printf("%lf\t", mat[i*N+j]);
        }
        printf("\n");
    }

}
/*double determinant(double *mat, int N){
    double det=1;
    for(int c=0; c<N; c++){
        double detplus=1;
        double detminus= 1;
        for(int i=0; i<N; i++){
            double detplus = detplus * mat[i*N+i+c];
            double detminus = detminus * mat[i*N+(N-i)-c];
        }
         det = det+detplus-detminus;
    }
    return det ;
}*/
void sorszor(double *mat, int N, int j,double q){
    for (int i=0; i<N; i++){
        mat[j*N+i] = q * mat[j*N+i];
    }
}


void eleminate(double *mat, int N, int i, int j,double *eye,int inv,double *vec){
    double q = mat[j*N+i]/ mat[N*i+i];
    for(int c=0;c<N; c++ ){
        eye[j*N+c]=eye[j*N+c]-q*eye[N*i+c];
        mat[j*N+c]=mat[j*N+c]-q*mat[N*i+c];
    }
       vec[j]=vec[j]-q*vec[i];
}
/*
int oszlopmax(double *mat, int N, int j){
    double max = mat[j];
    int hely = 1;
   // printf("%lf",max);
    for(int c=0; c < N;c++){
        if (fabs(mat[j+N*c]) > max){
            max=fabs(mat[j+N*c]);
            hely = c+1;}
    }
    return hely;
}
*/
void rowSwitch(double *mat, int l, int k, int N, double *eye,int inv,double *vec){
    for(int j=0; j<N; j++){
       double temp=mat[l*N+j];
        mat[l*N+j]=mat[(k)*N+j];
        mat[(k)*N+j]=temp;
        
        double demp=eye[l*N+j];
        eye[l*N+j]=eye[(k)*N+j];
        eye[(k)*N+j]=demp;}
         double gemp=vec[l];
        vec[l]=vec[k];
        vec[k]=gemp;
}

double sormax(double * mat, int N, int j){
    int c=0;
    double maximum;
    maximum = fabs(mat[j*N]);
    for (c = 0; c < N; c++)
    {
        if (fabs(mat[c+j*N]) > maximum)
        {
            maximum  = fabs(mat[c+j*N]);
                }
    }
    return maximum;
}
void Gauss(double * mat, int N, double * eye, int inv, int j, double *vec){
    //for (int j=0; j<N; j++) {
        double temp = mat[j*N+j];
    
        for (int i=0; i<N;i++){
        
         eye[j*N+i]=eye[j*N+i]/temp;
         mat[j*N+i]=mat[j*N+i]/temp;
        }
             vec[j]= vec[j]/temp;
}




int main(int argc, char *argv[])
{
    double *eye;
    double *matrix;
    double *vec;
    int N;
    int i,j;
    int inv;
    N = atoi(argv[1]);
    if (argc <4){
        inv =1;
    }
    if (argc >=4){
        inv =0;
        vec = (double *)calloc(N, sizeof(double));
        FILE *g;
        g=fopen(argv[3], "r");
                    for(i = 0; i < N; i++)
            {
                fscanf(g, "%lf" , &vec[i]);
            }
        
        fclose(g);
    }   else{vec = (double *)calloc(N, sizeof(double));}
    
    eye = (double *)calloc(N * N, sizeof(double));
    for(int i=0;i<N; i++){
        eye[i*N+i]=1;
    }
    matrix = (double*)calloc(N * N, sizeof(double));
    
    FILE *f;
    f=fopen(argv[2], "r");
    for (j = 0; j < N; j++)
    {
        for(i = 0; i < N; i++)
        {
            fscanf(f, "%lf" , &matrix[N * j + i]);
        }
    }
    fclose(f);
    printf("A bemeneti mátrix:\n");
    kiir(matrix,N);
    printf("\n");
  
    //  for(int j=0; j<N; j++){
    //    int f = oszlopmax(matrix,N,j);
      //  printf("%d\n", f);
   //}
  
 
    for(int i=0; i<N; i++){
        double f = sormax(matrix, N, i);
        for(int j=0;j<N;j++){
            matrix[i*N+j]/=f;
            eye[i*N+j]/=f;
        }
        vec[i]/=f;
    }
   //alsó eltüntetése
    
    for(int i=0; i<N-1;i++){
        if (matrix[i*N+i]!=0){Gauss(matrix,N,eye,inv,i,vec);
        for(int j=i+1;j<N; j++){
            eleminate(matrix, N, i, j,eye,inv,vec);
        }
        }   else{ int k=i+1;
            while (matrix[i*N+i]==0 || k==N){
            rowSwitch(matrix,i, k, N,eye,inv,vec);
            k++;}
                       }
            
            Gauss(matrix,N,eye,inv,i,vec);
            for(int j=i+1;j<N; j++){
                eleminate(matrix, N, i, j,eye,inv,vec);
        }
    }


    
   // felső eltüntetése
    for(int i=N-1; i>0;i--){
        for(int j=i-1; j > -1; j--){
            eleminate(matrix, N, i, j,eye,inv,vec);
        }
        Gauss(matrix,N,eye,inv,i,vec);
    }

    printf("A mátrix inverze:\n");
    kiir(eye,N);
    printf("\n");
    printf("Az egységmátrixxá eleminált mátrix:\n");

    kiir(matrix,N);
    printf("\n");
    if(argc == 4){
        printf("Az eredményvektor: \n");
        for(int k=0; k<N; k++){
        
        printf("%lf\n", vec[k]);
    }
}
    
    free(matrix);
    free(eye);
    free(vec);
    return 0;
   }
