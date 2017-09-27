//
//  main.c
//  poly
//
//  Created by bendeguz on 2017. 04. 07..
//  Copyright © 2017. bendeguz. All rights reserved.
//


void makeXTranspose(double * X, double * XT, int rowsOfData, int M, int independebtVariables,int order ){
                          for(int i=0; i< M; i++){
                              for(int j=0; j< rowsOfData; j++){
                                  XT[i*rowsOfData+j]=X[j*(independebtVariables * order+1)+i];
                              }
                          }
}
void makeXTX(double * XTX,double *XT, double *X, int M, int K){
    for(int i=0; i<M; i++){
        for(int j=0; j<M; j++){
            double sum = 0;
            for(int l=0; l < K; l++){
                
                sum = sum+XT[i*K+l]*X[l*M+j];
            }
            XTX[i*(M)+j]=sum;
        }
    }
}
void divideBySigma(double *sigma, double *resultValues, double *X,int M,int K){
    for(int i=0; i<M; i++){
        for (int j=0; j<K;j++){
            X[i*K+j]/=sigma[i];
            
        }
        resultValues[i]/=sigma[i];
    }
}
void multiplyBySigma(double *sigma, double *resultValues, double *X,int M,int K){
    for(int i=0; i<M; i++){
        for (int j=0; j<K;j++){
            X[i*K+j]*=sigma[i];
            
        }
        resultValues[i]*=sigma[i];
    }
}

void makeXTy(double *mult, double *vec,double *matrix, int cols, int rows){
    for(int j=0;j < cols; j++){
        for(int i = 0; i < rows; i++){
            mult[i] = mult[i]+vec[j] * matrix[cols * i + j];
        }
    }

}

void makeToPlot(double *mult, double *vec,double *matrix, int cols, int rows){
    for(int j=0;j <cols; j++){
        for(int i = 0; i < rows; i++){
            mult[i] = mult[i]+vec[j] * matrix[cols * i + j];
        }
    }
    
}


void Norm(double * mat, int N, int j, double *vec){
        double temp = mat[j*N+j];
    for (int i=0; i<N;i++){
        mat[j*N+i]=mat[j*N+i]/temp;
    }
    vec[j]= vec[j]/temp;
}
void eleminate(double *mat, int N, int i, int j,double *vec){
    double q = mat[j*N+i]/ mat[N*i+i];
    for(int c=0;c<N; c++ ){
        mat[j*N+c]=mat[j*N+c]-q*mat[N*i+c];
    }
    vec[j]=vec[j]-q*vec[i];
}
void rowSwitch(double *mat, int l, int k, int N, double *vec){
    for(int j=0; j<N; j++){
        double temp=mat[l*N+j];
        mat[l*N+j]=mat[(k)*N+j];
        mat[(k)*N+j]=temp;
        
        double gemp=vec[l];
        vec[l]=vec[k];
        vec[k]=gemp;
}

}

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


int main(int argc, const char * argv[]) {
    int order = atoi(argv[2]);
    int independentVariables = atoi(argv[3]);
    int rowsOfData = atoi(argv[4]);
    FILE * f;
    f = fopen(argv[1], "r");
    double * data;
    data = (double *)calloc((independentVariables+2)*rowsOfData, sizeof(double));
    for(int i=0; i<independentVariables+2; i++){
        for(int j=0; j<rowsOfData; j++){
            fscanf(f,"%lf", &data[i*rowsOfData+j]);
        }
    }
    fclose(f);
    double * sigma;
    sigma = (double *)calloc(rowsOfData, sizeof(double));
    for(int i=0; i < rowsOfData; i++){
        sigma[i] = data[independentVariables+1+i*(independentVariables+2)];
    }
    double * resultValues;
    resultValues = (double*)calloc(rowsOfData, sizeof(double));
    for(int i=0; i < rowsOfData; i++){
        resultValues[i] = data[independentVariables+i*(independentVariables+2)];
    }
    double *toplot;
    toplot = (double *)calloc(rowsOfData, sizeof(double));
    
    double * XT;
    XT = (double *)calloc(rowsOfData*(independentVariables * order+1), sizeof(double));
    
    double * XTX;
    XTX = (double *)calloc((independentVariables* order +1)*(independentVariables * order+1), sizeof(double));
    
    double *X;
    X=(double *)calloc(rowsOfData * (independentVariables * order+1), sizeof(double));
    
    double *XTy;
    XTy=(double *)calloc(independentVariables*order+1, sizeof(double));
    
    
    for(int j=0; j< rowsOfData; j++){
        X[j*(independentVariables*order +1)] = 1;
        for(int i=0; i < independentVariables; i++){
            for(int k=1; k < (order+1); k++ ){
            X[j*(independentVariables * order +1)+k+i*(order)] = pow(data[j*(independentVariables+2)+i],k);
            }
        }
    }
    
    divideBySigma(sigma, resultValues, X, rowsOfData, independentVariables*order+1);

  /*
   for(int i=0; i<rowsOfData; i++ ){
        for(int j=0; j<(independentVariables * order +1); j++ ){
            printf("%lf ", X[i*(independentVariables * order +1)+j]);
        }
        printf("\n");
    }
    printf("\n");*/
    
    makeXTranspose(X, XT, rowsOfData,independentVariables * order+1, independentVariables, order );
  /* for(int i=0; i<(independentVariables * order +1); i++ ){
       for(int j=0; j<rowsOfData; j++ ){
           printf("%lf ", XT[i*(rowsOfData)+j]);
      }
        printf("\n");
    }
   printf("\n");
    for(int i=0; i<rowsOfData; i++){
        printf("%lf\n", resultValues[i]);
    }*/
    
    makeXTX(XTX, XT, X, (independentVariables*order+1), rowsOfData);
    
    makeXTy(XTy, resultValues, XT,rowsOfData, (independentVariables*order+1));
    

  /*  for(int i=0; i<independentVariables*order+1; i++){
        for(int j=0; j<independentVariables*order+1; j++){
            printf("%lf  ", XTX[i*(independentVariables*order+1)+j]);
        }printf("\n");
    }
    for(int i=0; i<independentVariables*order+1; i++){
        printf("%lf\n", XTy[i]);
    }*/
    
    
    //alsó eltüntetése
   //itt
    int N=(independentVariables*order +1);
    for(int i=0; i<(independentVariables*order +1)-1;i++){
        if (XTX[i*(independentVariables*order +1)+i]!=0){
            Norm(XTX,(independentVariables*order +1),i,XTy);
            for(int j=i+1;j<N; j++){
                eleminate(XTX, N, i, j, XTy);
            }
        }   else{ int k=i+1;
                while (XTX[i*(independentVariables*order +1)+i]==0 || k==N){
                rowSwitch(XTX,i, k, (independentVariables*order +1),XTy);
                k++;
                }
            }
        
        Norm(XTX,(independentVariables*order +1),i,XTy);
        for(int j=i+1;j<(independentVariables*order +1); j++){
            eleminate(XTX, (independentVariables*order +1), i, j,XTy);
        }
    }

    
    // felső eltüntetése
    for(int i=(independentVariables*order +1)-1; i>0;i--){
        for(int j=i-1; j > -1; j--){
            eleminate(XTX, (independentVariables*order +1), i, j,XTy);
                        }
        Norm(XTX,(independentVariables*order +1),i,XTy);
    }

    for(int i=0; i<(independentVariables*order +1); i++){
        printf("%lf\n", XTy[i]);
    }
      //idáig
    multiplyBySigma(sigma, resultValues, X, rowsOfData, independentVariables*order+1);
   // makeXTy(toplot, XTy, X,(order*independentVariables+1) , rowsOfData);
    FILE * g;
    g = fopen("resultFile.txt", "w");
    makeToPlot(toplot, XTy, X,(order*independentVariables+1), rowsOfData );
    for(int i= 0; i< rowsOfData; i++){
        fprintf(g,"%lf %lf\n",resultValues[i], toplot[i]);
    }
    
    fclose(g);
    
    free(X);
    free(data);
    free(XT);
    free(XTX);
    free(sigma);
    free(resultValues);
    free(XTy);
    free(toplot);

    return 0;

}
