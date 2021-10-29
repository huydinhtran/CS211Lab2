#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 2;
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    int maxind;
    int temps;
    double max;
    double tempv[n];
    int i, t, j, k;
    int a, b, c;
    for (i = 0 ; i < n-1 ; i++){
        maxind = i; 
        max = fabs(A[i*n+i]); 
        for (t = i+1 ; t < n ; t++){
            if (fabs(A[t*n+i]) > max ) {
                maxind = t; 
                max = fabs(A[t*n+i]); 
            }
        }
        if (max==0) 
            return -1; 
        else if (maxind != i ) {
            temps = ipiv[i]; 
            ipiv[i] = ipiv[maxind]; 
            ipiv[maxind] = temps;
            for (a=0 ; a < n ; a++) tempv[a] = A[i*n+a]; 
            for (b=0 ; b < n ; b++) A[i*n+b] = A[maxind*n+b]; 
            for (c=0 ; c < n ; c++) A[maxind*n+c] = tempv[c];
        }
        for (j = i+1 ; j < n ; j++) { 
            A[j*n+i] = A[j*n+i]/A[i*n+i];
            for (k = i+1 ; k < n ; k++) 
                A[j*n+k] = A[j*n+k] - A[j*n+i] * A[i*n+k]; 
        } 
    }
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
    //forward substitution for lower triangular
    if (UPLO == 'L'){
        int i;
        int a, b;
        double* y; 
        y = (double*) malloc (n * sizeof(double));
        y[0] = B[ipiv[0]];
        for (i=1 ; i<n ; i++){
            y[i] = B[ipiv[i]];
            for (a=0 ; a <= i-1 ; a++){
                y[i] -= A[i*n+a] * y[a];
            }
        }
        for (b=0 ; b<n ; b++) B[b] = y[b];
    }
    
    //backward substitution for upper triangular
    if (UPLO == 'U'){
        int i;
        int a, b;
        double* y; 
        y = (double*) malloc (n * sizeof(double));
        for (a=0 ; a<n ; a++) y[a] = B[a];
        y[n-1] = B[n-1] / A[(n-1)*n+n-1];        
        
        for (i=n-2 ; i>=0 ; i--){
            for (a=i+1 ; a < n ; a++){
                y[i] -= A[i*n+a] * y[a];
            }
            y[i] = y[i] / A[i*n + i];
        }
        for (b=0 ; b<n ; b++) B[b] = y[b];
    }
    
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
    /* Multiply n x n matrices a and b  */
    int i1, j1, k1;
    for (i = 0; i < n; i+=b)
        for (j = 0; j < n; j+=b)
            for (k = 0; k < n; k+=b)
             /* B x B mini matrix multiplications */
                for (i1 = i; i1 < i+b; i1++)
                    for (j1 = j; j1 < j+b; j1++){
                        register double C1=C[i1*n + j1];
                        register double C2=C[(i1+1)*n + j1];                           
                        register double C3=C[i1*n + (j1+1)];
                        register double C4=C[(i1+1)*n + (j1+1)];
                        
                        register double A1=A[i1*n + k1];                           
                        register double A2=A[i1*n + k1+1];
                        register double A3=A[(i1+1)*n + k1];
                        register double A4=A[(i1+1)*n + k1+1];
                        
                        register double B1=B[k1*n + j1];                           
                        register double B2=B[(k1+1)*n + j1];
                        register double B3=B[k1*n + (j1+1)];                           
                        register double B4=B[(k1+1)*n + (j1+1)];                        
                        
                        for (k1 = k; k1 < k+b; k1++){                            
                            C1 = A1 * B1 + A2 * B2 + C1;                                       
                            C2 = A3 * B1 + A4 * B3 + C2;                    
                            C3 = A1 * B3 + A2 * B4 + C3;                    
                            C4 = A3 * B3 + A4 * B4 + C4;
                        }                        
                        C[i1*n + j1]        =C1;
                        C[(i1+1)*n + j1]    =C2;
                        C[i1*n + (j1+1)]    =C3;
                        C[(i1+1)*n + (j1+1)]=C4;                        
                        
                        A[i1*n + k1]        =A1;                           
                        A[i1*n + k1+1]      =A2;
                        A[(i1+1)*n + k1]    =A3;
                        A[(i1+1)*n + k1+1]  =A4;
                        
                        B[k1*n + j1]        =B1;                           
                        B[(k1+1)*n + j1]    =B2;
                        B[k1*n + (j1+1)]    =B3;                           
                        B[(k1+1)*n + (j1+1)]=B4;
                    }
    return;
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    int ib, end;
    int maxind;       
    int temps;        
    double max;        
    double tempv[n];        
    int i, t, j, k;        
    int x, y, z, q, w, e;
    for ( ib = 1 ; ib < n-1 ; ib += b){
        end = ib + b-1;         
//         //apply BLAS2 version of GEPP to  get A(ib:n , ib:end) = P’ * L’ * U’
        for (i = ib ; i < end ; i++){
            maxind = i; 
            max = fabs(A[i*n+i]); 
            for (t = i+1 ; t < end+1 ; t++){
                if (fabs(A[t*n+i]) > max ) {
                    maxind = t; 
                    max = fabs(A[t*n+i]); 
                }
            }
            if (max==0) 
                return -1; 
            else if (maxind != i ) {
                temps = ipiv[i]; 
                ipiv[i] = ipiv[maxind]; 
                ipiv[maxind] = temps;
                for (x=0 ; x < n ; x++) tempv[x] = A[i*n+x]; 
                for (y=0 ; y < n ; y++) A[i*n+y] = A[maxind*n+y]; 
                for (z=0 ; z < n ; z++) A[maxind*n+z] = tempv[z];
            }
            for (j = i+1 ; j < n ; j++) { 
                A[j*n+i] = A[j*n+i]/A[i*n+i];
                for (k = i+1 ; k < end ; k++) 
                    A[j*n+k] = A[j*n+k] - A[j*n+i] * A[i*n+k]; 
            } 
        }
        
        double L[n][n]; 
        for(i=0;i<n;i++) 
            L[i][1]=A[i][1];
        
        for(j=2;j<=n;j++)
            U[1][j]=A[1][j]/L[1][1];
        for(i=0;i<n;i++)
            U[i][i]=1;
        for(i=2;i<=n;i++)
            for(j=2;j<=n;j++)
                if(i>=j){ 
                    L[i][j]=A[i][j];
                    for(k=1;k<=j-1;k++)
                        L[i][j]-=L[i][k]*U[k][j];
                }else{             
                    U[i][j]=A[i][j];
                    for(k=1;k<=j-1;k++)
                        U[i][j] = -L[i][k]*U[k][j];
                    U[i][j] /= L[i][i];
                }
       
        // A(ib:end , end+1:n) = LL-1 * A(ib:end , end+1:n)
        for (i = ib ; i < end ; i++){
            for (j = end+1 ; j < ib+b ; j++){
                A[i*n+j] = L * A[i*n+j];
            }
        }
        
//         // A(end+1:n , end+1:n ) -= A(end+1:n , ib:end) * A(ib:end , end+1:n)  
        for (i = end+1 ; i < ib+b ; i++){
            for (j = ib ; j < end ; j++){
                A[i*n+j] = A[i*n+j] * A[j*n+i];
            }
        }                        
    }  
    return 0;
}
