#include "../include/for_you_to_do.h"



int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 1;
  
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
//     int maxind;
//     int temps;
//     double max;
//     double tempv[n];
//     int i, t, j, k;
//     int a, b, c;
//     for (i = 0 ; i < n-1 ; i++){
//         maxind = i; 
//         max = abs(A[i*n+i]); 
//         for (t = i ; t < n ; t++){
//             if (abs(A[t*n+i]) > max ) {
//                 maxind = t; 
//                 max = abs(A[t*n+i]); 
//             }
//         }
//         if (max==0) 
//             return -1; 
//         else if (maxind != i ) {
//             temps = ipiv[i]; 
//             ipiv[i] = ipiv[maxind]; 
//             ipiv[maxind] = temps;
//             for (a=0 ; a < n ; a++) tempv[a] = A[i*n+a]; 
//             for (b=0 ; b < n ; b++) A[i*n+b] = A[maxind*n+b]; 
//             for (c=0 ; c < n ; c++) A[maxind*n+c] = tempv[c];
//         }
//         for (j = i ; j < n ; j++) { 
//             A[j*n+i] = A[j*n+i]/A[i*n+i];
//             for (k = i ; k < n ; k++) 
//                 A[j*n+k] = A[j*n+k] - A[j*n+i] * A[i*n+k]; 
//         } 
//     }
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
        int a;
        double sum;
        double* y; 
        y = (double*) malloc (n * sizeof(double));
        y[0] = B[ipiv[0]];
        for (i=1 ; i<n ; i++){
            for (a=0 ; a < i-2 ; a++){
                sum += y[a] * A[i*n+a];
            }
            y[i] = B[ipiv[i]] - sum;
            sum = 0;
        }
    }
    //backward substitution for upper triangular
    if (UPLO == 'U'){
        int i;
        int a, b;
        double sum;
        double x[n];
        double* y; 
        y = (double*) malloc (n * sizeof(double));
        for (b=0 ; b<n ; b++) y[b] = x[b];
        x[n-1] = y[n-1] / A[(n-1)*n+n-1];
        for (i=n-1 ; i>=0 ; i--){
            for (a=i+1 ; a < n ; a++){
                sum += x[a] * A[i*n+a];
            }
            x[i] = (y[i] - sum) / A[i*n+i];
            sum = 0;
        }
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
    for (i = 0; i < n; i += 3) {
        for (j = 0; j < n; j += 3) {
            register int t   = i*n+j; 
            register int tt  = t+n; 
            register int ttt = tt+n; 
            register double x00 = C[t];
            register double x01 = C[t+1];
            register double x02 = C[t+2];
            register double x10 = C[tt];
            register double x11 = C[tt+1];
            register double x12 = C[tt+2];
            register double x20 = C[ttt];
            register double x21 = C[ttt+1];
            register double x22 = C[ttt+2];
            for (k = 0; k < n; k += 3) {
                register int ta   = i*n+k;
                register int tta  = ta+n;
                register int ttta = tta+n;
                register int tb   = k*n+j;
                register int ttb  = tb+n;
                register int tttb = ttb+n;

                register double R1 = A[ta]; 
                register double R2 = A[tta]; 
                register double R3 = A[ttta];
                register double R4 = B[tb]; 
                register double R5 = B[tb+1]; 
                register double R6 = B[tb+2]; 

                x00 += R1 * R4;
                x01 += R1 * R5;
                x02 += R1 * R6;
                x10 += R2 * R4;
                x11 += R2 * R5;
                x12 += R2 * R6;
                x20 += R3 * R4;
                x21 += R3 * R5;
                x22 += R3 * R6;

                R1 = A[ta+1];
                R2 = A[tta+1];
                R3 = A[ttta+1];
                R4 = B[ttb];
                R5 = B[ttb+1];
                R6 = B[ttb+2];

                x00 += R1 * R4;
                x01 += R1 * R5;
                x02 += R1 * R6;
                x10 += R2 * R4;
                x11 += R2 * R5;
                x12 += R2 * R6;
                x20 += R3 * R4;
                x21 += R3 * R5;
                x22 += R3 * R6;

                R1 = A[ta+2];
                R2 = A[tta+2];
                R3 = A[ttta+2];
                R4 = B[tttb];
                R5 = B[tttb+1];
                R6 = B[tttb+2];

                x00 += R1 * R4;
                x01 += R1 * R5;
                x02 += R1 * R6;
                x10 += R2 * R4;
                x11 += R2 * R5;
                x12 += R2 * R6;
                x20 += R3 * R4;
                x21 += R3 * R5;
                x22 += R3 * R6;
            
            }
            C[t]     = x00;
            C[t+1]   = x01;
            C[t+2]   = x02;
            C[tt]    = x10;
            C[tt+1]  = x11;
            C[tt+2]  = x12;
            C[ttt]   = x20;
            C[ttt+1] = x21;
            C[ttt+2] = x22;
        }
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
//     int ib, end;
//     for ( ib = 1 ; ib <=n-1 ; ib += b){
//         end = ib + b-1              
//         //apply BLAS2 version of GEPP to  get A(ib:n , ib:end) = P’ * L’ * U’
//         int maxind;
//         double max;
//         int temps;
//         double tempv[n];
//         int i, t, j, k;
//         int a, b, c;
//         for (i = 1 ; i < n-1 ; i++){
//             maxind = i; 
//             max = abs(A[i*n+i]); 
//             for (t = i+1 ; t <= n ; t++){
//                 if (abs(A[t*n+i]) > max ) {
//                     maxind = t; 
//                     max = abs(A[t*n+i]); 
//                 }
//             }
//         }
//         if (max==0) 
//             return -1; 
//         else 
//             if (maxind != i ) {
//                 temps = ipiv[i]; 
//                 ipiv[i] = ipiv[maxind]; 
//                 ipiv[maxind] = temps;
//                 for (a=0 ; a <=n ; a++) tempv[a] = A[i*n+a]; 
//                 for (b=0 ; b <=n ; b++) A[i*n+b] = A[maxind*n+b]; 
//                 for (c=0 ; c <=n ; c++) A[maxind*n+c] = tempv[c];
//             }
        
//         for (j = i+1 ; j <= n ; j++) {
//             A[j*n+i] = A[j*n+i]/A[i*n+i];
//             for (k = i+1 ; k <=n ; k++)
//             A[j*n+k] = A[j*n+k] - A[j*n+i] * A[i*n+k]; 
//         } 


//         //… let LL denote the strict lower triangular part of A(ib:end , ib:end) + I  //… update next b rows of U //… apply delayed updates with single matrix-multiply  //… with inner dimension b
        
//         A(ib:end , end+1:n) = LL-1 * A(ib:end , end+1:n)         
//         A(end+1:n , end+1:n ) = A(end+1:n , end+1:n ) - A(end+1:n , ib:end) * A(ib:end , end+1:n)                    
//     }     

    return 0;
}
