int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    int maxind;
    int temps;
    double max;
    double tempv[n];
    int i, t, j, k;
    int a, b, c;
    for (i = 1 ; i <= n-1 ; i++){
        maxind = i; 
        max = abs(A[i*n+i]); 
        for (t = i+1 ; t < n ; t++){
            if (abs(A[t*n+i]) > max ) {
                maxind = t; 
                max = abs(A[t*n+i]); 
            }
        }
        if (max==0) 
            return -1; 
        else if (maxind != i ) {
            temps = ipiv[i]; 
            ipiv[i] = ipiv[maxind]; 
            ipiv[maxind] = temps;
            for (a=0 ; a <n ; a++) tempv[a] = A[i*n+a]; 
            for (b=0 ; b <n ; b++) A[i*n+b] = A[maxind*n+b]; 
            for (c=0 ; c <n ; c++) A[maxind*n+c] = tempv[c];
        }
        for (j = i+1 ; j < n ; j++) {
            A[j*n+i] = A[j*n+i]/A[i*n+i];
            for (k = i+1 ; k <n ; k++)
                A[j*n+k] = A[j*n+k] - A[j*n+i] * A[i*n+k]; 
        } 
    }
    return 0;
}
