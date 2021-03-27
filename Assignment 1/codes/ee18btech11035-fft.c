# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <complex.h>
# include <math.h>

double pi =  3.141592653589793238;

void elementwise_multiply(double complex* Y,double complex* X,double complex* Hz,int n)
{
	for(int i=0;i<n;i++)
	{
		Y[i] = Hz[i]*X[i];
	}
}

	
void fft(double complex* X, int n){
    if (n == 1)
        return;

   
    double complex *Xeven = (double complex *) malloc(n/2 * sizeof(double complex));
    double complex *Xodd = (double complex *) malloc(n/2 * sizeof(double complex));
    for (int i = 0;i < n/2;i++) 
    {
        Xeven[i] = X[2*i];
        Xodd[i] = X[(2*i)+1];
    }

    fft(Xeven, n/2);
    fft(Xodd, n/2);

    
    double complex eraised;
    for (int i = 0;i < n/2;i++) 
    {
		eraised = CMPLX(cos(2*pi*i/n),-sin(2*pi*i/n));
		
	    X[i] = Xeven[i] + eraised * Xodd[i];
        X[i+n/2] = Xeven[i] - eraised * Xodd[i];
       
    }
    free(Xeven);
    free(Xodd);
}

void ifft(double complex* X, int n){
    if (n == 1)
        return;

   
    double complex *Xeven = (double complex *) malloc(n/2 * sizeof(double complex));
    double complex *Xodd = (double complex *) malloc(n/2 * sizeof(double complex));
    for (int i = 0;i < n/2;i++) 
    {
        Xeven[i] = X[2*i];
        Xodd[i] = X[(2*i)+1];
    }

    ifft(Xeven, n/2);
    ifft(Xodd, n/2);

    
    double complex eraised;
    for (int i = 0;i < n/2;i++) 
    {
		eraised = CMPLX(cos(2*pi*i/n),sin(2*pi*i/n));
		
	    X[i] = Xeven[i] + eraised * Xodd[i];
        X[i+n/2] = Xeven[i] - eraised * Xodd[i];
       
    }
    free(Xeven);
    free(Xodd);
}



int main(){
    int n = pow(2,21);
    
    double *x = (double *) malloc(n *sizeof(double));
    double complex *Hz = (double complex *) malloc(n * sizeof(double complex)); 
    double complex *Y = (double complex *) malloc(n * sizeof(double complex));
    double complex *y = (double complex *) malloc(n * sizeof(double complex));
    double complex *X = (double complex *) malloc(n * sizeof(double complex));
    FILE *fi1, *fo1,*fi2,*fo2,*fo3;

    fi1 = fopen("x_fft.dat", "r");
    
    /*int count = 0;*/
    for(int i=0;i<n;i++) 
    {
        double data;
        fscanf(fi1, "%lf", &data);
        x[i] = data;
    }

    fclose(fi1);
    for(int i = 0; i < n; i++)
    {
        X[i] = x[i];
    }

    fft(X, n);

    fo1  = fopen("X_fft.dat", "w");
    
    for(int i = 0; i < n; i++)
    {
        fprintf(fo1, "%lf+j%lf\n", creal(X[i]), cimag(X[i]));
    }
    
    fclose(fo1);
    
    fi2 = fopen("H_z.dat","r");
    for(int i=0;i<n;i++) 
    {
        double re,im;
        fscanf(fi2, "%lf" "%lf", &re,&im);
        Hz[i] = CMPLX(re,im);
    }
    fclose(fi2);
    elementwise_multiply(Y,X,Hz,n);
    
    fo2 = fopen("Y_fft.dat","w");
    for(int i = 0; i < n; i++)
    {
        fprintf(fo2, "%lf+j%lf\n", creal(Y[i]), cimag(Y[i]));
    }
    fclose(fo2);
    
    for(int i=0;i<n;i++)
    {
		y[i] = Y[i];
	}
	
    ifft(y,n);
    
    fo3 = fopen("y_fft.dat","w");
    for(int i=0;i<n;i++)
    {
		y[i] = y[i]/n;
		fprintf(fo3, "%lf %lf\n", creal(y[i]), cimag(y[i]));
	}
	fclose(fo3);

    return 0;
}
