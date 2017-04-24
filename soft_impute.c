#include <stdio.h>
#define M 10
#define N 10

void dgesvd_(char *JOBU, char *JOBVT, int *M, int *N, double *A, int *LDA, double *S, double *U, int *LDU, double *VT, int *LDVT, double *WORK, int *LWORK, int *INFO);

void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, double *A, int *LDA, double *B, int *LDB, double * BETA, double * C, int* LDC);

int min(int m, int n){
	if (m<n)
		return m;
	else
		return n;
}

void diag(double *s, int r, double *ds){
	int i,j;
	for (i = 0; i < r; i++){
		for (j = 0; j < r; j++){
			if (i == j)
				ds[i + j*r] = s[i];
			else
				ds[i + j*r] = 0;
		}
	}
}

double fnormsq(double *z, int m, int n){
	int i,j;
	double sum = 0;
	for (i=0; i< m; i++){
		for (j = 0; j < n; j++){
			sum = sum + z[i + j*m]*z[i + j*m];
		}
	}
	return sum;
}

double diff(double *zold, double *znew, int m, int n){
	double zold_minus_znew[m*n];
	int i,j;
	for (i = 0; i < m; i++){
		for (j = 0; j < m; j++){
			zold_minus_znew[i+j*m] = zold[i+j*m] - znew[i+j*m];
		}
	}
	return fnormsq(zold_minus_znew, m,n)/fnormsq(zold,m,n);
}

void soft_threshold(double *s, double lambda, int r){
	int i;
	for (i=0; i<r; i++){
		double s_minus_lambda = s[i] - lambda;
		if (s_minus_lambda < 0)
			s[i] = 0;
		else
			s[i] = s[i] - lambda;
	}
}

void P_omega_c(double *zold, double *P_omega_c_zold, double *omega, int m, int n){
	int i;
	for (i= 0; i < m*n; i++)
		P_omega_c_zold[i] = omega[i] * zold[i];
}

void get_znew(double *znew, double *z, double *zold, double *omega, double lambda, int m, int n){
	double P_Z_plus_P_c_Z_old[m*n];
	double P_c_Z_old[m*n];
	P_omega_c(zold, P_c_Z_old, omega, m, n);
	int i,j;
	for (i = 0; i<m; i++){
		for (j = 0; j<n; j++){
			P_Z_plus_P_c_Z_old[i+j*m] = z[i+j*m] + P_c_Z_old[i+j*m]; 
		}
	}

	int r;
	r = min(m,n);

	// svd
	char jobu = 'S';
    char jobvt = 'S';
    //int m = m;
    //int n = n;
    double *A = P_Z_plus_P_c_Z_old;
    double s[r];
    double u[m * r];
    double vt[r * n];
    double *work;
    int lwork = -1;
    double lworkopt;
    int info;

        dgesvd_(&jobu, &jobvt, &m, &n, A, &m, s, u, &m, vt, &n, &lworkopt, &lwork, &info);

    if (info != 0)
    {
        printf("The dgesvd error %d\n", info);
    }
    else
    {
        lwork = (int) lworkopt;
        work = (double *) malloc(lwork * sizeof(double));
        assert(work != NULL);

        dgesvd_(&jobu, &jobvt, &m, &n, A, &m, s, u, &m, vt, &n, work, &lwork, &info);

        if (info != 0)
        {
            printf("The dgesvd error %d\n", info);
        }
    }
    //svd end

    soft_threshold(s, lambda, r);

    double ds[r*r];
    diag(s,r,ds);

    double us[m*r];
    char transa = 'N', transb = 'N';
    double alpha = 1;
    double beta = 0;

    // U*S_{lambda}
    dgemm_(&transa, &transb, &m, &r, &r, &alpha, u, &m, ds, &r, &beta, us, &m);
    // U*S_{Lambda}*Vt
    dgemm_(&transa, &transb, &m, &n, &r, &alpha, us, &m, vt, &r, &beta, znew, &m);
}

void soft_impute(double *zhat, double *z, double *zold, double *omega, double lambda, int m, int n){
	get_znew(zhat, z, zold, omega, lambda, m, n);
	zold = zhat;
	while (1){
		get_znew(zhat, z, zold, omega, lambda, m, n);
		if (diff(zold, zhat, m, n) < 1e-5)
			break;
	}
}



int main(){
double Z_missing[M][N] =
{{0,-2.08,-1.67,-0.48,-0.52,-1.37,0.47,0.57,-0.2,-0.37},
{1.77,0,0.72,0,0,2.22,0.99,0,-2.54,3.76},
{-7.45,0,2.93,0,-0.68,0,-4.64,4.84,0,-7.13},
{3.26,0,-1.1,0,0,0,-1.34,-0.45,-2.14,0.73},
{3.46,0,0,-1.66,0,-0.83,3.65,-1.93,-1.8,3.19},
{-1.79,1.14,1.16,0,2.93,0,0,0,0,4.06},
{0,0.04,1.44,1.4,0,0,-2.14,1.52,0,-2.34},
{1.18,0.69,-0.49,0,0,0.45,4.33,-2.87,-1.26,0},
{0,0,0,-0.85,2.99,0,6.29,0,-1.53,6},
{2.16,0,0,-1.22,0,0,-1.05,0,0,0}};

double omega[M][N] = 
{{0,1,1,1,1,1,1,1,1,1},
{1,0,1,0,0,1,1,0,1,1},
{1,0,1,0,1,0,1,1,0,1},
{1,0,1,0,0,0,1,1,1,1},
{1,0,0,1,0,1,1,1,1,1},
{1,1,1,0,1,0,0,0,0,1},
{0,1,1,1,0,0,1,1,0,1},
{1,1,1,0,0,1,1,1,1,0},
{0,0,0,1,1,0,1,0,1,1},
{1,0,0,1,0,0,1,0,0,0}};

double Z_true[M][N] =
{{1.61,-2.08,-1.65,-0.48,-0.5,-1.36,0.49,0.58,-0.15,-0.42},
{1.79,4.1,0.69,-1.55,0.68,2.21,0.98,-2.93,-2.56,3.74},
{-7.43,-2.38,2.96,4.18,-0.71,-0.72,-4.64,4.82,5.29,-7.11},
{3.22,1.3,-1.05,-1.72,-1.19,0.38,-1.35,-0.43,-2.11,0.68},
{3.46,-1.02,-2.28,-1.68,0.67,-0.85,3.66,-1.94,-1.81,3.17},
{-1.76,1.18,1.16,0.47,2.92,1,5.04,-3.05,0.1,4.05},
{-2.66,0.05,1.44,1.38,-0.32,0.25,-2.13,1.51,1.6,-2.36},
{1.2,0.74,-0.47,-0.87,1.81,0.42,4.32,-2.87,-1.28,4.09},
{0.64,1.76,0.14,-0.84,3.01,1.12,6.3,-4.38,-1.52,6.05},
{2.13,1.53,-0.42,-1.23,-0.75,0.63,-1.06,-0.55,-1.62,0.73}};

double lambda[10];

// Let lambda be 1:100
for (i=0; i<10; i++){
	lambda[i] = i+1;
}

double z[M*N];
double zold[M*N];
double zhat[M*N];
//double res[M][N][10];
double error[10];
int i,j,k;

// Initialize zold as all zero
for (i=0; i<M; i++){
	for (j=0; j<N; j++){
		zold[i + M*j] = 0;
	}
}

// Initialize z as Z_missing, transpose for input
for (i=0; i<M; i++){
	for (j=0; j<N; j++){
		zold[i + M*j] = Z_missing[i][j];
	}
}


for (k = 0; k < 10; k++){
	soft_impute(zhat, z, zold, &(omega[0][0]), lambda[k], M, N);

	for (i = 0; i < M; i++){
		for (j = 0; j < M; j++){
			//res[i][j][k] = zhat[i + M*j];
			error[k] = diff(zhat, &(Z_true[0][0]), M, N);
			printf("%f\n", error[k]);
		}
		}
	

	for (i=0; i<M*N; i++){
		zold[i] = zhat[i];
	}
}





	

return 0;
}






