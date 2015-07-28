/*
//  C implementation of Markovian regimes switching Kalman filter
//  as specified by Kim and Nelson, extended to handle missing data.
//
//  author: Raoul Grasman
//  license: GNU GPL v2.0 (or higher at your choice)
//  copyright (c) 2008 Raoul Grasman
//
*/
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include <Rdefines.h>
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

#define alpha_x_At(alpha,A,p,q, C)         \
    for(k=0;k < (p); k++)                  \
    for(m=0;m < (q); m++)                  \
        (C)[m][k] = (alpha) * (A)[k][m];

#define a_alpha_x_At(alpha,A,p,q, C)       \
    for(k=0;k < (p); k++)                  \
    for(m=0;m < (q); m++)                  \
        (C)[m+(p)*k] = (alpha) * (A)[k+(p)*m];

#define a_alpha_x_A(alpha,A,p,q, C)        \
    for(k=0;k < (p); k++)                  \
    for(m=0;m < (q); m++)                  \
        (C)[m+(p)*k] = (alpha) * (A)[m+(p)*k];

#define matmulAB(A,p,q,B,r,C)              \
    for(k=0;k < (p);k++)                   \
    for(m=0;m < (r);m++) {                 \
        (C)[k][m] = 0.0;                   \
        for(n=0;n < q;n++)                 \
            (C)[k][m] += (A)[k][n] * (B)[n][m];  \
    }

#define a_matmulAB(A,p,q,B,r,C)            \
    for(k=0;k < (p);k++)                   \
    for(m=0;m < (r);m++) {                 \
        (C)[k+(p)*m] = 0.0;                \
        for(n=0;n < q;n++)                 \
            (C)[k+(p)*m] += (A)[k+(p)*n] * (B)[n+(q)*m];  \
    }

#define matmulABt(A,p,q,B,r,C)             \
    for(k=0;k < p;k++)                     \
    for(m=0;m < q;m++) {                   \
        (C)[k][m] = 0.0;                     \
        for(n=0;n < r;n++)                 \
            (C)[k][m] += (A)[k][n] * (B)[m][n];  \
    }

#define a_matmulABt(A,p,r,B,q,C)             \
    for(k=0;k < (p);k++)                     \
    for(m=0;m < (q);m++) {                   \
        (C)[k+(p)*m] = 0.0;                     \
        for(n=0;n < (r);n++)                 \
            (C)[k+(p)*m] += (A)[k+(p)*n] * (B)[m+(q)*n];  \
    }

#define matmulAtB(A,p,q,B,r,C)             \
    for(k=0;k < (q);k++)                   \
    for(m=0;m < (r);m++) {                 \
        C[k][m] = 0.0;                     \
        for(n=0;n < (p);n++)               \
            C[k][m] += A[n][k] * B[n][m];  \
    }

#define a_matmulAtB(A,p,q,B,r,C)                    \
    for(k=0;k < (q);k++)                            \
    for(m=0;m < (r);m++) {                          \
        C[k+(q)*m] = 0.0;                           \
        for(n=0;n < (p);n++)                        \
            C[k+(q)*m] += A[n+(p)*k] * B[n+(p)*m];  \
    }

#define matmulABAt(A,p,q,B,C)                          \
    for(k=0;k < (p);k++)                               \
    for(l=0;l < (p);l++) {                             \
        C[k][l] = 0.0;                                 \
        for(m=0;m < (q);m++)                           \
        for(n=0;n < (q);n++)                           \
            C[k][l] += A[k][m] * B[m][n] * A[l][n];    \
    }

#define a_matmulABAt(A,p,q,B,C)                          \
    for(k=0;k < (p);k++)                               \
    for(l=0;l < (p);l++) {                             \
        (C)[k+(p)*l] = 0.0;                                 \
        for(m=0;m < (q);m++)                           \
        for(n=0;n < (q);n++)                           \
            (C)[k+(p)*l] += (A)[k+(p)*m] * (B)[m+(q)*n] * (A)[l+(p)*n];    \
    }

#define matmulASAt(A,p,q,B,C)                          \
    for(k=0;k < (p);k++)                               \
    for(l=k;l < (p);l++) {                             \
        C[k][l] = 0.0;                                 \
        for(m=0;m < (q);m++)                           \
        for(n=0;n < (q);n++)                           \
            C[k][l] += A[k][m] * B[m][n] * A[l][n];    \
        C[l][k] = C[k][l];                             \
    }

#define a_matmulASAt(A,p,q,B,C)                        \
    for(k=0;k < (p);k++)                               \
    for(l=k;l < (p);l++) {                             \
        (C)[k+(p)*l] = 0.0;                            \
        for(m=0;m < (q);m++)                           \
        for(n=0;n < (q);n++)                           \
            (C)[k+(p)*l] += (A)[k+(p)*m] * (B)[m+(q)*n] * (A)[l+(p)*n]; \
        (C)[l+(p)*k] = (C)[k+(p)*l];                             \
    }

#define mataddAB(A,p,q,B,C)                            \
    for(k=0;k < (p);k++)                               \
    for(m=0;m < (q);m++)                               \
        (C)[k][l] = (A)[k][m] + (B)[k][m];

#define a_mataddAB(A,p,q,B,C)                            \
    for(k=0;k < (p)*(q);k++)                               \
        (C)[k] = (A)[k] + (B)[k];

#define mataddABt(A,p,q,B,C)                           \
    for(k=0;k < (p);k++)                                 \
    for(m=0;m < (q);m++)                                 \
        (C)[k][m] = (A)[k][m] + (B)[m][k];

#define a_mataddABt(A,p,q,B,C)                           \
    for(k=0;k < (p);k++)                                 \
    for(m=0;m < (q);m++)                                 \
        (C)[k+(p)*m] = (A)[k+(p)*m] + (B)[m+(q)*k];

#define x_plus_alpha_multAb(x,p,alpha,A,q,b,y)         \
    for(m=0;m < (p);m++) {                             \
        (y)[m] = 0.0;                                  \
        for(n=0;n < (q);n++)                           \
            (y)[m] += (A)[m][n] * (b)[n];              \
        (y)[m] *= (alpha);                             \
        (y)[m] += (x)[m];                              \
    }

#define a__x_plus_alpha_multAb(x,p,alpha,A,q,b,y)         \
    for(m=0;m < (p);m++) {                             \
        (y)[m] = 0.0;                                  \
        for(n=0;n < (q);n++)                           \
            (y)[m] += (A)[m+(p)*n] * (b)[n];              \
        (y)[m] *= (alpha);                             \
        (y)[m] += (x)[m];                              \
    }

#define X_plus_alpha_matmulAB(X,p,q,alpha,A,r,B, C)    \
    for(k=0;k < p;k++)                                 \
    for(l=0;l < q;l++) {                               \
        (C)[k+(p)*l] = 0.0;                                 \
        for(m=0;m < q;m++)                             \
            (C)[k+(p)*l] += (A)[k+(p)*m] * (B)[m+(r)*l];              \
        (C)[k+(p)*l] *= (alpha);                              \
        (C)[k+(p)*l] += (X)[k+(p)*l];                            \
    }

#define a__X_plus_alpha_matmulAB(X,p,q,alpha,A,r,B, C)    \
    for(k=0;k < (p);k++)                                  \
    for(l=0;l < (q);l++) {                                \
        (C)[k+(p)*l] = 0.0;                               \
        for(m=0;m < (r);m++)                              \
            (C)[k+(p)*l] += (A)[k+(p)*m] * (B)[m+(r)*l];  \
        (C)[k+(p)*l] *= (alpha);                          \
        (C)[k+(p)*l] += (X)[k+(p)*l];                     \
    }

#define a__X_plus_alpha_matmulABt(X,p,q,alpha,A,r,B, C)    \
    for(k=0;k < (p);k++)                                   \
    for(l=0;l < (q);l++) {                                 \
        (C)[k+(p)*l] = 0.0;                                \
        for(m=0;m < (r);m++)                               \
            (C)[k+(p)*l] += (A)[k+(p)*m] * (B)[l+(q)*m];   \
        (C)[k+(p)*l] *= (alpha);                           \
        (C)[k+(p)*l] += (X)[k+(p)*l];                      \
    }

#define new_matrix(A,p,q)       \
 (A) = Calloc((p), double*);    \
 for(m=0;m < (q); m++)          \
  (A)[m] = Calloc((q), double);

#define new_vector(A,p)       \
 (A) = Calloc((p), double);

#define a_new_matrix(A,p,q)       \
 (A) = Calloc((p)*(q), double);    \

#define free_matrix(A,p)        \
 for(m=0;m < (p); m++)          \
  Free(A[m]);                   \
  Free(A);

#define a_free_matrix(A,p)        \
  Free(A);

#define vec2mat(a,p,q,A)        \
 for(k=0;k < (p); k++)          \
 for(m=0;m < (q); m++)          \
   (A)[k][m] = (a)[k + (p)*m];

#define mat2vec(A,p,q,a)        \
 for(k=0;k < (p); k++)          \
 for(m=0;m < (q); m++)          \
   (a)[k + (p)*m] = (A)[k][m];


#define Rprintf_matrix(format,A,p,q)             \
  {                                              \
	for(k=0;k < (p); k++) {                        \
  for(m=0;m < (q); m++)                     		 \
      Rprintf((format), k, m, (A)[k+(p)*m]);     \
	    Rprintf("\n");                             \
	}                                              \
  }

#define __ 0
#define err(nme,idx,ival) (stop("oob %s, %s=%d",nme,idx,ival))
#define _H(i,j,k) (H[(i<ne?i:Err("H",1,i)) + ne*(j<ne?j:Err("H",2,j)) + ne*ne*(k<nm?k:Err("H",3,k))])
#define _K(i,j,k) (K[(i<ne?i:Err("K",1,i)) + ne*(j<ne?j:Err("K",2,j)) + ne*ne*(k<nm?k:Err("K",3,k))])
#define _G(i,j,k) (G[(i<ne?i:Err("G",1,i)) + ne*(j<ne?j:Err("G",2,j)) + ne*ne*(k<nm?k:Err("G",3,k))])
#define _B(i,j,k) (B[(i<ny?i:Err("B",1,i)) + ny*(j<nx?j:Err("B",2,j)) + ny*nx*(k<nm?k:Err("B",3,k))])
#define _c(i,j)   (c[(i<ne?i:Err("c",1,i)) + ne*(j<nm?j:Err("c",2,j))])
#define _R(i,j,k) (R[(i<ny?i:Err("R",1,i)) + ny*(j<ny?j:Err("R",2,j)) + ny*ny*(k<nm?k:Err("R",3,k))])
#define _W(i,j,k) (W[(i<ny?i:Err("W",1,i)) + ny*(j<ne?j:Err("W",2,j)) + ny*ne*(k<nm?k:Err("W",3,k))])
#define _p(i,j)   (p[(i<nm?i:Err("p",1,i)) + nm*(j<nm?j:Err("p",2,j))])
  // H = ne x ne x nm
  // K = ne x ne x nm
  // G = ne x ne x nm
  // B = ny x nx x nm
  // c = ne x nm
  // R = ny x ny x nm
  // W = ny x ne x nm
  // p = nm x nm

#define __v(i,j,k,l)    (   v[(i<ny?i:Err("v",   1,i))   + ny*(j<nm  ?j:Err("v",   2,j))   + ny*nm*(k<nm  ?k:Err("v",   3,k))   + ny*nm*nm*(l<nt  ?l:Err("v", 4,l))])
#define __F(i,j,k,l,m)  (   F[(i<ny?i:Err("F",   1,i))   + ny*(j<ny  ?j:Err("F",   2,j))   + ny*ny*(k<nm  ?k:Err("F",   3,k))   + ny*ny*nm*(l<nm  ?l:Err("F", 4,l)) + ny*ny*nm*nm*(m<nt?m:Err("F",5,m))])
#define __pr(i,j)       (  pr[(i<nm?i:Err("pr",  1,i))   + nm*(j<nt+1?j:Err("pr",  2,j))])
#define __tpr(i,j,k)    ( tpr[(i<nm?i:Err("tpr", 1,i))   + nm*(j<nm  ?j:Err("tpr", 2,j))   + nm*nm*(k<nt  ?k:Err("tpr", 3,k))])
#define __utpr(i,j,k)   (utpr[(i<nm?i:Err("utpr",1,i))   + nm*(j<nm  ?j:Err("utpr",2,j))   + nm*nm*(k<nt  ?k:Err("utpr",3,k))])
#define __jd(i,j,k)     (  jd[(i<nm?i:Err("jd",  1,i))   + nm*(j<nm  ?j:Err("jd",  2,j))   + nm*nm*(k<nt  ?k:Err("jd",  3,k))])
#define __md(i)         (  md[(i<nt?i:Err("md",  1,i))])
#define __auc(i,j,k)    ( auc[(i<ne?i:Err("auc", 1,i))   + ne*(j<nm  ?j:Err("auc", 2,j))   + ne*nm*(k<nt+1?k:Err("auc", 3,k))])
#define __au(i,j,k,l)   (  au[(i<ne?i:Err("au",  1,i))   + ne*(j<nm  ?j:Err("au",  2,j))   + ne*nm*(k<nm  ?k:Err("au",  3,k))  + ne*nm*nm*(l<nt  ?l:Err("au", 4,l))])
#define __ap(i,j,k,l)   (  ap[(i<ne?i:Err("ap",  1,i))   + ne*(j<nm  ?j:Err("ap",  2,j))   + ne*nm*(k<nm  ?k:Err("ap",  3,k))  + ne*ne*nm*(l<nt  ?l:Err("ap", 4,l))])
#define __Puc(i,j,k,l)  ( Puc[(i<ne?i:Err("Puc", 1,i))   + ne*(j<ne  ?j:Err("Puc", 2,j))   + ne*ne*(k<nm  ?k:Err("Puc", 3,k))  + ne*ne*nm*(l<nt+1?l:Err("Puc",4,l))])
#define __Pu(i,j,k,l,m) (  Pu[(i<ne?i:Err("Pu",  1,i))   + ne*(j<ne  ?j:Err("Pu",  2,j))   + ne*ne*(k<nm  ?k:Err("Pu",  3,k))  + ne*ne*nm*(l<nm  ?l:Err("Pu", 4,l)) + ne*ne*nm*nm*(m<nt?m:Err("Pu",5,m))])
#define __Pp(i,j,k,l,m) (  Pp[(i<ne?i:Err("Pp",  1,i))   + ne*(j<ne  ?j:Err("Pp",  2,j))   + ne*ne*(k<nm  ?k:Err("Pp",  3,k))  + ne*ne*nm*(l<nm  ?l:Err("Pp", 4,l)) + ne*ne*nm*nm*(m<nt?m:Err("Pp",5,m))])
#define __y(i,j)        (   y[(i<ny?i:Err("y",   1,i))   + ny*(j<nt  ?j:Err("y",   2,j))])
#define __x(i,j)        (   x[(i<nx?i:Err("x",   1,i))   + nx*(j<nt  ?j:Err("x",   2,j))])
#define __v1(i)         (  v1[(i<ny?i:Err("v1",  1,i))])

  // v =ny x nm x nm x nt
  // F = ny x ny x nm x nm x nt
  // pr = nm x (nt +1)
  // tpr = nm x nm x nt
  // utpr = nm x nm x nt
  // jd = nm x nm x nt
  // md = nt
  // auc = ne x nm x (nt+1)
  // au = ne x nm x nm x nt
  // ap = ne x nm x nm x nt
  // Puc = ne x ne x nm x (nt+1)
  // Pu = ne x ne x nm x nm x nt
  // Pp = ne x ne x nm x nm x nt
  // y = ny x nt (?)
  // x = nx x nt (?)


int Err(const char *nme, int idx, int ival){
	error("\nout of bound: %s, index %d = %d", nme, idx,ival);
	return 0;
}

/*** The following functions can be used to check the
 * correctness of the matrix multiplication macros.
 * The symbols are loaded in R; a set of R access
 * routines by the same name are defined in the mskf
 * package (not to be used by the user).
 */

void addmat(int *nr, int *nc, double *A, double *B, double *C){
	int k;
	a_mataddAB(A, *nr, *nc, B, C);
}
void mmatmulABAt(int *nr, int *nc, double *X, double *Y, double *Z){
	int k, l, m, n;
	a_matmulABAt(X, *nr, *nc, Y, Z);
}
void mmatmulASAt(int *nr, int *nc, double *X, double *Y, double *Z){
	int k, l, m, n;
	a_matmulASAt(X, *nr, *nc, Y, Z);
}
void mmatmulAB(int *nrx, int *nc, double *X, double *Y, int *ncy, double *Z){
	int k, m, n;
	a_matmulAB(X, *nrx, *nc, Y, *ncy, Z);
}
void mmatmulABt(int *nrx, int *nc, double *X, double *Y, int *nry, double *Z){
	int k, m, n;
	a_matmulABt(X, *nrx, *nc, Y, *nry, Z);
}

void mmatmulAtB(int *nrx, int *ncx, double *X, double *Y, int *ncy, double *Z){
	int k, m, n;
	a_matmulAtB(X, *nrx, *ncx, Y, *ncy, Z);
}

void mmataddAB(int *nrx, int *ncx, double *X, double *Y, double *Z){
	int k;
	a_mataddAB(X, *nrx, *ncx, Y, Z);
}
void mmataddABt(int *nrx, int *ncx, double *X, double *Y, double *Z){
	int k, m;
	a_mataddABt(X, *nrx, *ncx, Y, Z);
}

void mx_plus_alpha_multAb(double *x, int *nr, double *alpha, double *A, int *nc, double *b, double *y){
	int m, n;
	a__x_plus_alpha_multAb(x, *nr, *alpha, A, *nc, b, y);
}
void mX_plus_alpha_multAB(double *X, int *nr, int *nc, double *alpha, double *A, int *nca, double *B, double *Y){
	int k, l, m;
	a__X_plus_alpha_matmulAB(X, *nr, *nc, *alpha, A, *nca, B, Y);
}
void mX_plus_alpha_multABt(double *X, int *nr, int *nc, double *alpha, double *A, int *nca, double *B, double *Y){
	int k, l, m;
	a__X_plus_alpha_matmulABt(X, *nr, *nc, *alpha, A, *nca, B, Y);
}

double inv(double **help, int n, double eps)
{
    double deps, det, t;
    int j, k, l;
    if(n<=0) error("Invalid input in inv.");
    deps = (eps? eps: 1.49e-8);
    det = 1.0;
    for(j=0; j<n; j++){
        t = help[j][j];
        det *= t;
        help[j][j] = 1.0;
        for(k=0; k<n; k++){
            if(fabs(t)<deps)
                return det = 0.0;
            help[j][k] = help[j][k] / t;
        }
        for(k=0; k<n; k++){
            if(k==j) continue;
            t = help[k][j];
            help[k][j] = 0.0;
            for(l=0; l<n; l++)
                help[k][l] = help[k][l] - help[j][l]*t;
        }
    }
    return det;
}

void cholesk(const int n, double **a, const int ia, double *p, double *Det, int ifail)
{
	double tmp, det;
	ifail=1;
	if (a[0][0]<=0.0) return;
	p[0]=1.0/sqrt(a[0][0]);
	for(int i=1;i<n;i++)
		a[0][i]=p[0]*a[0][i];

	for(int i=2; i<=n; i++){
		tmp=0.0;
		for(int k=0; k<i-1;k++)
			tmp=tmp+a[k][i-1]*a[k][i-1];
		if (tmp>=a[i-1][i-1]) return;
		p[i-1]=1.0/sqrt(a[i-1][i-1]-tmp);
		for(int j=i; j<=n; j++){
			if (i==j) continue;
			a[i-1][j-1]=p[i-1];
			tmp=0.0;
			for(int k=0; k<i-1; k++)
				tmp=tmp+a[k][j-1]*a[k][i-1];
			a[i-1][j-1]=a[i-1][j-1]*(a[j-1][i-1]-tmp);
		}
	}
	det=1.0;
	for(int i=1; i<=n; i++){
		tmp = 1.0/p[i-1];
		det = det*tmp*tmp;
	}
	*Det = det;
	ifail=0;
	return;
}

void cholinv(double **a,const int ia, const int n, double *Det, int ifail, double *p)
{
	double det;
	ifail=0;
	det=0.0;
	cholesk(n, a, ia, p, &det, ifail);
	if (ifail != 0) return;
	for(int i=1; i<=n; i++){
		for(int j=1; j<=i; j++)
			a[i-1][j-1]=0.0;
		a[i-1][i-1]=p[i-1];
	}
	for(int j=1; j<=n-1; j++)
		for(int i=j+1; i<=n; i++){
			a[i-1][j-1]=0.0;
			for(int k=j; k <= i-1; k++){
				a[i-1][j-1]=a[i-1][j-1] - p[i-1]*(a[k-1][i-1]*a[k-1][j-1]);
		}
	}
	for(int j=1;j<=n-1;j++)
	for(int i=j+1; i<=n;i++){
		a[j-1][i-1]=0.0;
		for(int k=i; k<=n; k++)
			a[j-1][i-1] = a[j-1][i-1] + a[k-1][i-1] * a[k-1][j-1];
	}
	for(int i=1; i<=n; i++){
		p[i-1]=0.0;
		for(int k=i; k<=n; k++)
			p[i-1]=p[i-1]+a[k-1][i-1]*a[k-1][i-1];
		a[i-1][i-1]=p[i-1];
	}
	for(int i=1; i<=n; i++)
		for(int j=1; j<=i ; j++)
			a[i-1][j-1]=a[j-1][i-1];
	*Det = det;
	return;
}


void a_inv(int *n, double *A, double *eps, double *det)
{
	int k, m, ifail=0;
	double **help, *p;
	new_matrix(help, *n, *n);
	new_vector(p, *n);
	vec2mat(A, *n, *n, help);
	//(*det) = inv(help, *n, *eps);
	cholinv(help, *n, *n, det, ifail, p);
	mat2vec(help, *n, *n, A);
    free_matrix(help, *n);
	a_free_matrix(p, *n);
}

double _inv(int *n, double **A)
{
	int ifail=0;
	double *p, det;
	new_vector(p, *n);
	cholinv(A, *n, *n, &det, ifail, p);
	a_free_matrix(p, *n);
	return det;
}



int bla(int *n){
	Rprintf("binnen\n");
	return 1234;
}


void kfilter_timeloop(
		int *NT, int *NM, int *NE, int *NY, int *NX,
		double *c, double *H, double *B, double *K, double *G, double *R, double *W,
		double *v, double *F, double *pr, double *tpr, double *utpr, double *jd,
		double *md, double *auc, double *au, double *ap, double *p,
		double *Puc, double *Pu, double *Pp, double *y, double *x, double *Lp,
		int *Debug)
{
	(*Lp) = 1000.0;
	int nt = *NT, nm = *NM, ne = *NE, ny = *NY, nx = *NX;

	// declarations
	int t, i, j;
	int k, l, m, n;					// reserved for matrix calculation loops in macro's, do not use with these macro's!
	int debug = *Debug, missing = 0;
	double detF;



	/*** ----------------- debug code ------------------- ***/
	if(debug){
		//for(i=0;i<nx;i++)
		//	Rprintf("%d, ",x[i]); Rprintf("\n");
		Rprintf("\n=======================    INPUT    ===============================\n");
		Rprintf_matrix("x[%d,1] = %12.6g", &x[__], nx, 1);
		Rprintf("nt = %d\nne = %d\nnm = %d\nny = %d\nnx = %d\n", nt, ne, nm, ny, nx);
		Rprintf_matrix("p[%d,%d] = %g\t", &_p(__,__), nm, nm);
		for (i = 0; i < nm; i++) {
			Rprintf("\n\nRegime %d\n========\n\n\n\ta[t] = c + H a[t-1] + G u[t],     u[t] ~ N(0, K)\n\ty[t] = W a[t] + B x[t] + e[t],    e[t] ~ N(0, R)\n\n", i);
			Rprintf_matrix("H[%d,%d] = %12.6g\t", &_H(__, __, i), ne, ne); Rprintf("\n");
			Rprintf_matrix("K[%d,%d] = %12.6g\t", &_K(__, __, i), ne, ne); Rprintf("\n");
			Rprintf_matrix("G[%d,%d] = %12.6g\t", &_G(__, __, i), ne, ne); Rprintf("\n");
			Rprintf_matrix("c[%d,%d] = %12.6g\t", &_c(__, i),     ne,  1); Rprintf("\n");
			Rprintf_matrix("B[%d,%d] = %12.6g\t", &_B(__, __, i), ny, nx); Rprintf("\n");
			Rprintf_matrix("R[%d,%d] = %12.6g\t", &_R(__, __, i), ny, ny); Rprintf("\n");
			Rprintf_matrix("W[%d,%d] = %12.6g\t", &_W(__, __, i), ny, ne); Rprintf("\n");
		}
		Rprintf("-----------------------------------------------------------------------------------------------------------\n");
	}
	/*** ----------------- end  debug ------------------- ***/



	//
	// allocate temporary memory
	//

	double *HPucH, *GKG, *PpW, *PpWFinv, **_F, *Finv, *vtFinvv, *Ps; // matrices
	double  *v1, *as;           // vectors
	a_new_matrix(HPucH,   ne, ne);
	a_new_matrix(GKG,     ne, ne);
	a_new_matrix(PpW,     ne, ny);
	a_new_matrix(PpWFinv, ne, ny);
	new_matrix(_F,      ny, ny);
	a_new_matrix(Finv,    ny, ny);
	a_new_matrix(vtFinvv,  1, 1);
	a_new_matrix(Ps,      ne, ne);
	new_vector(v1,    ny);
	new_vector(as,    ne);


	/** ====================================  array sizes & definitions ================================= **/
	// H = ne x ne x nm
	// K = ne x ne x nm
	// G = ne x ne x nm
	// B = ny x nx x nm
	// c = ne x nm
	// R = ny x ny x nm
	// W = ny x ne x nm
	// p = nm x nm

	// v =ny x nm x nm x nt
	// F = ny x ny x nm x nm x nt
	// pr = nm x (nt +1)
	// tpr = nm x nm x nt
	// utpr = nm x nm x nt
	// jd = nm x nm x nt
	// md = nt
	// auc = ne x nm x (nt+1)
	// au = ne x nm x nm x nt
	// ap = ne x nm x nm x nt
	// Puc = ne x ne x nm x (nt+1)
	// Pu = ne x ne x nm x nm x nt
	// Pp = ne x ne x nm x nm x nt
	// y = ny x nt (?)
	// x = nx x nt (?)
	/** ================================================================================================= **/


	double jm = 0.0,                      /* joint marginal */
	L  = 0.0;                             /* log likelihood */

	for (t = 0; t < nt; t++) {
		/*** ----------------- debug code ------------------- ***/
		if(t >= debug) debug=0;
		if(debug){
			Rprintf("\nt = %d\n", t);
		}
		/*** ----------------- end debug ------------------- ***/

		// user interrupt
		R_CheckUserInterrupt(); // allow user to interrupt computations

		// //
		// test for missing values
		missing = 0;
		for (i = 0; i < ny; i++) {
			if(ISNA(y[i+ny*t])) // missing values indicated by NA
			{
				missing = 1;
				if(debug)
					Rprintf("missing value detected for t = %d\n", t);
				break;
			}
		}


		jm = 0.0;
		for(j=0; j<nm; j++)
			for(i=0; i<nm; i++)
			{
				///
				/// PREDICTION EQUATIONS
				///


				/*** ----------------- debug code ------------------- ***/
				if(debug && i==0 && j==0) {
					Rprintf_matrix("\t  p[%d,%d] = %12.6g", &_p(__, __),        nm, nm);  Rprintf("\n");
				}
				/*** ----------------- end debug ------------------- ***/


				// ap = c + H * auc
				a__x_plus_alpha_multAb(&_c(__, j), ne, 1.0, &_H(__, __,j), ne, &__auc(__, i, t), &__ap(__, i, j, t));                                     // H[][][j] = H[0+ne*0+ne*ne*j], auc[][i][t] = auc[0+ne*i+ne*nm*t], ap[][i][j][t] = ap[0+ne*i+ne*nm*j+ne*nm*nm*t]


				/*** ----------------- debug code ------------------- ***/
				if(debug) {
					Rprintf("ap = c + H %*% auc\n");
					Rprintf_matrix("\t  c[%d,%d] = %12.6g", &_c(__, j),         ne,  1);  Rprintf("\n");
					Rprintf_matrix("\t  H[%d,%d] = %12.6g", &_H(__, __, j),     ne, ne);  Rprintf("\n");
					Rprintf_matrix("\tauc[%d,%d] = %12.6g", &__auc(__, i, t),   ne,  1);  Rprintf("\n");
					Rprintf_matrix("\t ap[%d,%d] = %12.6g", &__ap(__, i, j, t), ne,  1);  Rprintf("\n");
				}
				/*** ----------------- end debug ------------------- ***/


				// HPucH = H * Puc * H'
				a_matmulABAt(&_H(__, __, j), ne, ne, &__Puc(__, __, i, t), HPucH );                                                             // H[][][j] = H[0+ne*0+ne*ne*j], Puc[][][i][t] = Puc[0+ne*0+ne*ne*i+ne*ne*nm*t]
				// GKG = G * K * G'
				a_matmulABAt(&_G(__, __, j), ne, ne, &_K(__, __, j), GKG);                                                                // G[][][j] = G[0 + ne*0 + ne*ne*j], K[][][j] = K[0 + ne*0 + ne*ne*j]
				// Pp = H * Puc * H'  + G * K * G'
				a_mataddAB(HPucH, ne, ne, GKG, &__Pp(__, __, i, j, t));


				/*** ----------------- debug code ------------------- ***/
				if(debug) Rprintf("\n S(i=%d) --> S(j=%d) at t = %d\t\t***********************************************\n\nPp = H Puc H' + G K G'", i, j, t);
				if(debug){
					Rprintf("\n");
					Rprintf_matrix("\tPuc[%d,%d] = %12.6g", &__Puc(__, __, j, t),    ne, ne); Rprintf("\n");
					Rprintf_matrix("\t  G[%d,%d] = %12.6g", &   _G(__, __, j),       ne, ne); Rprintf("\n");
					Rprintf_matrix("\t  K[%d,%d] = %12.6g", &   _K(__, __, j),       ne, ne); Rprintf("\n");
					Rprintf_matrix("\t Pp[%d,%d] = %12.6g", & __Pp(__, __, i, j, t), ne, ne);
				}
				/*** ----------------- end debug ------------------- ***/


				/// ONE-STEP-AHEAD PREDICTION ERROR
				// v1 = y - W * ap
				if(!missing) {
					a__x_plus_alpha_multAb(&__y(__, t), ny, -1.0, &_W(__, __, j), ne, &__ap(__, i, j, t), &__v1(__));                                                                 // y[][t] = y[0+ny*t], W[][][j] = W[0+ny*0+ny*ne*j], ap[][i][j][t] = ap[0+ne*i+ne*nm*j+ne*nm*nm*t],


					/*** ----------------- debug code ------------------- ***/
					if(debug){
						Rprintf("v1 = y - W ap\n");
						Rprintf_matrix("\t  y[%d,%d] = %12.6g",  & __y(__,  t),       ny,  1); Rprintf("\n");
						Rprintf_matrix("\t  W[%d,%d] = %12.6g",  &  _W(__, __, j),    ny, ne); Rprintf("\n");
						Rprintf_matrix("\t ap[%d,%d] = %12.6g", &__ap(__,  i, j, t),  ne,  1); Rprintf("\n");
						Rprintf_matrix("\t v1[%d,%d] = %12.6g", &__v1(__),            ny,  1);
					}
					/*** ----------------- end debug ------------------- ***/


					// v = y - W * ap - B * x = v1 - B * x
					a__x_plus_alpha_multAb(&__v1(__), ny, -1.0, &_B(__, __, j), nx, &__x(__, t), &__v(__, i, j, t));                                      // B[][][j] = B[0+ny*0+ny*nx*j], x[][t] = x[0+nx*t], v[][i][j][t] = v[0+ny*i+ny*nm*j+ny*nm*nm*t]


					/*** ----------------- debug code ------------------- ***/
					if(debug){
						Rprintf("v = v1 - B x\n");
						Rprintf_matrix("\t  B[%d,%d] = %12.6g", & _B(__, __, j),    ny, nx);  Rprintf("\n");
						Rprintf_matrix("\t  x[%d,%d] = %12.6g", &__x(__, t),        nx,  1);  Rprintf("\n");
						Rprintf_matrix("\t  v[%d,%d] = %12.6g", &__v(__,  i, j, t), ny,  1);
					}
					/*** ----------------- end debug ------------------- ***/
				}
				else {
					if(debug){
						Rprintf("\n\tmissing values for y ==> au = ap\n");
					}
				}


				/// COVARIANCE UPDATE
				// PpW = Pp * W'
				a_matmulABt(&__Pp(__, __, i, j, t), ne, ne, &_W(__, __, j), ny, PpW);                                                                // Pp[][][i][j][t] = Pp[0+ne*0+ne*ne*i+ne*ne*nm*j+ne*ne*nm*nm*t], W[][][j] = W[0+ny*0+ny*ne*j]
				// mmatmulABt(&ne, &ne, &__Pp(__, __, i, j, t), &_W(__, __, j), &ny, PpW);
				// F = W Pp W' + R
				a__X_plus_alpha_matmulAB(&_R(__, __, j), ny, ny, 1.0, &_W(__, __, j), ne, PpW, &__F(__, __, i, j, t));                        // R[][][j] = R[0+ny*0+ny*ny*j], W[][][j] = W[0+ny*0+ny*ne*j], F[][][i][j][t] = F[0+ny*0+ny*ny*i+ne*ne*nm*j+ne*ne*nm*nm*t]


				/*** ----------------- debug code ------------------- ***/
				if(debug){
					Rprintf("F = W Pp W' + R = W PpW + R\n");
					Rprintf_matrix("\t Pp[%d,%d] = %12.6g", &__Pp(__, __, i, j, t), ne, ne); Rprintf("\n");
					Rprintf_matrix("\tPpW[%d,%d] = %12.6g",   PpW,                  ne, ny); Rprintf("\n");
					//for (int ii=0; ii < ne*ny; ii++) {Rprintf("%12.6g ", PpW[ii]);}; Rprintf("\n");
					Rprintf_matrix("\t  W[%d,%d] = %12.6g", &  _W(__, __, j),       ny, ne); Rprintf("\n");
					Rprintf_matrix("\t  R[%d,%d] = %12.6g", &  _R(__, __, j),       ny, ny); Rprintf("\n");
					Rprintf_matrix("\t  F[%d,%d] = %12.6g", & __F(__, __, i, j, t), ny, ny); Rprintf("\n");
				}
				/*** ----------------- end debug ------------------- ***/


				///
				/// UPDATE EQUATIONS
				///


				if(missing) { // updated state with expected value under Martingale-condition
					for(m=0;m<ne;m++){
						__au(m, i, j, t) = __ap(m, i, j, t);
					}
				}
				else {
					// inverse and determinant of covariance
					vec2mat(&__F(__, __, i, j, t), ny, ny, _F); // double* -> double**
					detF = _inv(&ny, _F); // inverse of _F stored in _F and returns determinant of F
					mat2vec(_F, ny, ny, Finv); // double** -> double*

					/*** ----------------- debug code ------------------- ***/
					if(debug){
						Rprintf("\n");
						Rprintf_matrix("\tFinv[%d,%d] = %12.6g", &Finv[0], ny, ny); Rprintf("\n");
						Rprintf("\t|F| = %12.6g\n", detF);
					}
					/*** ----------------- end debug ------------------- ***/

					// PpWFinv = Pp W' Finv = PpW Finv
					a_matmulAB(PpW, ne, ny, Finv, ny, PpWFinv);
					// au = ap + Pp W Finv v
					a__x_plus_alpha_multAb(&__ap(__, i, j, t), ne, 1.0, PpWFinv, ny, &__v(__, i, j, t), &__au(__, i, j, t));                                    // ap[][i][j][t] = ap[0+ne*i+ne*nm*j+ne*nm*nm*t], v[][i][j][t] = v[0+ny*i+ny*nm*j+ny*nm*nm*t], au[][i][j][t] = au[0+ne*i+ne*nm*j+ne*nm*nm*t];
				}

				// Pu = Pp + Pp W Finv W' Pp // was Pp - Pp W Finv W' Pp in earlier version, original R implementation specifies +
				a__X_plus_alpha_matmulABt(&__Pp(__, __, i, j, t), ne, ne, 1.0, PpWFinv, ny, PpW, &__Pu(__, __, i, j, t));                      // Pp[][][i][j][t] = Pp[0+ne*0+ne*ne*i+ne*ne*nm*j+ne*ne*nm*nm*t], Pu[][][i][j][t] = Pu[0+ne*0+ne*ne*i+ne*ne*nm*j+ne*ne*nm*nm*t]

				/*** ----------------- debug code ------------------- ***/
				if(debug){
					double *myPu;
					a_new_matrix(myPu,   ne, ne);
					a__X_plus_alpha_matmulABt(&__Pp(__, __, i, j, t), ne, ne, 1.0, PpWFinv, ny, PpW, myPu);
					// Rprintf("\tPp = %g, PpWFinv = %g, PpW = %g, Pu = Pp + Pp * W' * Finv * W * Pp = %g\n", Pp[0+ne*0+ne*ne*i+ne*ne*nm*j+ne*ne*nm*nm*t], PpWFinv[0], PpW[0], Pu[0+ne*0+ne*ne*i+ne*ne*nm*j+ne*ne*nm*nm*t]);
					Rprintf("au = ap + Pp W Finv v = ap + PpW Finv v\n");
					Rprintf_matrix("\tPpWFinv[%d,%d] = %12.6g",  PpWFinv,               ne, ny); Rprintf("\n");
					Rprintf_matrix("\t     au[%d,%d] = %12.6g", &__au(__, i, j, t),     ne,  1); Rprintf("\n");
					Rprintf("Pu = Pp + Pp W Finv W' Pp\n"); // before this was Rprintf("Pu = Pp - Pp W Finv W' Pp\n");
					Rprintf_matrix("\t     Pu[%d,%d] = %12.6g", &__Pu(__, __, i, j, t), ne, ne); Rprintf("\n");
					Rprintf_matrix("\t   myPu[%d,%d] = %12.6g", myPu,                   ne, ne); Rprintf("\n");
					a_free_matrix(myPu, ne);
					Rprintf("freed myPu\n");
				}
				/*** ----------------- end debug ------------------- ***/


				///
				/// HAMILTON FILTER STEP
				///


				// ,,,,,,  begin hamilton (end is outsite this regime loop)  ,,,,,,


				// conditional transition probability S(t-1)=i --> S(t) = j
				__tpr(j, i, t) = _p(j, i) * __pr(i, t); //  tpr[j][i][t] = tpj[j+nm*i+nm*nm*t], p[j][i] = p[j+nm*i],  pr[i][t] = pr[i+nm*t]

				if(!missing){

					// conditional density
					a_matmulASAt(&__v(__, i, j, t), 1, ny, &Finv[0], vtFinvv);            // v[][i][j][t] = v[0+ny*i+ny*nm*j+ny*nm*nm*nt]
											/************
											 * MOET HIER BOVEN NIET &__v, 1 EN ny OMGEDRAAID WORDEN MET &Finv[0], EN 1, ny EIGENLIJK ny, 1 ZIJN????
											************/

					/*** ----------------- debug code ------------------- ***/
					if(debug)
						Rprintf("\n\t        v = %12.6g\n\tv' Finv v = %12.6g\n", __v(0, i, j, t), vtFinvv[0]);
					/*** ----------------- end debug ------------------- ***/


					__jd(j, i, t) = __tpr(j, i, t) * pow(M_2PI ,-((double)ny)/2.0) * exp(-0.5*vtFinvv[0]) / sqrt(detF);                          // jd[j][i][t] = jd[j+nm*i+nm*nm*t], tpr[j][i][t] = tpr[j+nm*i+nm*nm*t]


					/*** ----------------- debug code ------------------- ***/
					if(debug) {
						Rprintf("\njd[%d,%d] = tpr[%d,%d] * (2 * pi)^(ny/2) * exp(-0.5 * v' Finv v) / sqrt(det(F))", j, i, j, i);
						Rprintf("\n\t             p = %12.6g\n\t            pr = %12.6g\n\t           tpr = %12.6g\n\t(2*pi)^(-ny/2) = %12.6g\n\t    vtFinvv[0] = %12.6g\n\t          detF = %12.6g\n\t       jd[%d,%d] = %12.6g\n", _p(j, i), __pr(i, t), __tpr(j,i,t), pow(M_2PI,-(double)ny/2.0), vtFinvv[0], detF, j, i, __jd(j, i, t));
					}
					/*** ----------------- end debug ------------------- ***/


					// marginal density
					jm += __jd(j, i, t);                                                      // jd[j][i][t] = jd[j+nm*i+nm*nm*t]

				}
			}


			if(!missing)
				L += log(jm);


			/*** ----------------- debug code ------------------- ***/
			if(debug)
				Rprintf("\n\n              +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\t       missing = %12s\n\tlog-Likelihood = %12.6g\n\tjm             = %12.6g\n\n", (missing?"TRUE":"FALSE"), L, jm);
			/*** ----------------- end debug ------------------- ***/


			// transition probability update

			if(missing){
				// updated transition probabilities remain unchanged
				for( m=0; m < nm; m++)
					for( k=0; k < nm; k++)
						__utpr(m, k, t) = __tpr(k, m, t);
			}
			else {
				// utpr[i][j] = jd[j][i] / jm
				a_alpha_x_At(1.0/jm, &__jd(__, __, t), nm, nm, &__utpr(__, __, t));          // jd[][][t] = jd[0+nm*0+nm*nm*t], utpr[][][t] = utpr[0+nm*0+nm*nm*t] /* ik heb bij utpr de kolommen en rijen (i en j) omgedraait  */
			}


			/*** ----------------- debug code ------------------- ***/
			if(debug){
				Rprintf("utpr = jd * (1/jm)\n");
				Rprintf_matrix("\tutpr[%d,%d] = %12.6g", &__utpr(__, __, t), nm, nm); Rprintf("\n");
			}
			/*** ----------------- end debug ------------------- ***/



			// `````` end of hamilton ``````



			///
			/// COLLAPSE STEP (to reduce the number of posteriors)
			///

			for(j=0; j<nm; j++){ // loop over regimes

				// updated regime probabilities

				__pr(j,t+1) = 0.0;                                                          // pr[j][(t+1)] = pr[j+nm*(t+1)]
				for(i=0; i<nm; i++)
					__pr(j, t+1) += __utpr(i, j, t);                                    // pr[j][(t+1)] = pr[j+nm*(t+1)], utpr[i][j][t] = utpr[i+nm*j+nm*nm*t]


				/*** ----------------- debug code ------------------- ***/
				if(debug){
					Rprintf("\n");
					Rprintf_matrix("\t  pr[%d,%d] = %12.6g", &__pr(__, t+1), nm, 1); Rprintf("\n");
				}
				/*** ----------------- end debug ------------------- ***/


				for(k=0; k<ne; k++) {
					as[k] = 0.0;
					for(m=0; m<ne; m++)
						Ps[m+k*ne] = 0.0;
				}

				for(i=0; i<nm; i++){
					for(k=0; k<ne; k++)
						as[k] = as[k] + __au(k, i, j, t) * __utpr(i, j, t) / __pr(j, t+1);
				}

				/*** ----------------- debug code ------------------- ***/
				if(debug){
					Rprintf("as = au utpr / pr \n");
					Rprintf_matrix("\tas[%d,%d] = %12.6g", as, ne, 1); Rprintf("\n");
					Rprintf_matrix("\tau[%d,%d] = %12.6g", &__au(__, __, j, t), ne, nm);
				}
				/*** ----------------- end debug ------------------- ***/

				for(k=0; k<ne; k++)
					__auc(k, j, t+1) = as[k];

				for(i=0; i<nm; i++){

					for(k=0; k<ne; k++)
						for(m=0; m<ne; m++)
							Ps[m+k*ne] = Ps[m+k*ne]
						+ (__Pu(m, k, i, j, t)
           + (__auc(m, j, t+1)-__au(m, i, j, t)) * (__auc(k, j, t+1)-__au(k,i,j,t))) *
           __utpr(i, j, t) / __pr(j, t+1);

					/*** ----------------- debug code ------------------- ***/
					if(debug){
						Rprintf("\n\t auc(j,t+1) - au(i,j,t) = %12.6g", __auc(0, j, t+1)-__au(0, i, j, t));
						Rprintf("\n\tutpr(i,j,t) / pr(j,t+1) = %12.6g", __utpr(i, j, t) / __pr(j,t+1)); Rprintf("\n");
					}
					/*** ----------------- end debug ------------------- ***/
				}


				/*** ----------------- debug code ------------------- ***/
				if(debug){
					Rprintf("\n");
					Rprintf("\tj= %d\n", j);
					Rprintf_matrix("\tPu[%d,%d] = %12.6g", &__Pu(__,__, 0, j, t), ne, ne); Rprintf("\n");
					Rprintf_matrix("\tPs[%d,%d] = %12.6g", Ps, ne, ne);
				}
				/*** ----------------- end debug ------------------- ***/


				for(k=0; k<ne; k++)
					for(m=0; m<ne; m++)
						__Puc(m,k,j,t+1) = Ps[m+k*ne];

			} // end loop over regimes
	} // end of time loop
	if(*Debug) Rprintf("\n\t\t<><><><><><><><><><><><>        E N D   O F   L O O P      <><><><><><><><><><><><>\n");

	a_free_matrix(HPucH, ne);
	if(*Debug) Rprintf("freed HPucH\n");
	a_free_matrix(GKG,   ne);
	if(*Debug) Rprintf("freed GKG\n");
	a_free_matrix(PpW,   ne);
	if(*Debug) Rprintf("freed PpW\n");
	a_free_matrix(Finv,  ny);
	if(*Debug) Rprintf("freed Finv\n");

	a_free_matrix(PpWFinv, ne);
	if(*Debug) Rprintf("freed PpWFinv\n");
	a_free_matrix(vtFinvv,  1);
	if(*Debug) Rprintf("freed vtFinvv\n");
	a_free_matrix(Ps,    ne);
	if(*Debug) Rprintf("freed Ps\n");
	a_free_matrix(v1,    ny);
	if(*Debug) Rprintf("freed v1\n");
	a_free_matrix(as,    ne);
	if(*Debug) Rprintf("freed as\n");

	free_matrix(_F,    ny);
	if(*Debug) Rprintf("freed _F\n");
	(*Lp) = L;
	if(*Debug) Rprintf("Exiting kfilter_timeloop\n");
}



