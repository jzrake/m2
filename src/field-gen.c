#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "jsw_rand.h"


typedef complex double Complex;

static double BESSJ(int N, double X);
//static void SPHJ(int N, double X, int *NM, double *SJ, double *DJ);

typedef struct fourier_mode
{
  double k[4];
  Complex A[4];
} fourier_mode;



double m2_force_free_vector_potential(double x[4], double n[4], int model, int k2)
{
#define RAND jsw_random_double(&rand, -1, 1)
  int m,i,j,k;
  jsw_rand_t rand;
  jsw_seed(&rand, model);


  int k2_sphere = k2;//26;
  int k_cube = floor(sqrt(k2_sphere)) + 1;
  fourier_mode *modes = NULL;
  int num_modes = 0;
  Complex A[4] = {0, 0, 0, 0};
  double amp;


  for (i=-k_cube; i<=k_cube; ++i) {
    for (j=-k_cube; j<=k_cube; ++j) {
      for (k=-k_cube; k<=k_cube; ++k) {
  	fourier_mode M;

  	if (i*i + j*j + k*k != k2_sphere) {
  	  continue;
  	}
  	else {
  	  M.k[0] = 0.0;
  	  M.k[1] = i;
  	  M.k[2] = j;
  	  M.k[3] = k;
  	  M.A[0] = 0.0;
  	  M.A[1] = RAND + RAND*I;
  	  M.A[2] = RAND + RAND*I;
  	  M.A[3] = RAND + RAND*I;
  	  amp = sqrt(M.A[1]*conj(M.A[1]) + M.A[2]*conj(M.A[2]) + M.A[3]*conj(M.A[3]));
  	  M.A[1] /= amp;
  	  M.A[2] /= amp;
  	  M.A[3] /= amp;
  	  num_modes += 1;
  	  modes = (fourier_mode *) realloc(modes, num_modes * sizeof(fourier_mode));
  	  modes[num_modes-1] = M;

  	  /* printf("k[%d] = [%d %d %d] is on shell\n", num_modes, i, j, k); */
  	}
      }
    }
  }


  for (m=0; m<num_modes; ++m) {
    fourier_mode M = modes[m];
    double a = sqrt(M.k[1]*M.k[1] + M.k[2]*M.k[2] + M.k[3]*M.k[3]);
    Complex K[4] = {0, I*M.k[1], I*M.k[2], I*M.k[3]};
    Complex Ikx  = (K[1]*x[1] + K[2]*x[2] + K[3]*x[3]) * M_PI;
    Complex P[4] = {0, M.A[1], M.A[2], M.A[3]}; /* a times psi */

    Complex T[4] = {0, /* T = K cross (a psi) */
		    K[2]*P[3] - K[3]*P[2],
		    K[3]*P[1] - K[1]*P[3],
		    K[1]*P[2] - K[2]*P[1]};

    Complex S[4] = {0, /* S = K cross T / alpha */
		    (K[2]*T[3] - K[3]*T[2])/a,
		    (K[3]*T[1] - K[1]*T[3])/a,
		    (K[1]*T[2] - K[2]*T[1])/a};

    A[1] += (S[1] + T[1])/a * cexp(Ikx);
    A[2] += (S[2] + T[2])/a * cexp(Ikx);
    A[3] += (S[3] + T[3])/a * cexp(Ikx);
  }


  free(modes);
  return creal(A[1]*n[1] + A[2]*n[2] + A[3]*n[3]);

#undef RAND
}


double m2_force_free_magnetic_field(double x[4], double n[4], int model, int k2)
{
  /* correct up to factor of alpha */
  return m2_force_free_vector_potential(x, n, model, k2);
}


double m2_spheromak_vector_potential(double x[4], double n[4])
{
#define c cos
#define s sin
  double r0 = 4.49340945791; /* magnetic surface */
  double Ar, At, Ap;

  double R = 1e-16 + sqrt(x[1]*x[1] + x[2]*x[2]);
  double r = 1e-16 + sqrt(x[1]*x[1] + x[2]*x[2] + x[3]*x[3]);
  double t = acos(x[3] / r);
  //double p = atan2(x[2], x[1]);
  double rhat[4] = {0, x[1]/r, x[2]/r, x[3]/r};
  double that[4] = {0, x[3]/r * x[1]/R, x[3]/r * x[2]/R, -R/r};
  double phat[4] = {0,-x[2]/R, x[1]/R, 0};

  if (r < r0) {
    Ar = 1.0 * pow(r,-3) * c(t) * (r*c(r) - s(r));
    At = 0.5 * pow(r,-3) * s(t) * (r*c(r) + s(r) * (r*r - 1));
    Ap = 0.5 * pow(r,-2) * s(t) * (r*c(r) - s(r));

  }
  else {
    double A = 2./3 * (r * cos(r) + (r*r - 1) * sin(r)) * pow(r, -3);
    double B = 1./3 * (r * cos(r) + (r*r - 1) * sin(r));

    Ar = -0.5 * (A - 2 * B * pow(r, -3)) * c(t);
    At = +0.5 * (A + 1 * B * pow(r, -3)) * s(t);
    Ap = +0.0;
  }

  double A[4];
  A[1] = Ar * rhat[1] + At * that[1] + Ap * phat[1];
  A[2] = Ar * rhat[2] + At * that[2] + Ap * phat[2];
  A[3] = Ar * rhat[3] + At * that[3] + Ap * phat[3];

  return A[1]*n[1] + A[2]*n[2] + A[3]*n[3];
#undef c
#undef x
}


double m2_tubomak_vector_potential(double x[4], double n[4])
{
  double a = 3.83170597021; /* first zero of A-phi will be at r=1 */
  double r = sqrt(x[1]*x[1] + x[2]*x[2]);
  if (fabs(r) < 1e-12) r = 1e-12;
  double Af =   a*BESSJ(1, a*r);
  double Az = ((a*BESSJ(1, a*r))/r + (pow(a,2)*(BESSJ(0, a*r) - BESSJ(2, a*r)))/2.)/a;
  double Ax = Af * (-x[2]/r);
  double Ay = Af * (+x[1]/r);
  if (r < 1) {
    return Ax*n[1] + Ay*n[2] + Az*n[3];
  }
  else {
    r = 1;
    double Az = ((a*BESSJ(1, a*r))/r + (pow(a,2)*(BESSJ(0, a*r) - BESSJ(2, a*r)))/2.)/a;
    return Az * n[3];
  }
}






/***********************************************************************
 *                                                                      *
 *    Program to calculate the first kind Bessel function of integer    *
 *    order N, for any REAL X, using the function BESSJ(N,X).           *
 *                                                                      *
 * -------------------------------------------------------------------- *
 *                                                                      *
 *    SAMPLE RUN:                                                       *
 *                                                                      *
 *    (Calculate Bessel function for N=2, X=0.75).                      *
 *                                                                      *
 *    Bessel function of order  2 for X =  0.7500:                      *
 *                                                                      *
 *         Y =  0.06707400                                              *
 *                                                                      *
 * -------------------------------------------------------------------- *
 *   Reference: From Numath Library By Tuan Dang Trong in Fortran 77.   *
 *                                                                      *
 *                               C++ Release 1.0 By J-P Moreau, Paris.  *
 *                                        (www.jpmoreau.fr)             *
 ***********************************************************************/
#include <math.h>
#include <stdio.h>

static double Sign(double X, double Y);
static double BESSJ0(double X);
static double BESSJ1(double X);

double Sign(double X, double Y) {
  if (Y<0.0) return (-fabs(X));
  else return (fabs(X));
}

double BESSJ0 (double X) {
  /***********************************************************************
      This subroutine calculates the First Kind Bessel Function of
      order 0, for any real number X. The polynomial approximation by
      series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
      REFERENCES:
      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
      VOL.5, 1962.
  ************************************************************************/
  const double
    P1=1.0, P2=-0.1098628627E-2, P3=0.2734510407E-4,
    P4=-0.2073370639E-5, P5= 0.2093887211E-6,
    Q1=-0.1562499995E-1, Q2= 0.1430488765E-3, Q3=-0.6911147651E-5,
    Q4= 0.7621095161E-6, Q5=-0.9349451520E-7,
    R1= 57568490574.0, R2=-13362590354.0, R3=651619640.7,
    R4=-11214424.18, R5= 77392.33017, R6=-184.9052456,
    S1= 57568490411.0, S2=1029532985.0, S3=9494680.718,
    S4= 59272.64853, S5=267.8532712, S6=1.0;
  double
    AX,FR,FS,Z,FP,FQ,XX,Y, TMP;

  if (X==0.0) return 1.0;
  AX = fabs(X);
  if (AX < 8.0) {
    Y = X*X;
    FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
    FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
    TMP = FR/FS;
  }
  else {
    Z = 8./AX;
    Y = Z*Z;
    XX = AX-0.785398164;
    FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
    FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
    TMP = sqrt(0.636619772/AX)*(FP*cos(XX)-Z*FQ*sin(XX));
  }
  return TMP;
}

double BESSJ1 (double X) {
  /**********************************************************************
      This subroutine calculates the First Kind Bessel Function of
      order 1, for any real number X. The polynomial approximation by
      series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
      REFERENCES:
      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
      VOL.5, 1962.
  ***********************************************************************/
  const double  
    P1=1.0, P2=0.183105E-2, P3=-0.3516396496E-4, P4=0.2457520174E-5,
    P5=-0.240337019E-6,  P6=0.636619772,
    Q1= 0.04687499995, Q2=-0.2002690873E-3, Q3=0.8449199096E-5,
    Q4=-0.88228987E-6, Q5= 0.105787412E-6,
    R1= 72362614232.0, R2=-7895059235.0, R3=242396853.1,
    R4=-2972611.439,   R5=15704.48260,  R6=-30.16036606,
    S1=144725228442.0, S2=2300535178.0, S3=18583304.74,
    S4=99447.43394,    S5=376.9991397,  S6=1.0;

  double AX,FR,FS,Y,Z,FP,FQ,XX, TMP;

  AX = fabs(X);
  if (AX < 8.0) {
    Y = X*X;
    FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
    FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
    TMP = X*(FR/FS);
  }
  else {
    Z = 8.0/AX;
    Y = Z*Z;
    XX = AX-2.35619491;
    FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
    FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
    TMP = sqrt(P6/AX)*(cos(XX)*FP-Z*sin(XX)*FQ)*Sign(S6,X);
  }
  return TMP;
}

double BESSJ (int N, double X) {
  /************************************************************************
      This subroutine calculates the first kind modified Bessel function
      of integer order N, for any REAL X. We use here the classical
      recursion formula, when X > N. For X < N, the Miller's algorithm
      is used to avoid overflows.
      ----------------------------- 
      REFERENCE:
      C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
      MATHEMATICAL TABLES, VOL.5, 1962.
  *************************************************************************/
  const int IACC = 40; 
  const double BIGNO = 1e10,  BIGNI = 1e-10;
    
  double TOX,BJM,BJ,BJP,SUM,TMP;
  int J, JSUM, M;

  if (N == 0) return BESSJ0(X);
  if (N == 1) return BESSJ1(X);
  if (X == 0.0) return 0.0;

  TOX = 2.0/X;
  if (X > 1.0*N) {
    BJM = BESSJ0(X);
    BJ  = BESSJ1(X);
    for (J=1; J<N; J++) {
      BJP = J*TOX*BJ-BJM;
      BJM = BJ;
      BJ  = BJP;
    }
    return BJ;
  }
  else {
    M = (int) (2*((N+floor(sqrt(1.0*(IACC*N))))/2));
    TMP = 0.0;
    JSUM = 0;
    SUM = 0.0;
    BJP = 0.0;
    BJ  = 1.0;
    for (J=M; J>0; J--) {
      BJM = J*TOX*BJ-BJP;
      BJP = BJ;
      BJ  = BJM;
      if (fabs(BJ) > BIGNO) {
	BJ  = BJ*BIGNI;
	BJP = BJP*BIGNI;
	TMP = TMP*BIGNI;
	SUM = SUM*BIGNI;
      }
      if (JSUM != 0)  SUM += BJ;
      JSUM = 1-JSUM;
      if (J == N)  TMP = BJP;
    }
    SUM = 2.0*SUM-BJ;
    return (TMP/SUM);
  }
}





/**************************************************************
 *     Purpose: This program computes the spherical Bessel      *
 *              functions jn(x) and jn'(x) using subroutine     *
 *              SPHJ                                            *
 *     Input :  x --- Argument of jn(x)                         *
 *              n --- Order of jn(x)  ( n = 0 to 250 )          * 
 *     Output:  SJ(n) --- jn(x)                                 *
 *              DJ(n) --- jn'(x)                                * 
 *     Example:   x =10.0                                       *
 *                n          jn(x)              jn'(x)          *
 *              --------------------------------------------    *
 *                0    -.5440211109D-01    -.7846694180D-01     * 
 *                1     .7846694180D-01    -.7009549945D-01     *
 *                2     .7794219363D-01     .5508428371D-01     *
 *                3    -.3949584498D-01     .9374053162D-01     *
 *                4    -.1055892851D+00     .1329879757D-01     *
 *                5    -.5553451162D-01    -.7226857814D-01     *
 * ------------------------------------------------------------ *
 * REFERENCE: "Fortran Routines for Computation of Special      *
 *             Functions,                                       *
 *             jin.ece.uiuc.edu/routines/routines.html".        *
 *                                                              *
 *                     Visual C++ Release By J-P Moreau, Paris. *
 *                                 (www.jpmoreau.fr)            *
 ***************************************************************/
#include <stdio.h>
#include <math.h>

static int MSTA1(double, int);
static int MSTA2(double,int,int);
static double ENVJ(int N, double X);



void SPHJ(int N, double X, int *NM, double *SJ, double *DJ) {
  /*     =======================================================
	 !      Purpose: Compute spherical Bessel functions jn(x) and
	 !               their derivatives
	 !      Input :  x --- Argument of jn(x)
	 !               n --- Order of jn(x)  ( n = 0,1,úúú )
	 !      Output:  SJ(n) --- jn(x)
	 !               DJ(n) --- jn'(x)
	 !               NM --- Highest order computed
	 !      Routines called:
	 !               MSTA1 and MSTA2 for computing the starting
	 !               point for backward recurrence
	 !      ======================================================= */
  int K, M; double CS, F, F0, F1, SA, SB;
        
  *NM=N;
  if (fabs(X) < 1e-100) {
    for (K=0; K<=N; K++) {
      SJ[K]=0.0;
      DJ[K]=0.0;
    } 
    SJ[0]=1.0;
    DJ[1]=0.333333333333333;
    return;
  }
  SJ[0]=sin(X)/X;
  SJ[1]=(SJ[0]-cos(X))/X;
  if (N >= 2) {
    SA=SJ[0];
    SB=SJ[1];
    M=MSTA1(X,200);
    if (M < N)
      *NM=M;
    else
      M=MSTA2(X,N,15);
    F0=0.0;
    F1=1.0-100;
    for (K=M; K>-1; K--) {
      F=(2.0*K+3.0)*F1/X-F0;
      if (K <= *NM)  SJ[K]=F;
      F0=F1;
      F1=F;
    }
    if (fabs(SA) > fabs(SB))  CS=SA/F;
    if (fabs(SA) <= fabs(SB)) CS=SB/F0;
    for (K=0; K<=*NM; K++)  SJ[K] *= CS;
  }      
  DJ[0]=(cos(X)-sin(X)/X)/X;
  for (K=1; K<=*NM; K++)
    DJ[K]=SJ[K-1]-(K+1.0)*SJ[K]/X;
}

int MSTA1(double X, int MP) {

  /*     ===================================================
	 !      Purpose: Determine the starting point for backward  
	 !               recurrence such that the magnitude of    
	 !               Jn(x) at that point is about 10^(-MP)
	 !      Input :  x     --- Argument of Jn(x)
	 !               MP    --- Value of magnitude
	 !      Output:  MSTA1 --- Starting point   
	 !      =================================================== */

  double A0,F,F0,F1; int IT,NN,N0,N1;
  A0=fabs(X);
  N0=floor(1.1*A0)+1;
  F0=ENVJ(N0,A0)-MP;
  N1=N0+5;
  F1=ENVJ(N1,A0)-MP;
  for (IT=1; IT<=20; IT++) {             
    NN=N1-(N1-N0)/(1.0-F0/F1);                  
    F=ENVJ(NN,A0)-MP;
    if (abs(NN-N1) < 1) goto e20;
    N0=N1;
    F0=F1;
    N1=NN;
    F1=F;
  }
 e20:	return NN;
}

int MSTA2(double X, int N, int MP) {

  /*     ===================================================
	 !      Purpose: Determine the starting point for backward
	 !               recurrence such that all Jn(x) has MP
	 !               significant digits
	 !      Input :  x  --- Argument of Jn(x)
	 !               n  --- Order of Jn(x)
	 !               MP --- Significant digit
	 !      Output:  MSTA2 --- Starting point
	 !      =================================================== */

  double A0,EJN,F,F0,F1,HMP,OBJ; int IT,N0,N1,NN;
  A0=fabs(X);
  HMP=0.5*MP;
  EJN=ENVJ(N,A0);
  if (EJN <= HMP) {
    OBJ=MP;
    N0=floor(1.1*A0);
  }
  else {
    OBJ=HMP+EJN;
    N0=N;
  }
  F0=ENVJ(N0,A0)-OBJ;
  N1=N0+5;
  F1=ENVJ(N1,A0)-OBJ;
  for (IT=1; IT<=20; IT++) {
    NN=N1-(N1-N0)/(1.0-F0/F1);
    F=ENVJ(NN,A0)-OBJ;
    if (abs(NN-N1) < 1) goto e20;
    N0=N1;
    F0=F1;
    N1=NN;
    F1=F;
  }
 e20:    return NN+10;
}


double ENVJ(int N, double X) {
  return (0.5*log(6.28*N)-N*log(1.36*X/N));
}
