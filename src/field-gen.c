#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>


typedef complex double Complex;
static void get_constants(int model_num, Complex C[2][2][2][4]);


double m2_force_free_magnetic_field(double x[4], double n[4], int model)
{
  Complex C[2][2][2][4];
  get_constants(model, C);
  Complex B[4] = {0, 0, 0, 0};
  double F = 1.0;
  double a = F * sqrt(3.0);
  int i, j, k;
  for (i=0; i<=1; i+=1) {
    for (j=0; j<=1; j+=1) {
      for (k=0; k<=1; k+=1) {
        Complex K[4] = {0, I*(i?F:-F), I*(j?F:-F), I*(k?F:-F)};
        Complex Ikx = K[1]*x[1] + K[2]*x[2] + K[3]*x[3];

        Complex D[4] = {0,
          C[i][j][k][1] + conj(C[!i][!j][!k][1]),
          C[i][j][k][2] + conj(C[!i][!j][!k][2]),
          C[i][j][k][3] + conj(C[!i][!j][!k][3])};

        Complex Q[4] = {0, /* Q = K cross D */
          K[2]*D[3] - K[3]*D[2],
          K[3]*D[1] - K[1]*D[3],
          K[1]*D[2] - K[2]*D[1]};

        Complex P[4] = {0, /* P = K cross Q */
          K[2]*Q[3] - K[3]*Q[2],
          K[3]*Q[1] - K[1]*Q[3],
          K[1]*Q[2] - K[2]*Q[1]};

        Complex H[4] = {0, /* B_k */
          (P[1]/a + Q[1]),
          (P[2]/a + Q[2]),
          (P[3]/a + Q[3])};

        B[1] += H[1] * cexp(Ikx);
        B[2] += H[2] * cexp(Ikx);
        B[3] += H[3] * cexp(Ikx);
      }
    }
  }
  return creal(B[1]*n[1] + B[2]*n[2] + B[3]*n[3]);
}


double m2_force_free_vector_potential(double x[4], double n[4], int model)
{
  Complex C[2][2][2][4];
  get_constants(model, C);
  Complex A[4] = {0, 0, 0, 0};
  double F = 1.0;
  double a = F * sqrt(3.0);
  int i, j, k;
  for (i=0; i<=1; i+=1) {
    for (j=0; j<=1; j+=1) {
      for (k=0; k<=1; k+=1) {
        Complex K[4] = {0, I*(i?F:-F), I*(j?F:-F), I*(k?F:-F)};
        Complex Ikx = K[1]*x[1] + K[2]*x[2] + K[3]*x[3];

        Complex D[4] = {0,
          C[i][j][k][1] + conj(C[!i][!j][!k][1]),
          C[i][j][k][2] + conj(C[!i][!j][!k][2]),
          C[i][j][k][3] + conj(C[!i][!j][!k][3])};

        Complex Q[4] = {0, /* Q = K cross D */
          K[2]*D[3] - K[3]*D[2],
          K[3]*D[1] - K[1]*D[3],
          K[1]*D[2] - K[2]*D[1]};

        Complex H[4] = {0, /* A_k */
          (Q[1]/a + D[1]),
          (Q[2]/a + D[2]),
          (Q[3]/a + D[3])};

        A[1] += H[1] * cexp(Ikx);
        A[2] += H[2] * cexp(Ikx);
        A[3] += H[3] * cexp(Ikx);
      }
    }
  }
  return creal(A[1]*n[1] + A[2]*n[2] + A[3]*n[3]);
}


void get_constants(int num, Complex C[2][2][2][4])
{
  static Complex constants_0[2][2][2][4] =
  {{{
      {0.0,
        +1.2710228918e-01 +5.7420168075e-01*I,
        -4.6624765524e-02 +2.0973467515e-01*I,
        +4.5569268171e-01 -6.3270776431e-01*I},
      {0,
        +5.3536443444e-01 +4.9804111967e-01*I,
        +5.5727913986e-01 -7.8426510608e-02*I,
        +2.6786374266e-01 +2.7726922572e-01*I},
    },
  {
    {0.0,
      -1.8452668691e-01 -6.6961182649e-01*I,
      +2.9657058428e-01 -1.2798393230e-01*I,
      +4.9648869719e-01 -4.0833182440e-01*I},
    {0.0,
      +3.5188694559e-01 -4.7455968488e-01*I,
      -5.3421858192e-01 +5.9338367744e-01*I,
      +6.0667618625e-02 +9.8966868931e-02*I},
  }},
  {{
     {0.0,
       +2.0437027241e-01 -1.6246977297e-01*I,
       -4.5273940537e-01 -3.2542376123e-01*I,
       +6.4582072058e-01 -4.5152892260e-01*I},
     {0.0,
       +4.1197124153e-01 +3.9049710575e-01*I,
       -3.8520090673e-01 -4.9823819948e-01*I,
       -1.3936598699e-01 -5.1161292621e-01*I},
   },
  {
    {0.0,
      +4.1746663956e-01 +4.1609708675e-03*I,
      -5.9875896183e-01 -3.6359721550e-01*I,
      +1.4494217055e-01 -5.6033992262e-01*I},
    {0.0,
      -4.0304620616e-01 -3.8160138898e-01*I,
      -6.4120887911e-02 +6.3660922060e-01*I,
      -2.3610096193e-01 -4.7624330272e-01*I},
  }}};

  static Complex constants_1[2][2][2][4] =
  {{{
      {0.0,
        -2.0752133480e-01 -5.0019504770e-01*I,
        -5.9616114241e-01 -2.8203014301e-01*I,
        -4.5326076673e-01 +2.5757596534e-01*I},
      {0,
        +3.7206686174e-01 +6.6239803560e-01*I,
        +2.1813513904e-01 +1.2042397368e-01*I,
        +2.7025420620e-01 -5.3635145617e-01*I},
    },
  {
    {0.0,
      +4.0419288883e-01 +1.9901435209e-01*I,
      -1.8343814118e-01 +2.6316934297e-01*I,
      +6.4581201302e-01 -5.2634644998e-01*I},
    {0.0,
      +4.0364180864e-02 -3.3114435900e-01*I,
      -1.4281123899e-01 -3.6045088102e-01*I,
      -8.5754287680e-01 +5.4904218787e-02*I},
  }},
  {{
     {0.0,
       +3.4530092572e-01 +4.5542156483e-01*I,
       +2.1093275689e-01 -5.8631971489e-01*I,
       -2.3822458881e-03 -5.3393759738e-01*I},
     {0.0,
       +2.0652792801e-01 +2.1239033484e-01*I,
       +3.3586042155e-01 +3.7218697644e-01*I,
       -5.4549420997e-01 +6.0278292885e-01*I},
   },
  {
    {0.0,
      +3.0185535927e-02 +6.9364816142e-01*I,
      +4.2426174166e-01 +5.5014926606e-01*I,
      +1.8478846934e-01 +3.3645843531e-02*I},
    {0.0,
      +2.4335158764e-01 +7.1136310324e-02*I,
      -4.4885395584e-01 -4.1897345893e-01*I,
      +1.2576597962e-01 -7.3681335191e-01*I},
  }}};

  static Complex constants_2[2][2][2][4] =
  {{{
      {0.0,
        -6.5103812205e-01 +2.8174385625e-01*I,
        +5.4362307961e-01 -2.3820472097e-01*I,
        +8.3242933112e-02 -3.7090812261e-01*I},
      {0,
        -3.6461976053e-01 +1.2382952809e-01*I,
        +2.6884537894e-01 +8.0987062126e-01*I,
        -3.0195054528e-01 +1.7993411421e-01*I},
    },
  {
    {0.0,
      +8.6808728274e-02 -3.3887519154e-01*I,
      +8.1402619849e-02 -6.9709167880e-01*I,
      +6.2046123737e-01 -9.6180565258e-03*I},
    {0.0,
      +5.9750038732e-01 +3.8142192754e-01*I,
      +4.5965959960e-02 +4.7535187780e-01*I,
      +1.6567836816e-01 +4.9192377607e-01*I},
  }},
  {{
     {0.0,
       +4.6531534844e-01 +2.5896232150e-01*I,
       -4.8385703971e-01 -1.5749920773e-01*I,
       +5.3725875557e-01 +4.1091305261e-01*I},
     {0.0,
       -5.5031986862e-03 -1.1571833791e-01*I,
       -6.8163481456e-01 +2.6428649724e-01*I,
       -6.5954763039e-01 +1.3077664651e-01*I},
   },
  {
    {0.0,
      +4.2803365212e-01 -1.7746653654e-01*I,
      -5.7426420519e-01 +4.3413799460e-01*I,
      +4.2785881786e-01 +2.8978350079e-01*I},
    {0.0,
      -4.8518977636e-01 +2.6277894179e-01*I,
      -2.2901016862e-01 -4.5405564290e-01*I,
      +5.0400325231e-01 -4.2767586576e-01*I},
  }}};

  static Complex constants_3[2][2][2][4] =
  {{{
      {0.0,
        +1.7633667487e-02 +4.3243571618e-01*I,
        -3.8915534043e-01 +3.3354618293e-01*I,
        +6.3965669514e-01 +3.7527694090e-01*I},
      {0,
        -2.5699498018e-01 +1.6710617802e-02*I,
        +6.2463892985e-01 -6.4230834524e-01*I,
        -1.9306698313e-01 -3.0604848047e-01*I},
    },
  {
    {0.0,
      +4.2827987379e-01 -4.9708202108e-01*I,
      -3.8187059254e-01 +1.3518194070e-01*I,
      +4.6217387817e-01 -4.3792900546e-01*I},
    {0.0,
      +1.3055128720e-01 +3.3098439098e-01*I,
      -5.1506193177e-01 -4.6572885494e-01*I,
      -3.4944851776e-02 +6.2449370836e-01*I},
  }},
  {{
     {0.0,
       -5.5581720963e-01 +2.6581132295e-01*I,
       -4.1699524152e-01 -5.6877882792e-01*I,
       +3.2458079326e-01 +1.3290783339e-01*I},
     {0.0,
       +5.4327204639e-01 -4.1701968941e-01*I,
       -4.7805912201e-01 +3.2458998116e-01*I,
       -6.1486737027e-02 -4.3962513965e-01*I},
   },
  {
    {0.0,
      -4.9576622936e-01 +3.5215941018e-01*I,
      +3.9566324372e-01 +2.1983167515e-01*I,
      -3.7300035385e-01 +5.3497192809e-01*I},
    {0.0,
      -3.2571587684e-01 +8.2011584324e-01*I,
      +2.4587723387e-01 +3.2986920806e-01*I,
      -1.2504170168e-01 +1.9082566789e-01*I},
  }}};

  int i, j, k, q;
  for (i=0; i<=1; i+=1) {
    for (j=0; j<=1; j+=1) {
      for (k=0; k<=1; k+=1) {
        for (q=0; q<=3; q+=1) {
          switch (num) {
            case 0: C[i][j][k][q] = constants_0[i][j][k][q]; break;
            case 1: C[i][j][k][q] = constants_1[i][j][k][q]; break;
            case 2: C[i][j][k][q] = constants_2[i][j][k][q]; break;
            case 3: C[i][j][k][q] = constants_3[i][j][k][q]; break;
            case 4: C[i][j][k][q] = 0.0; break;
          } 
        }
      }
    }
  }
}
