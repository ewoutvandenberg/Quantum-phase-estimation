/*
# Copyright 2021 IBM Inc. All rights reserved
# SPDX-License-Identifier: Apache2.0

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This file is part of the code to reproduce the results in the paper:
# E. van den Berg, "Efficient Bayesian phase estimation using mixed priors"
# arXiv:2007.11629.
*/

#include <sys/stat.h> /* mkdir */
#include <errno.h>    /* EEXIST */
#include <stdlib.h>   /* malloc, free */
#include <string.h>   /* memset, memcpy */
#include <time.h>

#include "bpe.h"


/* Debug settings */
#ifdef VALIDATE_MEAN
FILE *fp_validate_mean = NULL;
#endif



/* ------------------------------------------------------------------------ */
void initDistribution(TDistribution *distribution)
/* ------------------------------------------------------------------------ */
{
   distribution -> gridpoints = NULL;
   distribution -> mass       = NULL;
   distribution -> n          = 0;
   distribution -> capacity   = 0;
}


/* ------------------------------------------------------------------------ */
void freeDistribution(TDistribution *distribution)
/* ------------------------------------------------------------------------ */
{
   if (distribution -> gridpoints) free(distribution -> gridpoints);
   if (distribution -> mass      ) free(distribution -> mass);

   initDistribution(distribution);
}


/* ------------------------------------------------------------------------ */
int resizeDistribution(TDistribution *distribution, int size)
/* ------------------------------------------------------------------------ */
{
   if (distribution == NULL) return -1;

   if (size > distribution -> capacity)
   {  freeDistribution(distribution);

      distribution -> gridpoints = (double *)malloc(sizeof(double) * (size+1));
      distribution -> mass       = (double *)malloc(sizeof(double) * size);
      distribution -> capacity   = size;

      if ((distribution -> gridpoints == NULL) || (distribution -> mass == NULL))
      {  freeDistribution(distribution);
         return -1;
      }
   }

   distribution -> n = size;

   return 0;
}


/* ------------------------------------------------------------------------ */
int linearGrid(TDistribution *distribution, int n)
/* ------------------------------------------------------------------------ */
{  double a, b = 0;
   int i;

   if (n <= 0) return -1;
   if (resizeDistribution(distribution, n) != 0) return -1;


   distribution -> gridpoints[0] = 0;
   for (i = 1; i < n; i++)
   {  distribution -> gridpoints[i] = (i * 2 * M_PI) / n;
   }
   distribution -> gridpoints[n] = 2*M_PI;

   return 0;
}


/* ------------------------------------------------------------------------ */
void uniformPrior(TDistribution *distribution)
/* ------------------------------------------------------------------------ */
{  double a,b;
   int i, n;

   n = distribution -> n; if (n <= 0) return ;

   b = distribution -> gridpoints[0];
   for (i = 0; i < n; )
   {  i ++;
      a = b;
      b = distribution -> gridpoints[i];
      distribution -> mass[i] = (b - a);
   }

   normalizeDistribution(distribution);
}


/* ------------------------------------------------------------------------ */
int uniformPriorAround(TDistribution *distribution, double phi, double delta, double gridspacing)
/* ------------------------------------------------------------------------ */
{  double a, b, theta;
   int idx;

   if (gridspacing <= 0) return -1;
   if (gridspacing > 3) gridspacing = 1.;

   /* Set a uniform prior over a region around phi               */
   /* This assumes we know which of the peaks we should focus on */
   if (resizeDistribution(distribution, ceil(2*M_PI / gridspacing + 10)) != 0)
   {  printf("Error resizing the distribution\n");
      return -1;
   }

   idx = 0;
   distribution -> gridpoints[0] = 0.;

   /* Add grid for interval [a, b) before phi */
   a = 0;
   b = (phi - delta); if (b < 0) b = 0;
   theta = a;
   while (theta < b)
   {  theta += gridspacing; if (theta > b) theta = b;
      distribution -> gridpoints[idx+1] = theta;
      distribution -> mass[idx]         = 0.;
      idx += 1;
   }

   /* Add grid for interval [a, b) around phi*/
   a = b;
   b = phi + delta; if (b > 2*M_PI) b = 2*M_PI;
   theta = a;
   while (theta < b)
   {  theta += gridspacing; if (theta > b) theta = b;
      distribution -> gridpoints[idx+1] = theta;
      distribution -> mass[idx]         = theta - distribution -> gridpoints[idx];
      idx += 1;
   }

   /* Add grid for interval [a, b) after phi */
   a = b;
   b = (2*M_PI + phi - delta); if (b > 2*M_PI) b = 2*M_PI;
   theta = a;
   while (theta < b)
   {  theta += gridspacing; if (theta > b) theta = b;
      distribution -> gridpoints[idx+1] = theta;
      distribution -> mass[idx]         = 0.;
      idx += 1;
   }

   /* Add grid for interval [a, b) for wrap-around of phi*/
   a = b;
   b = 2*M_PI;
   theta = a;
   while (theta < b)
   {  theta += gridspacing; if (theta > b) theta = b;
      distribution -> gridpoints[idx+1] = theta;
      distribution -> mass[idx]         = theta - distribution -> gridpoints[idx];
      idx += 1;
   }

   distribution -> n = idx;
   normalizeDistribution(distribution);

   return 0;
}


/* ------------------------------------------------------------------------ */
void normalizeDistribution(TDistribution *distribution)
/* ------------------------------------------------------------------------ */
{  double mass = 0;
   int i, n;

   n = distribution -> n; if (n <= 0) return;

   for (i = 0; i < n; i++) mass += distribution -> mass[i];
   for (i = 0; i < n; i++) distribution -> mass[i] /= mass;
}


/* ------------------------------------------------------------------------ */
double integrateDistribution(TDistribution *distribution, double x, double y)
/* ------------------------------------------------------------------------ */
{  double a, b, p;
   int i, n;

   n = distribution -> n; if (n <= 0) return 0;

   /* Swap x and y if needed, do not negate the probability */
   if (x > y) { a = x; x = y; y = a; }

   if (x < 0)
   {  x += 2*M_PI;
      y += 2*M_PI;
   }
   else if (x >= 2*M_PI)
   {  x -= 2*M_PI;
      y -= 2*M_PI;
   }

   if (y > 2*M_PI)
   {  return integrateDistribution(distribution, x, 2*M_PI) +
             integrateDistribution(distribution, 0, y - 2*M_PI);
   }

   /* Avoid summing up from 0 to x, just in case the integral value is actually very small. */
   p = 0; b = 0;
   for (i = 0; i < n; i++)
   {  a = b; b = distribution -> gridpoints[i+1];
 
      if (x < b)
      {  if (y < b)
         {  /* x and y lie in the same segment */
            p = distribution->mass[i] * (y - x) / (b - a);
            return p;
         }
         else
         {  p = distribution->mass[i] * (b - x) / (b - a);
            i ++;
            break;
         }
      }
   }

   b = distribution -> gridpoints[i];
   for (; i < n; i++)
   {  a = b; b = distribution -> gridpoints[i+1];

      if (y < b)
      {  p += distribution->mass[i] * (y - a) / (b - a);
         break;
      }
      else
      {  p += distribution->mass[i];
      }
   }

   return p;
}


/* ------------------------------------------------------------------------ */
double getMean(TDistribution *distribution, double *resultMSE)
/* ------------------------------------------------------------------------ */
{  double mean, mse, mseBest;
   double mu, mubar, delta;
   double d0,d1, da,db, dap,dbp, a,b,p,s;
   double polycoef[3], root;
   int    idx, idxOpposite;
   int    mode, i;

   /* Initialize the mean squared error */
   mean = 0;
   mseBest = getErrorSquared(distribution, mean);

   /* We need to ensure that the segments are sufficiently small */

   /* Traverse the grid */
   idx = 0;
   for (idxOpposite = 0; idxOpposite < distribution -> n; idxOpposite++)
   {  if (distribution->gridpoints[idxOpposite+1] > M_PI) break;
   }

   /* Traverse the segments */
   mu = 0; mubar = M_PI; mode = 0;
   while (1)
   {  
      /* Distance to end of current and opposite interval */
      d0 = distribution->gridpoints[idx+1] - mu;
      d1 = distribution->gridpoints[idxOpposite+1] - mubar;
      if (d0 < d1)
      {  delta = d0; mode = 0; }
      else
      {  delta = d1; mode = 1; }

      /* Initialize the polynomial coefficients */
      polycoef[0] = 0.;
      polycoef[1] = 0.;
      polycoef[2] = 0.;

      b = 0;
      for (i = 0; i < distribution -> n; i++)
      {  /* Determine the segment and probability density */
         a = b;
         b = distribution -> gridpoints[i+1];
         p = distribution -> mass[i] / (b - a);

         /* Determine the signed distance*/
         da = a - mu; db = b - mu;
         if (i != idxOpposite)
         {  if (idx < idxOpposite)
            {  if (i > idxOpposite)
               {  da = da - 2*M_PI;
                  db = db - 2*M_PI;
               }
            }
            else if (idx > idxOpposite)
            {  if (i < idxOpposite)
               {  da = 2*M_PI + da;
                  db = 2*M_PI + db;
               }
            }
         }
         else /* (i == idxOpposite) */
         {  if (idx < idxOpposite)
            {  db -= 2*M_PI;
            }
            else
            {  da += 2*M_PI;
            }
         }

         /* The integral of [d - delta]^2 from d=da to d=db */
         /* The distance da should always be be less than db */
         dap = p * da;
         dbp = p * db;
         polycoef[2] += (dbp - dap);
         dap *= da;
         dbp *= db;
         polycoef[1] -= (dbp - dap);
         dap *= da;
         dbp *= db;
         polycoef[0] += (dbp - dap) / 3.;

         if (i == idxOpposite) polycoef[0] += p*(2*M_PI * M_PI * M_PI / 3.);
      }

      /* Verification: output the interval and polynomial */
      #ifdef VALIDATE_MEAN
         fprintf(fp_validate_mean, "S %.8f %.8f %.8e %.8e %.8e\n",
                 mu, delta, polycoef[0], polycoef[1], polycoef[2]);
      #endif

      /* Check for local maxima */
      if (polycoef[2] > 0)
      {  root = -polycoef[1] / (2 * polycoef[2]);
         if ((root >= 0) && (root <= delta))
         {  mse = polycoef[0] + polycoef[1] * root + polycoef[2] * (root*root);
            if (mse < mseBest)
            {  mseBest = mse;
               mean = mu + root;
            }

            /* Verification: output the root */
            #ifdef VALIDATE_MEAN
               fprintf(fp_validate_mean, "R %.8f %.8e\n", mu+root, mse);
            #endif
         }
      }
      
      /* Update the segment */
      if (mode == 0)
      {  idx += 1; if (idx >= distribution -> n) break;
         mu = distribution -> gridpoints[idx];
         mubar = (mu < M_PI) ? mu + M_PI : mu - M_PI;
      }
      else
      {  idxOpposite = (idxOpposite + 1) % (distribution -> n);
         mubar = distribution -> gridpoints[idxOpposite];
         mu = (mubar <= M_PI) ? mubar + M_PI : mubar - M_PI;
      }
   }

   /* Verification: output the mean */
   #ifdef VALIDATE_MEAN
      fprintf(fp_validate_mean, "M %.8f %.8e\n", mean, mseBest);
   #endif

   /* Return the mean-squared error */
   if (resultMSE) *resultMSE = mseBest;

   return mean;
}


/* ------------------------------------------------------------------------ */
double getErrorAbs(TDistribution *distribution, double ref)
/* ------------------------------------------------------------------------ */
{  double rbar;
   double error = 0;
   double a,b,da,db,p,t;
   int opposite;
   int i,j;

   /* Opposite of the mean in range [0,2*pi) */
   /* Note that when ref is pi-epsilon, we can have ref+PI = 2pi.  */
   /* This caused difficulties when evaluating the integral over   */
   /* the interval [0,a), since rbar now lies in the interval, but */
   /* the angle difference to 0 is -(pi-epsilon), while a positive */
   /* value is expected. We therefore evaluate rbar as given below.*/
   rbar = (ref < M_PI) ? (ref+M_PI) : (ref-M_PI);

   /* We need to ensure that the segments are sufficiently small */

   b = 0;
   for (i = 0; i < distribution -> n; i++)
   {  /* Determine the segment and probability density */
      a = b;
      b = distribution -> gridpoints[i+1];
      p = distribution -> mass[i] / (b - a); /* Probability density */

      /* Integral of |d|*p dx */

      /* Compute the signed angle differences */
      da = a - ref; db = b - ref;
      opposite = ((a <= rbar) && (rbar < b));
      if (opposite)
      {  if (ref < M_PI)
              db -= 2*M_PI;
         else da += 2*M_PI;
      }
      else
      {  if (da < -M_PI)
         {  da += 2*M_PI; db += 2*M_PI;
         }
         if (da > M_PI)
         {  da -= 2*M_PI;
            db -= 2*M_PI;
         } 
      }

      if (opposite)
      {  /* Integrate over values [da,pi] and [-pi, dbi] */
         error += p*((M_PI * M_PI) - ((da*da) + (db*db)) / 2.);
      }
      else if ((a <= ref) && (ref < b))
      {  /* Integrate over values [0,da] and [0,db] */
         error += p*((da*da) + (db*db)) / 2.;
      }
      else
      {  /* Integrate over values [da,db] */
         if (da > 0)
              error += p*((db*db) - (da*da)) / 2.;
         else error += p*((da*da) - (db*db)) / 2.;
      }
   }

   return error;
}


/* ------------------------------------------------------------------------ */
double getErrorSquared(TDistribution *distribution, double ref)
/* ------------------------------------------------------------------------ */
{  double rbar;
   double error = 0;
   double a,b,da,db,p;
   int opposite;
   int i,j;

   /* Opposite of the mean in range [0,2*pi) */
   /* Note that when ref is pi-epsilon, we can have ref+PI = 2pi.  */
   /* This caused difficulties when evaluating the integral over   */
   /* the interval [0,a), since rbar now lies in the interval, but */
   /* the angle difference to 0 is -(pi-epsilon), while a positive */
   /* value is expected. We therefore evaluate rbar as given below.*/
   rbar = (ref < M_PI) ? (ref+M_PI) : (ref-M_PI);

   /* We need to ensure that the segments are sufficiently small */

   b = 0;
   for (i = 0; i < distribution -> n; i++)
   {  /* Determine the segment and probability density */
      a = b;
      b = distribution -> gridpoints[i+1];
      p = distribution -> mass[i] / (b - a); /* Probability density */

      /* Compute the signed angle differences */
      da = a - ref; db = b - ref;
      opposite = ((a <= rbar) && (rbar < b));
      if (opposite)
      {  if (ref < M_PI)
              db -= 2*M_PI;
         else da += 2*M_PI;
      }
      else
      {  if (da < -M_PI)
         {  da += 2*M_PI; db += 2*M_PI;
         }
         if (da > M_PI)
         {  da -= 2*M_PI;
            db -= 2*M_PI;
         } 
      }

      /* Integral of d^2*p dx */
      error += p * (db*db*db - da*da*da) / 3.;

      if ((a <= rbar) && (rbar < b))
      {  /* Add correction term */
         error += p * (2 * M_PI * M_PI * M_PI) / 3.;
      }
   }

   return error;
}


/* ------------------------------------------------------------------------ */
void updateNoiseless(TDistribution *distribution,
                     double phi, double s, double theta,
                     int measurement)
/* ------------------------------------------------------------------------ */
{  double a,b, pMeasurement;
   int i;

   b = 0;
   for (i = 0; i < distribution -> n; i++)
   {  
      a = b;
      b = distribution -> gridpoints[i+1];

      /* Evaluate the measurement probability at the center   */
      /* of the box. We omit the division by a factor of two; */
      /* this will be done during normalization.              */
      pMeasurement = (1 + (1-2*measurement)*cos(s * ((a+b)/2.) + theta));
      distribution -> mass[i] *= pMeasurement; /* Multiply by prior */
   }

   /* Normalize the distribution */
   normalizeDistribution(distribution);
}


/* ------------------------------------------------------------------------ */
int refineDistribution(TDistribution *src, TDistribution *dst, double pMerge, double pSplit)
/* ------------------------------------------------------------------------ */
{  int i,j;
   double a,b,c, mass;

   if (resizeDistribution(dst, 2 * src -> n) != 0) return -1;

   dst -> gridpoints[0] = 0; b = 0;
   for (i = 0, j = 0; i < src -> n; i++)
   {  a = b;
      b = src -> gridpoints[i+1];

      if (src -> mass[i] > pSplit)
      {  c    = (a + b) / 2.;
         mass = (src -> mass[i]) / 2.;

         dst -> gridpoints[j+1] = c;
         dst -> gridpoints[j+2] = b;
         dst -> mass[j]         = mass;
         dst -> mass[j+1]       = mass;
         j += 2;
      }
      else if ((src -> mass[i  ] < pMerge) && ((i+1) < src -> n) &&
               (src -> mass[i+1] < pMerge) && ((src -> gridpoints[i+2] - a) < 3))
      {  dst -> gridpoints[j+1] = src -> gridpoints[i+2];
         dst -> mass[j]         = src -> mass[i] + src -> mass[i+1];
         j++; i++;
      }
      else
      {  dst -> gridpoints[j+1] = src -> gridpoints[i+1];
         dst -> mass[j]         = src -> mass[i];
         j++;
      }
   }

   dst -> n = j;

   return 0;
}






#ifdef DEBUG_RANDOM
double __sampleRand(void)
{
   return (double)rand() / RAND_MAX;
}

double (*sampleRand)(void) = __sampleRand;
double *randomBuffer = NULL;
int     randomSize = 0;
int     randomIndex = 0;

double sampleRandomBuffer(void)
{
   if (randomIndex >= randomSize)
   {  printf("Ran out of random numbers!\n");
      exit(-1);
   }

   return randomBuffer[randomIndex++];
}

void initRandomBuffer(void)
{  FILE *fp;
   int i,n;

   fp = fopen("random_0.txt","r");
   if (fp == NULL) { printf("Error opening random number file\n"); exit(1); }

   if (fscanf(fp,"%d",&n) != 1) { printf("Error reading random file\n"); exit(1); }
   randomBuffer = (double *)malloc(sizeof(double) * n);
   if (randomBuffer == NULL) { printf("Error allocating random buffer\n"); exit(1); }

   for (i = 0; i < n; i++)
   {
      if (fscanf(fp,"%lf",&(randomBuffer[i])) != 1) { printf("Error reading random file\n"); exit(1); }
   }

   sampleRand = sampleRandomBuffer;
   randomSize = n;
   randomIndex = 0;

   fclose(fp);   
}


void freeRandomBuffer(void)
{
   if (randomBuffer) free(randomBuffer);
}
#endif


/* ------------------------------------------------------------------------ */
inline double sampleUniform(void)
/* ------------------------------------------------------------------------ */
{
#ifdef DEBUG_RANDOM
   return sampleRand();
#else
   return (double)rand() / RAND_MAX;
#endif
}


/* ------------------------------------------------------------------------ */
void sampleAngles(double *ptr, int n)
/* ------------------------------------------------------------------------ */
{  int i;

   for (i = 0; i < n; i++)
   {
      #ifdef DEBUG_RANDOM
         ptr[i] = (2*M_PI) * sampleRand();
      #else
         ptr[i] = (2*M_PI*rand()) / RAND_MAX;
      #endif
   }
}



/* ------------------------------------------------------------------------ */
double sampleDistribution(TDistribution *distribution)
/* ------------------------------------------------------------------------ */
{  double a,b,p,m;
   int i;

   /* Note: The distribution must be normalized! */

   p = sampleUniform();

   for (i = 0; i < distribution -> n; i++)
   {
      m = distribution -> mass[i];
      if (p < m)
      {  a = distribution -> gridpoints[i];
         b = distribution -> gridpoints[i+1];
         return a + (b-a) * p / m;
      }
      else
      {  p -= m;
      }
   }

   return 0;
}


/* ------------------------------------------------------------------------ */
int measureNoiseless(double phi, double s, double theta)
/* ------------------------------------------------------------------------ */
{  double sample;

   /* Sample the uniform distribution */
   sample = (double)rand() / RAND_MAX;
   return (sample < (1 + cos(s*phi + theta))/2.) ? 0 : 1;
}


/* ------------------------------------------------------------------------ */
double angleDiff(double alpha, double beta)
/* ------------------------------------------------------------------------ */
{  double diff;

   /* Normalize the angles */
   if ((alpha < 0) || (alpha >= 2*M_PI))
      alpha -= floor(alpha / (2*M_PI)) * (2*M_PI);
   if ((beta < 0) || (beta >= 2*M_PI))
      beta  -= floor(beta  / (2*M_PI)) * (2*M_PI);

   /* Compute and normalize the difference */
   diff = (alpha < beta) ? (beta - alpha) : (alpha - beta);
   if (diff > M_PI) diff = 2*M_PI - diff;

   return diff;
}


/* ------------------------------------------------------------------------ */
void outputDistribution(FILE *fp, TDistribution *distribution)
/* ------------------------------------------------------------------------ */
{  int i, n = distribution -> n;

   fprintf(fp, "%d ", n);
   for (i = 0; i <= n; i++) fprintf(fp," %.8f", distribution -> gridpoints[i]);
   for (i = 0; i <  n; i++) fprintf(fp," %.8e", distribution -> mass[i]); 
   fprintf(fp,"\n");   
}


/* ------------------------------------------------------------------------ */
void ensureDirectory(const char *str)
/* ------------------------------------------------------------------------ */
{  
   if (mkdir(str, 0700) == 0) return ;
   if (errno == EEXIST) return ;

   printf("Error creating directory '%s'\n", str);
   exit(-1);
}



/* ------------------------------------------------------------------------ */
void sampleRounds(int nEig, const double *weights, const double *phi,
                  int nRounds, const int *k, const double *beta,
                  double *weightsBuffer, int *measurements)
/* ------------------------------------------------------------------------ */
{  double *w1 = (double *)weights;
   double *w2 = weightsBuffer;
   double pZero, s,v;
   int round,i;

   /* ---------------------------------------------------- */
   /* The weightsBuffer must be of size at least nEig, the */
   /* measurement buffer must be of size at least nRounds. */
   /* ---------------------------------------------------- */

   if (nRounds > 0)
   {  round = 0;
      while (1)
      {  pZero = 0;
         for (i = 0; i < nEig; i++)
         {  pZero += w1[i] * ((1 + cos(k[round] * phi[i] + beta[round])) / 2.);
         }
         measurements[round] = (sampleUniform() > pZero);

         /* Avoid spurious updates of weigtsBuffer */
         if (round+1 >= nRounds) break;

         /* Update weightsBuffer */
         s = 0;
         for (i = 0; i < nEig; i++)
         {  v = ((1 + cos(k[round]*phi[i] + beta[round] - measurements[round]*M_PI))/2.);
            v *= w1[i];
            s += v;
            w2[i] = v;
         }

         /* Normalize weightsBuffer and update pointers*/
         for (i = 0; i < nEig; i++) w2[i] /= s;
         w1 = w2;

         /* Update the round */
         round++;
      }
   }
}


/* ------------------------------------------------------------------------ */
void sampleFourierRounds(int nEig, const double *weights, const double *phi,
                         int nRounds, double *weightsBuffer, int *measurements)
/* ------------------------------------------------------------------------ */
{  double *w1 = (double *)weights;
   double *w2 = weightsBuffer;
   double pZero, s,v;
   double k,beta = 0;
   int round,i;

   /* ---------------------------------------------------- */
   /* The weightsBuffer must be of size at least nEig, the */
   /* measurement buffer must be of size at least nRounds. */
   /* ---------------------------------------------------- */
   for (k=1, i=1; i < nRounds; i++) k *= 2;

   if (nRounds > 0)
   {  round = 0;
      while (1)
      {  pZero = 0;
         for (i = 0; i < nEig; i++)
         {  pZero += w1[i] * ((1 + cos(k * phi[i] + beta)) / 2.);
         }
         measurements[round] = (sampleUniform() > pZero);

         /* Avoid spurious updates of weigtsBuffer */
         if (round+1 >= nRounds) break;

         /* Update beta */
         beta -= measurements[round]*M_PI;

         /* Update weightsBuffer */
         s = 0;
         for (i = 0; i < nEig; i++)
         {  v = ((1 + cos(k*phi[i] + beta))/2.);
            v *= w1[i];
            s += v;
            w2[i] = v;
         }

         /* Normalize weightsBuffer and update pointers*/
         for (i = 0; i < nEig; i++) w2[i] /= s;
         w1 = w2;

         /* Update beta and k*/
         beta /= 2;
         k    /= 2;

         /* Update the round */
         round++;
      }
   }
}


/* ======================================================================== */
/* Fourier representation                                                   */
/* ======================================================================== */

/* ------------------------------------------------------------------------ */
TFourier *newFourier(void)
/* ------------------------------------------------------------------------ */
{  TFourier *result;
   result = (TFourier *)malloc(sizeof(TFourier));
   if (result == NULL)
   {  printf("Error allocating Fourier representation");
      exit(-1);
   }

   /* Initialize and return the result */
   initFourier(result);
   return result;
}


/* ------------------------------------------------------------------------ */
void deleteFourier(TFourier *dist)
/* ------------------------------------------------------------------------ */
{
   if (dist)
   {  freeFourier(dist);
      free(dist);
   }
}


/* ------------------------------------------------------------------------ */
void freeFourier(TFourier *dist)
/* ------------------------------------------------------------------------ */
{
   if (dist)
   {  if (dist -> buffer) free(dist -> buffer);
      initFourier(dist);
   }
}


/* ------------------------------------------------------------------------ */
void initFourier(TFourier *dist)
/* ------------------------------------------------------------------------ */
{
   /* The buffer holds both the cosine and sine coefficients. Moreover, */
   /* it contains space for manipulation of negative index k values.    */
   /* For evaluate only the k values >= 0 are used.                     */
   /*                                                                   */
   /* buffer = [kCapacity, 1, kCapacity, kCapacity, 1, kCapacity]       */
   /*                      ^                            ^               */
   /* coefCos -------------+                            |               */
   /* coefSin ------------------------------------------+               */

   dist -> coefCos    = NULL;
   dist -> coefSin    = NULL;
   dist -> buffer     = NULL;
   dist -> kMax       = -1;
   dist -> kCapacity  = -1;   /* Largest available kMax */
}



/* ------------------------------------------------------------------------ */
void initFourierZero(TFourier *dist)
/* ------------------------------------------------------------------------ */
{
   resizeFourier(dist, 0, 0);

   dist -> coefCos[0] = 0;
   dist -> coefSin[0] = 0;
}


/* ------------------------------------------------------------------------ */
void initFourierOne(TFourier *dist)
/* ------------------------------------------------------------------------ */
{
   resizeFourier(dist, 0, 0);

   dist -> coefCos[0] = 1;
   dist -> coefSin[0] = 0;
}


/* ------------------------------------------------------------------------ */
void initFourierUniform(TFourier *dist)
/* ------------------------------------------------------------------------ */
{
   resizeFourier(dist, 0, 0);

   dist -> coefCos[0] = 1. / (2*M_PI);
   dist -> coefSin[0] = 0;
}


/* ------------------------------------------------------------------------ */
void initFourierCosine(TFourier *dist, int nCoef, int *k, double *beta, double *weights)
/* ------------------------------------------------------------------------ */
{  double w = 1.;
   int i, kMax;

   /* Note: k values are assumed to be nonnegative */

   if (nCoef <= 0)
   {  dist -> kMax = -1;
      return ;
   }

   /* Determine the maximum k value */
   kMax = k[0];
   for (i = 1; i < nCoef; i++)
   {  if (k[i] > kMax) kMax = k[i];
   }

   /* Resize the array and reset coefficients */
   resizeFourier(dist, kMax, 1);

   /* Update the coefficients */
   for (i = 0; i < nCoef; i++)
   {  if (weights) w = weights[i];
      dist->coefCos[k[i]] += w * cos(beta[i]);
      dist->coefSin[k[i]] -= w * sin(beta[i]);
   }

   /* Standardize */
   if (kMax >= 0) dist -> coefSin[0] = 0;
}


/* ------------------------------------------------------------------------ */
void initFourierSine(TFourier *dist, int nCoef, int *k, double *beta, double *weights)
/* ------------------------------------------------------------------------ */
{  double w = 1.;
   int i, kMax;

   /* Note: k values are assumed to be nonnegative */

   if (nCoef <= 0)
   {  dist -> kMax = -1;
      return ;
   }

   /* Determine the maximum k value */
   kMax = k[0];
   for (i = 1; i < nCoef; i++)
   {  if (k[i] > kMax) kMax = k[i];
   }

   /* Resize the array and reset coefficients */
   resizeFourier(dist, kMax, 1);

   /* Update the coefficients */
   for (i = 0; i < nCoef; i++)
   {  if (weights) w = weights[i];
      dist->coefCos[k[i]] += w * sin(beta[i]);
      dist->coefSin[k[i]] += w * cos(beta[i]);
   }

   /* Standardize */
   if (kMax >= 0) dist -> coefSin[0] = 0;
}


/* ------------------------------------------------------------------------ */
void initFourierProb(TFourier *dist, int nCoef, int *k, double *beta, int *measurements)
/* ------------------------------------------------------------------------ */
{  TFourier *scratch, *distSrc, *distDst, *distSwap;
   double    betaHat, cosBeta, sinBeta, v1, v2;
   int i,j, kMax, kVal;

   /* Probability based on no rounds */
   if (nCoef <= 0)
   {  initFourierOne(dist);
      return ;
   }

   /* Determine the sum of k values */
   kMax = 0;
   for (i = 0; i < nCoef; i++) kMax += k[i];

   /* Ensure the capacity of the distribution */
   setFourierCapacity(dist, kMax);

   /* Special case of a single term */
   if (nCoef == 1)
   {  kVal = k[0];

      /* Set the size and clear the coefficients */
      dist -> kMax = kMax;
      clearFourier(dist);

      /* Set the coefficients */
      betaHat = beta[0]; if (measurements!=NULL) betaHat -= measurements[0] * M_PI;
      if (kVal == 0)
      {  dist -> coefCos[0] = (1 + cos(betaHat)) / 2;
      }
      else
      {  dist -> coefCos[0]    = 0.5;
         dist -> coefCos[kVal] = cos(betaHat) /  2;
         dist -> coefSin[kVal] = sin(betaHat) / -2;
      }
      return ;
   }

   /* Initialize the distribution to 1 */
   initFourierOne(dist);

   /* Create a scratch buffer */
   scratch = newFourier();
   setFourierCapacity(scratch, kMax);

   /* Set the distributions */
   distSrc = dist;
   distDst = scratch;

   /* Update the distributions */
   for (i = 0; i < nCoef; i++)
   {
      /* Get the k and beta values */
      kVal = k[i];
      betaHat = beta[i]; if (measurements!=NULL) betaHat -= measurements[i] * M_PI;
      cosBeta = cos(betaHat) / 4;
      sinBeta = sin(betaHat) / 4;

      /* Set the new size */
      resizeFourier(distDst, distSrc -> kMax + kVal, 1);

      /* Copy the initial coefficients and scale by 1/2 */
      for (j = 0; j <= distSrc->kMax; j++)
      {  distDst -> coefCos[j] = distSrc -> coefCos[j] / 2.;
         distDst -> coefSin[j] = distSrc -> coefSin[j] / 2.;
      }

      /* Update the coefficients */
      for (j = 0; j <= distSrc->kMax; j++)
      {  v1 = cosBeta * distSrc -> coefCos[j];
         v2 = sinBeta * distSrc -> coefSin[j];
         distDst -> coefCos[j-kVal] += v1 - v2;
         distDst -> coefCos[j+kVal] += v1 + v2;

         v1 = sinBeta * distSrc -> coefCos[j];
         v2 = cosBeta * distSrc -> coefSin[j];
         distDst -> coefSin[j-kVal] += v1 + v2;
         distDst -> coefSin[j+kVal] += v2 - v1;
      }

      /* Make sure sine coefficients is zero */
      distDst -> coefSin[0] = 0;

      /* Map all negative indices to positive ones */
      foldFourier(distDst);

      /* Exchange the two distributions */
      distSwap = distSrc; distSrc = distDst; distDst = distSwap;
   }

   /* Swap the distribution contents if needed */
   if (nCoef % 2 != 0) swapFourier(dist, scratch);

   /* Free the scratch space */
   deleteFourier(scratch);
}


/* ------------------------------------------------------------------------ */
void initFourierNormal(TFourier *dist, int kMax, double mu, double sigma)
/* ------------------------------------------------------------------------ */
{  double v1,v2;
   int k;

   /* Initialize the coefficients */
   resizeFourier(dist, kMax, 1);

   for (k = 0; k <= kMax; k++)
   {  v1 = sigma * k;
      v1 = exp(-((v1*v1)/2)) / M_PI;
      v2 = mu * k;
      dist -> coefCos[k] = cos(v2) * v1;
      dist -> coefSin[k] = sin(v2) * v1;
   }
   dist -> coefCos[0] /= 2;
}


/* ------------------------------------------------------------------------ */
void setFourierCapacity(TFourier *dist, int kCapacity)
/* ------------------------------------------------------------------------ */
{  double *buffer;
   double *coefCos, *coefSin;
   int     kMax = dist -> kMax;

   if (kCapacity <= dist -> kCapacity) return ;

   buffer = (double *)malloc(sizeof(double) * (4*kCapacity + 2));
   if (buffer == 0)
   {  printf("Error allocating buffer for Fourier representation");
      exit(-1);
   }

   /* Determine the new pointers */
   coefCos = buffer + kCapacity;
   coefSin = buffer + (3*kCapacity+1);

   /* Copy existing buffers */
   if (kMax >= 0)
   {  memcpy(coefCos - kMax, dist->coefCos - kMax, (2*kMax + 1) * sizeof(double));
      memcpy(coefSin - kMax, dist->coefSin - kMax, (2*kMax + 1) * sizeof(double));
   }

   /* Update the pointers and capacity */
   if (dist -> buffer) free(dist -> buffer);
   dist -> buffer    = buffer;
   dist -> coefCos   = coefCos;
   dist -> coefSin   = coefSin;
   dist -> kCapacity = kCapacity;
}


/* ------------------------------------------------------------------------ */
void resizeFourier(TFourier *dist, int kMax, int clearMode)
/* ------------------------------------------------------------------------ */
{  int k = dist -> kMax;

   /* --------------------------------------------- */
   /* clearMode = 0   Do not zero out anything      */
   /*           = 1   Zero out all coefficients     */
   /*           = 2   Zero out added coefficients   */
   /* --------------------------------------------- */

   /* Ensure the capacity is large enough */
   setFourierCapacity(dist, kMax);

   /* Update kMax */
   dist->kMax = kMax;

   /* Zero out added parts */
   if (clearMode == 1)
   {  /* Zero out all coeffients */
      clearFourier(dist);
   }
   else if (clearMode == 2)
   {  /* Zero out added coefficients */
      if (kMax > k)
      {  memset(dist->coefCos - kMax, 0, sizeof(double) * (kMax - k));
         memset(dist->coefCos + k+1,  0, sizeof(double) * (kMax - k));
         memset(dist->coefSin - kMax, 0, sizeof(double) * (kMax - k));
         memset(dist->coefSin + k+1,  0, sizeof(double) * (kMax - k));
      }
   }
}


/* ------------------------------------------------------------------------ */
void clearFourier(TFourier *dist)
/* ------------------------------------------------------------------------ */
{  int k = dist -> kMax;

   if (k >= 0)
   {  memset(dist->coefCos - k, 0, sizeof(double) * (2*k+1));
      memset(dist->coefSin - k, 0, sizeof(double) * (2*k+1));
   }
}


/* ------------------------------------------------------------------------ */
void __swapInt(int *p1, int *p2)
/* ------------------------------------------------------------------------ */
{  int v;
   v = *p1; *p1 = *p2; *p2 = v;
}


/* ------------------------------------------------------------------------ */
void __swapDoublePtr(double **p1, double **p2)
/* ------------------------------------------------------------------------ */
{  double *v;
   v = *p1; *p1 = *p2; *p2 = v;
}


/* ------------------------------------------------------------------------ */
void swapFourier(TFourier *distA, TFourier *distB)
/* ------------------------------------------------------------------------ */
{
   /* Exchange the pointers and dimensions */
   __swapDoublePtr(&(distA -> coefCos), &(distB -> coefCos));
   __swapDoublePtr(&(distA -> coefSin), &(distB -> coefSin));
   __swapDoublePtr(&(distA -> buffer ), &(distB -> buffer ));
   __swapInt(&(distA -> kMax     ), &(distB -> kMax     ));
   __swapInt(&(distA -> kCapacity), &(distB -> kCapacity));
}


/* ------------------------------------------------------------------------ */
void copyFourier(TFourier *distA, TFourier *distB)
/* ------------------------------------------------------------------------ */
{  int k = distB -> kMax;

   /* Update the size of A */
   resizeFourier(distA, k, 0);

   if (k >= 0)
   {  memcpy(distA -> coefCos - k, distB -> coefCos - k, sizeof(double) * (2*k+1));
      memcpy(distA -> coefSin - k, distB -> coefSin - k, sizeof(double) * (2*k+1));
   }
}


/* ------------------------------------------------------------------------ */
void addFourier(TFourier *distA, TFourier *distB)
/* ------------------------------------------------------------------------ */
{  int i, k = distB -> kMax;

   /* Update the size of A -- add zero padding if needed */
   if (k > distA -> kMax)
      resizeFourier(distA, k, 2);

   for (i = 0; i <= k; i++)
   {  distA -> coefCos[i] += distB -> coefCos[i];
      distA -> coefSin[i] += distB -> coefSin[i];
   }
}


/* ------------------------------------------------------------------------ */
void subtractFourier(TFourier *distA, TFourier *distB)
/* ------------------------------------------------------------------------ */
{  int i, k = distB -> kMax;

   /* Update the size of A -- add zero padding if needed */
   if (k > distA -> kMax)
      resizeFourier(distA, k, 2);

   for (i = 0; i <= k; i++)
   {  distA -> coefCos[i] -= distB -> coefCos[i];
      distA -> coefSin[i] -= distB -> coefSin[i];
   }
}


/* ------------------------------------------------------------------------ */
void multiplyFourier(TFourier *distA, TFourier *distB)
/* ------------------------------------------------------------------------ */
{  TFourier *distC;
   double v, *src, *dst;
   int i,j,k;

   /* As an optimization we could check the multiplication order */
   /* of distA and distB and select the one that is fastest.     */

   /* Deal with the case where one of the representations is empty */
   if ((distA -> kMax < 0) || (distB -> kMax < 0))
   {  distA -> kMax = -1;
      return ;
   }

   /* Create a new distribution for the result */
   distC = newFourier();
   resizeFourier(distC, distA -> kMax + distB -> kMax, 1);

   /* Initialize array by dealing with the constant term */
   v = distA -> coefCos[0];
   for (j = 0; j <= distB -> kMax; j++)
   {  distC -> coefCos[j] = distB -> coefCos[j] * v;
      distC -> coefSin[j] = distB -> coefSin[j] * v;
   }

    /* Deal with the remaining terms */
   for (j = 1; j <= distA -> kMax; j++)
   {  v = distA -> coefCos[j] / 2.;
      if (v != 0)
      {  for (k = 0; k <= distB -> kMax; k++)
         {  distC -> coefCos[j+k] += distB -> coefCos[k] * v;
            distC -> coefSin[j+k] += distB -> coefSin[k] * v;
         }
         for (i = 0; i < j; i++)
         {  distC -> coefCos[j-i] += distB -> coefCos[i] * v;
            distC -> coefSin[j-i] -= distB -> coefSin[i] * v;
         }
         for (i = 0; i <= distB -> kMax - j; i++)
         {  distC -> coefCos[i] += distB -> coefCos[i+j] * v;
            distC -> coefSin[i] += distB -> coefSin[i+j] * v;
         }
      }

      v = distA -> coefSin[j] / 2.;
      if (v != 0)
      {  for (k = 0; k <= distB -> kMax; k++)
         {  distC -> coefSin[j+k] += distB -> coefCos[k] * v;
            distC -> coefCos[j+k] -= distB -> coefSin[k] * v;
         }
         for (i = 0; i < j; i++)
         {  distC -> coefSin[j-i] += distB -> coefCos[i] * v;
            distC -> coefCos[j-i] += distB -> coefSin[i] * v;
         }
         for (i = 0; i <= distB -> kMax - j; i++)
         {  distC -> coefSin[i] -= distB -> coefCos[i+j] * v;
            distC -> coefCos[i] += distB -> coefSin[i+j] * v;
         }
      }
   }

   /* Zero out the zeroth sine coefficients */
   distC -> coefSin[0] = 0;

   /* Swap the result */
   swapFourier(distA, distC);

   /* Free the intermediate distribution */
   deleteFourier(distC);   
}


/* ------------------------------------------------------------------------ */
void scaleFourier(TFourier *dist, double scale)
/* ------------------------------------------------------------------------ */
{  int i,k;

   /* Note: only the nonnegative coefficients are normalized */
   k = dist -> kMax;
   for (i = 0; i <= k; i++)
   {  dist -> coefCos[i] *= scale;
      dist -> coefSin[i] *= scale;
   }
}


/* ------------------------------------------------------------------------ */
void normalizeFourier(TFourier *dist)
/* ------------------------------------------------------------------------ */
{
   /* Note: only the nonnegative coefficients are normalized */
   scaleFourier(dist, 1. / (dist -> coefCos[0] * (2*M_PI)));
}


/* ------------------------------------------------------------------------ */
void foldFourier(TFourier *dist)
/* ------------------------------------------------------------------------ */
{  double *ptr;
   int i, j, k;

   /* Initialize k and offset j */
   k = dist -> kMax;
   if (k <= 0) return ;
   j = 2 * k;

   /* Fold the cosine coefficients */
   ptr = dist -> coefCos - k;
   for (i = 0; i < k; i++)
   {  ptr[j-i] += ptr[i];
      ptr[i]    = 0.;
   }

   /* Fold the sine coefficients */
   ptr = dist -> coefSin - k;
   for (i = 0; i < k; i++)
   {  ptr[j-i] -= ptr[i];
      ptr[i]    = 0.;
   }
}


/* ------------------------------------------------------------------------ */
void truncateFourier(TFourier *dist, int kMax)
/* ------------------------------------------------------------------------ */
{
   if (dist -> kMax > kMax) dist -> kMax = (kMax >= 0) ? kMax : -1;
}


/* ------------------------------------------------------------------------ */
double evaluateFourier(TFourier *dist, double phi)
/* ------------------------------------------------------------------------ */
{  double result;
   int k;

   /* Note: only the nonnegative coefficients are used */
   result = 0;
   for (k = 0; k <= dist -> kMax; k++)
   {  result += dist->coefCos[k] * cos(k*phi);
      result += dist->coefSin[k] * sin(k*phi);
   }

   return result;
}


/* ------------------------------------------------------------------------ */
double getFourierMean(TFourier *dist)
/* ------------------------------------------------------------------------ */
{  double angle;

   if (dist -> kMax <= 0) return 0;

   angle = atan2(dist->coefSin[1], dist->coefCos[1]);
   if (angle < 0) angle += 2*M_PI;
   return angle;
}


/* ------------------------------------------------------------------------ */
double getFourierVar(TFourier *dist)
/* ------------------------------------------------------------------------ */
{  double s, v1, v2;

   if (dist -> kMax < 1) return -1;

   v1 = dist -> coefCos[1]; v1 *= v1;
   v2 = dist -> coefSin[1]; v2 *= v2;
   
   s = (1 / ((M_PI*M_PI) * (v1+v2))) - 1.;
   return s;
}


/* ------------------------------------------------------------------------ */
double getFourierStd(TFourier *dist)
/* ------------------------------------------------------------------------ */
{  double s, v1, v2;

   if (dist -> kMax < 1) return -1;

   v1 = dist -> coefCos[1]; v1 *= v1;
   v2 = dist -> coefSin[1]; v2 *= v2;
   
   s = (1 / ((M_PI*M_PI) * (v1+v2))) - 1.;
   return (s >= 0) ? sqrt(s) : -1;
}


/* ------------------------------------------------------------------------ */
void printFourierCoef(TFourier *dist)
/* ------------------------------------------------------------------------ */
{  double *ptr;
   int i,j;

   for (j = 0; j < 2; j++)
   {
      ptr = (j == 0) ? dist -> coefCos : dist -> coefSin;

      for (i = 0; i <= dist -> kMax; i++)
      {  printf("%.12e   ", ptr[i]);
      }
      printf("\n");
   }
}


/* ------------------------------------------------------------------------ */
double fourierNormalError(int kMax, double sigma)
/* ------------------------------------------------------------------------ */
{
   return erfc(kMax*sigma/sqrt(2)) / (sigma * sqrt(2 * M_PI));
}


/* ------------------------------------------------------------------------ */
double fourierCriticalSigma(int kMax, double epsilon)
/* ------------------------------------------------------------------------ */
{  double low = 0.0, high = 1.0;
   double scale = kMax / sqrt(2.);
   double sigma, error;

   /* Scale the error */ 
   epsilon *= sqrt(2*M_PI);


   /* Ensure the interval contains the desired value */
   while (erfc(high*scale) > epsilon * high)
   {  low = high;
      high *= 2;
   }

   /* Bisection search until sigma is sufficiently accurate */   
   while (high-low > 1e-7)
   {  sigma = (low+high) / 2;
      error = erfc(sigma * scale);
      if (error < epsilon * sigma)
           high = sigma;
      else low  = sigma;
   }

   /* Return sigma for which error < epsilon */
   return high;
}




/* ======================================================================== */
/* MAIN ALGORITHM                                                           */
/* ======================================================================== */

/* ------------------------------------------------------------------------ */
void initSolveWeights(TSolveWeights *data, int mMax, int n, double *x)
/* ------------------------------------------------------------------------ */
{
   /* Initialize the data structure */
   data -> C     = (double *)exitMalloc(sizeof(double) * mMax * n);
   data -> v     = (double *)exitMalloc(sizeof(double) * mMax);
   data -> g     = (double *)exitMalloc(sizeof(double) * n);;
   data -> xbar  = (double *)exitMalloc(sizeof(double) * n);;
   data -> buf   = (double *)exitMalloc(sizeof(double) * n);;
   data -> n     = n;
   data -> mMax  = mMax;

   if (x)
   {  data -> x      = x;
      data -> xOwner = 0;
   }
   else
   {  data -> x      = (double *)exitMalloc(sizeof(double) * n);
      data -> xOwner = 1;
   }
}


/* ------------------------------------------------------------------------ */
void freeSolveWeights(TSolveWeights *data)
/* ------------------------------------------------------------------------ */
{
   /* Free all dynamically allocated memory */
   if (data -> C)    { free(data -> C);    data -> C    = NULL; }
   if (data -> v)    { free(data -> v);    data -> v    = NULL; }
   if (data -> g)    { free(data -> g);    data -> g    = NULL; }
   if (data -> xbar) { free(data -> xbar); data -> xbar = NULL; }
   if (data -> buf)  { free(data -> buf);  data -> buf  = NULL; }

   if ((data -> xOwner) && (data -> x)) free(data -> x);
   data -> x = NULL;
}


/* ------------------------------------------------------------------------ */
void solveWeights(TSolveWeights *data, int m, int verbose)
/* ------------------------------------------------------------------------ */
{  double  alpha, beta = 0.001;
   double  f, fbar, df;
   double  v,s;
   double *C = data -> C;
   int     iteration = 0, iterationMax = 1000;
   int     lsIter, lsIterMax = 20;
   int     i,j,idx,n;

   if ((m > data -> mMax)) return ;
   n = data -> n;

   /* Trivial solution */
   if (n == 1) { data -> x[0] = 1; return; }

   /* Zero out the gradient */
   memset(data -> g, 0, sizeof(double) * n);

   /* We compute v and g simultaneously: first we compute  */
   /* the inner product v = <C_i, x>, and then add vector  */
   /* v*C_i to the gradient. Coefficients Cij and elements */
   /* in x should be nonnegative, and therefore v >= 0.    */
   f = 0;                   /* Initialize the objective value */
   idx = 0;                 /* Initialize the index into C    */
   for (i = 0; i < m; i++)
   {  v = 0;
      for (j = 0; j < n; j++) v += C[idx++] * data -> x[j];
      idx -= n;
      for (j = 0; j < n; j++) data -> g[j] -= C[idx++] / v;
      data -> v[i] = v;
      f -= log(v);
   }

   /* Normalize f */
   f /= m;

   /* Main loop */
   while (1)
   {  
      /* Normalize g */
      for (j = 0; j < n; j++) data -> g[j] /= m;

      /* Line search */
      alpha = 1.;
      for (lsIter = 0; lsIter < lsIterMax; lsIter ++)
      {  /* Compute xbar = x - alpha * g */
         for (j = 0; j < n; j++) data -> xbar[j] = data -> x[j] - alpha * data -> g[j];

         /* Compute the projection of xbar onto the unit-norm simplex */
         solveWeightsProject(data -> xbar, data -> buf, n);

         /* Compute the projected difference for the first full step */
         if (lsIter == 0)
         {  s = 0; df = 0;
            for (j = 0; j < n; j++)
            {  v   = data -> xbar[j] - data -> x[j];
               df += data -> g[j] * v;
               s  += v*v;
            }

            if (s < 1e-8) /* sqrt(s) = ||x - xbar||_2 < 1e-4 */
            {  /* Optimality tolerance reached */
               break;
            }
         }
         else
         {  /* Compute <xbar-x, g> */
            df = 0;
            for (j = 0; j < n; j++)
            {  v   = data -> xbar[j] - data -> x[j];
               df += data -> g[j] * v;
            }
         }

         if (df > 0)
         {  /* Non-descent direction */
            break;
         }

         /* Update v = C*xbar and compute fbar */
         fbar = 0; idx = 0;
         for (i = 0; i < m; i++)
         {  v = 0;
            for (j = 0; j < n; j++) v += C[idx++] * data -> xbar[j];
            data -> v[i] = v;
            fbar -= log(v);
         }
         fbar /= m;

         /* Check Armijo line-search condition */
         if (fbar < f + beta * df)
         {  /* Accept the step */
            break;
         }
         else
         {  alpha /= 2;
         }
      }

      /* Check for optimality */
      if (s < 1e-8)
      {  /* Optimal solution */
         if (verbose) printf("+ %d 0 %.12e\n", lsIter, sqrt(s));
         return ;
      }

      /* Check for line search error */
      if ((lsIter == 20) || (df > 0))
      {  /* Line-search error */
         if (verbose) printf("+ %d 1\n", lsIter);
         return ;
      }

      /* Update x */
      memcpy(data -> x, data -> xbar, sizeof(double) * n);

      /* Increment the iteration */
      iteration ++;
      if (iteration >= iterationMax) break;

      /* Update the gradient */
      memset(data -> g, 0, sizeof(double) * n); /* Zero out the gradient */
      idx = 0;                    /* Initialize the index into C */
      for (i = 0; i < m; i++)
      {  v = data -> v[i];
         for (j = 0; j < n; j++) data -> g[j] -= C[idx++] / v;
      }

      /* Set the new objective value */
      f = fbar;
   }

   /* Maximum number of iterations reached */
   if (verbose) printf("+ %d 2\n", lsIter);
}



/* ------------------------------------------------------------------------ */
void solveWeightsBacktrack(TSolveWeights *data, int m, int verbose)
/* ------------------------------------------------------------------------ */
{  double  alpha, beta = 0.001;
   double  f, fbar, df;
   double  v,s;
   double *C = data -> C;
   int     iteration = 0, iterationMax = 1000;
   int     lsIter, lsIterMax = 20;
   int     i,j,idx,n;

   if ((m > data -> mMax)) return ;
   n = data -> n;

   /* Trivial solution */
   if (n == 1) { data -> x[0] = 1; return; }

   /* Zero out the gradient */
   memset(data -> g, 0, sizeof(double) * n);

   /* We compute v and g simultaneously: first we compute  */
   /* the inner product v = <C_i, x>, and then add vector  */
   /* v*C_i to the gradient. Coefficients Cij and elements */
   /* in x should be nonnegative, and therefore v >= 0.    */
   f = 0;                   /* Initialize the objective value */
   idx = 0;                 /* Initialize the index into C    */
   for (i = 0; i < m; i++)
   {  v = 0;
      for (j = 0; j < n; j++) v += C[idx++] * data -> x[j];
      idx -= n;
      for (j = 0; j < n; j++) data -> g[j] -= C[idx++] / v;
      data -> v[i] = v;
      f -= log(v);
   }

   /* Normalize f */
   f /= m;

   /* Main loop */
   while (1)
   {  
      /* Normalize g */
      for (j = 0; j < n; j++) data -> g[j] /= m;

      /* Compute xbar = x - alpha * g */
      for (j = 0; j < n; j++) data -> xbar[j] = data -> x[j] - 1.0 * data -> g[j];

      /* Compute the projection of xbar onto the unit-norm simplex */
      solveWeightsProject(data -> xbar, data -> buf, n);

      /* Line search */
      alpha = 1.;
      for (lsIter = 0; lsIter < lsIterMax; lsIter ++)
      {  
         /* Compute the projected difference for the first full step */
         if (lsIter == 0)
         {  s = 0; df = 0;

            for (j = 0; j < n; j++)
            {  v   = data -> xbar[j] - data -> x[j];
               df += data -> g[j] * v;
               s  += v*v;
               data -> buf[j] = v;
            }
            if (s < 1e-8) /* sqrt(s) = ||x - xbar||_2 < 1e-4 */
            {  /* Optimality tolerance reached */
               /* Note that xbar is x + alpha * data -> buf */
               break;
            }
         }

         /* We have xbar(alpha) = x + alpha * buf and already computed  */
         /* df = <xbar(1)-x, g> in the very first iteration. It follows */
         /* that <xbar(alpha)-x, g> = alpha * df. */
         if (df > 0)
         {  /* Non-descent direction */
            break;
         }

         /* Update v = C*xbar and compute fbar */
         fbar = 0; idx = 0;
         for (i = 0; i < m; i++)
         {  v = 0;
            for (j = 0; j < n; j++)
            {  data -> xbar[j] = data -> x[j] + alpha * data -> buf[j];
               v += C[idx++] * data -> xbar[j];
            }
            data -> v[i] = v;
            fbar -= log(v);
         }
         fbar /= m;

         /* Check Armijo line-search condition */
         if (fbar < f + beta * (alpha * df))
         {  /* Accept the step */
            break;
         }
         else
         {  alpha /= 2;
         }
      }

      /* Check for optimality */
      if (s < 1e-8)
      {  /* Optimal solution */
         if (verbose) printf("+ %d 0 %.12e\n", lsIter, sqrt(s));
         return ;
      }

      /* Check for line search error */
      if ((lsIter == 20) || (df > 0))
      {  /* Line-search error */
         if (verbose) printf("+ %d 1\n", lsIter);
         return ;
      }

      /* Update x */
      memcpy(data -> x, data -> xbar, sizeof(double) * n);

      /* Increment the iteration */
      iteration ++;
      if (iteration >= iterationMax) break;

      /* Update the gradient */
      memset(data -> g, 0, sizeof(double) * n); /* Zero out the gradient */
      idx = 0;                    /* Initialize the index into C */
      for (i = 0; i < m; i++)
      {  v = data -> v[i];
         for (j = 0; j < n; j++) data -> g[j] -= C[idx++] / v;
      }

      /* Set the new objective value */
      f = fbar;
   }

   /* Maximum number of iterations reached */
   if (verbose) printf("+ %d 2\n", lsIter);
}


/* ------------------------------------------------------------------------ */
void solveWeightsProject(double *vec, double *buf, int n)
/* ------------------------------------------------------------------------ */
{  double x,y;
   int i,j;

   /* Copy the elements of vec into buf for sorting */
   memcpy(buf, vec, sizeof(double) * n);

   /* Sort the entries of x in decreasing order using insertion sort */
   for (i = 1; i < n; i++)
   {  j = i;
      x = buf[i];
      while (j > 0)
      {  y = buf[j-1];
         if (y > x) break;
         buf[j] = y;
         j--;
      }
      buf[j] = x;
   }

   /* Search for thresholding parameter y */
   x = buf[0];
   for (i = 1; i < n; i++)
   {  /* Try value y = buf[i] */
      y = buf[i];
      if ((x - i*y) >= 1) break;
      x += y;
   }
   y = (x-1) / i;

   /* Threshold the original values */
   for (i = 0; i < n; i++)
   {  vec[i] = (vec[i] > y) ? vec[i]-y : 0;
   }
}



/* ------------------------------------------------------------------------ */
void BayesianPhaseEstimation(BPEMeasurementScheme *scheme, BPEConfig *config)
/* ------------------------------------------------------------------------ */
{  /* Density distributions */
   BPEDensity     density;
   TFourier     **densityInterim;
   TFourier      *densityCos;
   TFourier      *densitySin;
   TSolveWeights  dataSolveWeights;

   /* Weight updates */
   int    updateEvery;
   int    updateDouble;
   double updateFactor = 1.5;   /* Scaling factor of updateDouble.        */
   int    updateWeights;

   /* Display settings */
   int    lastDisplayed= -1;
   int    displayIndex  = -1;
   double displaySpacing = 0.01; /* On the log scale (100 points per interval) */

   /* Settings */
   double criticalSigma;

  /* Miscellaneous */
   clock_t start, end;
   double *coefCos, *coefSin;
   double *wHat, *C;
   double  BC, scale, scaleA, scaleB;
   double  mu, sigma;
   double  real, imag;
   double  c, s, v, r;
   int     iter;
   int     i, j, k;
   int     trial;
   int     verbose;
   int     nEig;

   /* Extract configuration parameters */
   nEig = config -> nEig;


   void  (*fptrInitDensity)(BPEConfig *config, BPEDensity **);
   void  (*fptrInitWeights)(BPEConfig *config, double *);


   /* Initialize the distributions */
   initDensity(&density, config -> mode, nEig);
   if (config -> pluginInitDensity != NULL)
   {  config -> pluginInitDensity -> fptrInitDensity(config -> pluginInitDensity, &density);
   }
   else
   {  for (i = 0; i < nEig; i++)
      {  /* Set the desired initial distribution */
         density.mu[i]    = (2*M_PI*i)/(nEig+1);
         density.sigma[i] = 3;

         if (config -> mode != BPE_MODE_NORMAL_ONLY)
         {  initFourierNormal(density.fourier[i], 100, density.mu[i], density.sigma[i]);
         }
         else
         {  density.fourier[i] = NULL;
         }
      }
   }

   /* Initialize the internal distributions */
   densityInterim = exitMalloc(nEig * sizeof(TFourier *));
   densityCos     = newFourier();
   densitySin     = newFourier();
   for (i = 0; i < nEig; i++) densityInterim[i] = newFourier();

   /* Initialize the weight solver information */
   initSolveWeights(&dataSolveWeights, config -> nIter, config -> nEig, density.weights);
   wHat = dataSolveWeights.x;
   C    = dataSolveWeights.C;

   /* Copy weights from the distribution */
   if (config -> fptrInitWeights)
   {  config -> fptrInitWeights(config, density.weights);
   }
   else
   {  for (i = 0; i < nEig; i++) wHat[i] = density.weights[i];
   }

   /* Compute the critical sigma value */
   if (config -> mode == BPE_MODE_MIXED)
   {  criticalSigma = fourierCriticalSigma(config->kMax, config->epsilon);
   }
   else
   {  criticalSigma = -1;
   }

   /* Initialize update double */
   updateEvery  = config -> updateDouble - 1;
   updateDouble = config -> updateDouble;

   /* Print start of new experiment */
/*TOOD*/
   printf("w %7d",0);
   for (i = 0; i < nEig; i++) printf(" %.12e", wHat[i]);
   printf("\n");
/*TOOD*/

   /* Record the start time */
   start = clock();


   /* =========================================== */
   /*                  MAIN LOOP                  */
   /* =========================================== */
   for (iter = 0; iter < config->nIter; iter ++)
   {
      /* Generate measurements */
      if (scheme -> fptrSetK) scheme -> fptrSetK(scheme, &density, iter);
      if (scheme -> fptrSetBeta) scheme -> fptrSetBeta(scheme, &density, iter);
      scheme -> fptrMeasure(scheme, &density, iter);

      /* Compute C values */
      for (i = 0; i < nEig; i++)
      {
         /* Product of P(phi_j) * product P_{kr,beta_r}(m | phi_j) */
         copyFourier(densityInterim[i], scheme -> probability);

         if (density.fourier[i])
         {  multiplyFourier(densityInterim[i], density.fourier[i]);
            C[iter*nEig + i] = densityInterim[i] -> coefCos[0] * (2*M_PI);
         }
         else
         {  /* Normal distribution */
            mu      = density.mu[i];
            sigma   = density.sigma[i];
            k       = densityInterim[i] -> kMax;
            coefCos = densityInterim[i] -> coefCos;
            coefSin = densityInterim[i] -> coefSin;

            s = 0;
            for (j = 0; j <= k; j++)
            {  v = j*sigma;
               s+= (coefCos[j] * cos(j*mu) + coefSin[j]*sin(j*mu)) * exp((v*v)/-2);
            }
            C[iter*nEig + i] = s;
         }
      }

      /* Update weights */
      updateWeights = 0;
      if ((iter < updateEvery) || (iter == updateDouble) || (iter == config->nIter-1))
          updateWeights = 1;


      if (updateWeights)
      {  if (iter == updateDouble)
         {  if (updateDouble*updateFactor < updateDouble + 1)
                 updateDouble += 1;
            else updateDouble = (int)(updateDouble *updateFactor);

            /* Always show iterates where the frequency doubles */
            verbose = 1;
         }
         else
         {  /* Determine the verbosity (Originally had +100 instead of *1.02) */
            verbose = ((iter <= 100) || (iter > lastDisplayed*1.02) || (iter == config->nIter-1));
         }

         if (verbose) lastDisplayed = iter;

         /* Solve the subproblem */
         if (config->curvilinear)
              solveWeights         (&dataSolveWeights, iter+1, verbose && config->displaySolve);
         else solveWeightsBacktrack(&dataSolveWeights, iter+1, verbose && config->displaySolve);

         /* Print updated weights */
         if (verbose)
         {  printf("w %7d", iter+1);
            for (i = 0; i < nEig; i++) printf(" %.12e", wHat[i]);
            printf("\n");
         }
      }

      /* Compute the sum of B*C, with B=wHat */
      for (i = 0, BC = 0; i < nEig; i++) BC += wHat[i] * C[iter*nEig+i];

      /* Update the distributions */
      for (i = 0; i < nEig; i++)
      {  /* Compute the sum of products except index i */
         scale = BC - wHat[i] * C[iter*nEig+i];

         if (density.fourier[i])
         {  /* Compute scale * densityFourier[i] + W[i] * densityInterim[i] */
            if (scale > 0) scaleFourier(density.fourier[i], scale);
            scaleFourier(densityInterim[i], wHat[i]);
            if (scale > 0) addFourier(densityInterim[i], density.fourier[i]);
            normalizeFourier(densityInterim[i]);
            truncateFourier(densityInterim[i], config->kMax);

            /* Compute the standard deviation and compare to the threshold    */
            /* In case the variance is negative we return -1, which will also */
            /* trigger the conversion to a normal distribution, albeit with   */
            /* the distribution as it appeared in the previous iteration.     */
            mu = getFourierMean(densityInterim[i]);
            sigma = getFourierStd(densityInterim[i]);

            if ((sigma < 0) || (isnan(mu)))
            {  /* Do nothing */
               if (scale > 0) scaleFourier(density.fourier[i],1/scale);
               mu    = getFourierMean(density.fourier[i]);
               sigma = getFourierStd(density.fourier[i]);
            }
            else if (sigma < criticalSigma)
            {  /* Convert to a normal distribution */
               mu = getFourierMean(densityInterim[i]);
               if (mu < 0) mu += 2*M_PI;

               /* Delete the Fourier distribution */
               deleteFourier(density.fourier[i]);
               density.fourier[i] = NULL;

               /* Initialize the normal distribution */
               density.mu[i] = mu;
               density.sigma[i] = sigma;

               printf("s %7d %d\n", iter+1, i);
            }
            else
            {  /* Update the distribution */
               swapFourier(density.fourier[i], densityInterim[i]);
            }
         }
         else
         {  /* Update the normal distribution */

            /* Sum{ell != j}C[ell]W[ell] + Prob*W[j] */
            mu      = density.mu[i];
            sigma   = density.sigma[i];
            scaleFourier(densityInterim[i], wHat[i]);
            densityInterim[i] -> coefCos[0] += BC - (wHat[i]*C[iter*nEig + i]);
            scaleFourier(densityInterim[i], 1./BC);

            resizeFourier(densityCos, 1, 1);
            resizeFourier(densitySin, 1, 1);

            densityCos -> coefCos[1] = 1;
            densitySin -> coefSin[1] = 1;
            multiplyFourier(densityCos, densityInterim[i]);
            multiplyFourier(densitySin, densityInterim[i]);

            real = 0;
            imag = 0;
            k = densityCos -> kMax;
            for (j = 0; j <= k; j++)
            {  c = cos(j*mu);
               s = sin(j*mu);
               v = sigma*j;
               v = exp((v*v)/-2);

               real += (densityCos -> coefCos[j] * c + densityCos -> coefSin[j] * s) * v;
               imag += (densitySin -> coefCos[j] * c + densitySin -> coefSin[j] * s) * v;
            }


            /* Compute the new mean and variance */
            mu = atan2(imag, real);
            if (mu < 0) mu += 2*M_PI;

            r = (real*real + imag*imag);
            if (r > 1e-6)
                 v = (1. / r) - 1;
            else v = 1e6;

            if (v > 0)
            {  sigma = sqrt(v);
               if (sigma > 20) sigma = 20;

               density.mu[i]    = mu;
               density.sigma[i] = sigma;
            }
         }
      }      

      /* Display the status */
      s = log10(iter+1) / displaySpacing;
      if ((s >= displayIndex) || (iter == config->nIter - 1))
      {  displayIndex += 1;
         printf("* %7d", iter+1);
         for (i = 0; i < nEig; i++)
         {  if (density.fourier[i])
            {  density.mu[i]    = getFourierMean(density.fourier[i]);
               density.sigma[i] = getFourierStd(density.fourier[i]);
            }
            printf(" %.12e %.12e", density.mu[i],density.sigma[i]);
         }
         printf("\n");
      }

   }

   /* Free distributions */
   freeDensity(&density);

   /* Free internal distributions */
   for (i = 0; i < nEig; i++) deleteFourier(densityInterim[i]);
   free(densityInterim);
   deleteFourier(densityCos);
   deleteFourier(densitySin);

   /* Free weight solver information */
   freeSolveWeights(&dataSolveWeights);

   /* Record the end time */
   end = clock();

   printf("! %.8f\n", (float)(end - start) / CLOCKS_PER_SEC);
}



/* ======================================================================== */
/* HELPER DATA STRUCTURES                                                   */
/* ======================================================================== */

/* ------------------------------------------------------------------------ */
void newPhases(BPEPhases **phases, int nEig)
/* ------------------------------------------------------------------------ */
{  BPEPhases *ptr;

   ptr = *phases = exitMalloc(sizeof(BPEPhases));
   ptr -> phi    = exitMalloc(nEig * sizeof(double));
   ptr -> w      = exitMalloc(nEig * sizeof(double));
   ptr -> nEig   = nEig;
}


/* ------------------------------------------------------------------------ */
void deletePhases(BPEPhases *phases)
/* ------------------------------------------------------------------------ */
{
   if (phases)
   {  nullFree(phases -> phi);
      nullFree(phases -> w);
      nullFree(phases);
   }
}


/* ------------------------------------------------------------------------ */
void initDensity(BPEDensity *density, BPEMode mode, int nEig)
/* ------------------------------------------------------------------------ */
{  int i;

   /* Allocate memory for the distributions */
   density -> fourier = exitMalloc(nEig * sizeof(TFourier *));
   density -> mu      = exitMalloc(nEig * sizeof(double));
   density -> sigma   = exitMalloc(nEig * sizeof(double));
   density -> weights = exitMalloc(nEig * sizeof(double));
   
   /* Initialize the distributions */
   for (i = 0; i < nEig; i++)
   {  if (mode != BPE_MODE_NORMAL_ONLY)
      {  density -> fourier[i] = newFourier();
         initFourierNormal(density -> fourier[i], 100, (2*M_PI*i)/(nEig+1),3);
      }
      else
      {  density -> fourier[i] = NULL;
         density -> mu[i]      = (2*M_PI*i)/(nEig+1);
         density -> sigma[i]   = 3;
      }
   }

   /* Initialize the weights */
   for (i = 0; i < nEig; i++) density -> weights[i] = 1. / nEig;

   /* Set the number of eigenvalues */
   density -> nEig = nEig;
}


/* ------------------------------------------------------------------------ */
void freeDensity(BPEDensity *density)
/* ------------------------------------------------------------------------ */
{  int i;

   if (density)
   {  nullFree(density -> mu);
      nullFree(density -> sigma);
      nullFree(density -> weights);

      if (density -> fourier)
      {  for (i = 0; i < density -> nEig; i++) deleteFourier(density -> fourier[i]);
         nullFree(density -> fourier);
      }
   }
}


/* ------------------------------------------------------------------------ */
void newMeasurementScheme(BPEMeasurementScheme **scheme, int capacity)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   /* Allocate memory for the structure */
   *scheme = ptr = exitMalloc(sizeof(BPEMeasurementScheme));

   ptr -> k            = NULL;
   ptr -> nRounds      = 0;
   ptr -> beta         = NULL;   
   ptr -> betaHat      = NULL;
   ptr -> measurements = NULL;
   ptr -> data         = NULL;
   ptr -> probability  = newFourier();;
   ptr -> phases       = NULL;

   ptr -> noise1       = 0;
   ptr -> noise2       = 0;
   ptr -> state        = NULL;

   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = NULL;
   ptr -> fptrFreeData = NULL;

   /* Set capacity */
   resizeMeasurementScheme(ptr, capacity);
}


/* ------------------------------------------------------------------------ */
void deleteMeasurementScheme(BPEMeasurementScheme *scheme)
/* ------------------------------------------------------------------------ */
{
   if (scheme -> probability) deleteFourier(scheme -> probability);
   if (scheme -> fptrFreeData) scheme -> fptrFreeData(scheme -> data);
   resizeMeasurementScheme(scheme, -1);
   deletePhases(scheme -> phases);
   free(scheme);
}


/* ------------------------------------------------------------------------ */
void resizeMeasurementScheme(BPEMeasurementScheme *scheme, int capacity)
/* ------------------------------------------------------------------------ */
{
   nullFree(scheme -> k);
   nullFree(scheme -> beta);
   nullFree(scheme -> betaHat);
   nullFree(scheme -> measurements);
   nullFree(scheme -> state);

   if (capacity > 0)
   {  scheme -> k            = malloc(capacity * sizeof(int));
      scheme -> beta         = malloc(capacity * sizeof(double));
      scheme -> betaHat      = malloc(capacity * sizeof(double));
      scheme -> measurements = malloc(capacity * sizeof(int));
      scheme -> state        = malloc(capacity * sizeof(int));
      scheme -> nRounds      = capacity;
   }
   else
   {  scheme -> nRounds      = 0;
   }
}


/* ------------------------------------------------------------------------ */
void *exitMalloc(size_t size)
/* ------------------------------------------------------------------------ */
{  void *ptr = malloc(size);
   if (ptr == NULL)
   {  printf("Error allocating memory");
      exit(-1);
   }
   return ptr;
}
