#include "bpe.h"
#include "bpe_plugins.h"
#include <stdlib.h>  /* malloc */


/* ======================================================================== */
/* Note: The functions in this file are as developed during the various     */
/*       simulations. They could be cleaned up substantially by organizing  */
/*       them as classes with parameters. However, to avoid introducing     */
/*       discrepances, they are left as is.                                 */
/* ======================================================================== */


/* ======================================================================== */
/* GROUND TRUTH FUNCTIONS                                                   */
/* ======================================================================== */

/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruth1(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  BPEPhases *ptr;

   newPhases(phases, 1); ptr = *phases;
   ptr -> phi[0] = 5;  ptr -> w[0] = 1.0;
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruth2(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  BPEPhases *ptr;

   newPhases(phases, 2); ptr = *phases;
   ptr -> phi[0] = 2;  ptr -> w[0] = 0.7;
   ptr -> phi[1] = 4;  ptr -> w[1] = 0.3;
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruth3(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  BPEPhases *ptr;

   newPhases(phases, 3); ptr = *phases;
   ptr -> phi[0] = 2;  ptr -> w[0] = 0.5;
   ptr -> phi[1] = 4;  ptr -> w[1] = 0.3;
   ptr -> phi[2] = 5;  ptr -> w[2] = 0.2;
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruth4(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  BPEPhases *ptr;

   newPhases(phases, 1); ptr = *phases;
   ptr -> phi[0] = 2*M_PI * (1/2. + 1/8. + 1/16.);
   ptr -> w[0]   = 1.0;
}


/* ------------------------------------------------------------------------ */
void _pluginCreateGroundTruthMultiple(BPEPhases **phases, int idx, int n)
/* ------------------------------------------------------------------------ */
{  BPEPhases *ptr;
   double totalWeight;
   int i;

   srand(idx);
   newPhases(phases, n); ptr = *phases;
   for (i = 0; i < n; i++)
   {  /* We allow up to 12 phases; for each value of the  */
      /* number of phases, n, we use the same distance    */
      /* between successive phases. We offset the first   */
      /* phase by 2*pi/24 to avoid wrapping around 0.     */
      /* The final phases are obtained by adding a random */
      /* perturbation uniformly in the range [-0.05,0.05].*/
      ptr -> phi[i]  = (2*M_PI / 12.) * i;
      ptr -> phi[i] += (2*M_PI / 24.);
      ptr -> phi[i] += (0.1 * rand()) / RAND_MAX - 0.05;
      ptr -> w[i] = 0.5 + (0.5 * rand()) / RAND_MAX; /* Between 0.5 and 1 */
   }


   totalWeight = 0;
   for (i = 0; i < n; i++) totalWeight += ptr -> w[i];
   for (i = 0; i < n; i++) ptr -> w[i] /= totalWeight;
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthMultiple1(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  _pluginCreateGroundTruthMultiple(phases, idx, 1);
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthMultiple2(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  _pluginCreateGroundTruthMultiple(phases, idx, 2);
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthMultiple3(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  _pluginCreateGroundTruthMultiple(phases, idx, 3);
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthMultiple4(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  _pluginCreateGroundTruthMultiple(phases, idx, 4);
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthMultiple5(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  _pluginCreateGroundTruthMultiple(phases, idx, 5);
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthMultiple6(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  _pluginCreateGroundTruthMultiple(phases, idx, 6);
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthMultiple7(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  _pluginCreateGroundTruthMultiple(phases, idx, 7);
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthMultiple8(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  _pluginCreateGroundTruthMultiple(phases, idx, 8);
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthMultiple9(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  _pluginCreateGroundTruthMultiple(phases, idx, 9);
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthMultiple10(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  _pluginCreateGroundTruthMultiple(phases, idx, 10);
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthMultiple11(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  _pluginCreateGroundTruthMultiple(phases, idx, 11);
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthMultiple12(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  _pluginCreateGroundTruthMultiple(phases, idx, 12);
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthFixed7(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  BPEPhases *ptr;

   newPhases(phases, 7); ptr = *phases;
   ptr -> phi[0] = 0.26216515;
   ptr -> phi[1] = 0.83374334;
   ptr -> phi[2] = 1.34389931;
   ptr -> phi[3] = 1.78313676;
   ptr -> phi[4] = 2.30625159;
   ptr -> phi[5] = 2.91174359;
   ptr -> phi[6] = 3.44401116;

   ptr -> w[0] = 0.12101411;
   ptr -> w[1] = 0.12789448;
   ptr -> w[2] = 0.14541791;
   ptr -> w[3] = 0.17970617;
   ptr -> w[4] = 0.12482831;
   ptr -> w[5] = 0.18947023;
   ptr -> w[6] = 0.11166879;
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthFixed10(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  BPEPhases *ptr;

   newPhases(phases, 10); ptr = *phases;

   ptr -> phi[0] = 0.30446753;
   ptr -> phi[1] = 0.76945638;
   ptr -> phi[2] = 1.31161402;
   ptr -> phi[3] = 1.85355769;
   ptr -> phi[4] = 2.32111714;
   ptr -> phi[5] = 2.85182662;
   ptr -> phi[6] = 3.41998333;
   ptr -> phi[7] = 3.96716253;
   ptr -> phi[8] = 4.47011048;
   ptr -> phi[9] = 4.96090710;

   ptr -> w[0] = 0.10535719;
   ptr -> w[1] = 0.06587977;
   ptr -> w[2] = 0.10960839;
   ptr -> w[3] = 0.10362178;
   ptr -> w[4] = 0.10033285;
   ptr -> w[5] = 0.12502626;
   ptr -> w[6] = 0.09406571;
   ptr -> w[7] = 0.11082952;
   ptr -> w[8] = 0.09238624;
   ptr -> w[9] = 0.09289229;
}


/* ------------------------------------------------------------------------ */
void pluginCreateGroundTruthFixed12(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  BPEPhases *ptr;

   newPhases(phases, 12); ptr = *phases;
   ptr -> phi[ 0] = 0.25500390;
   ptr -> phi[ 1] = 0.74516641;
   ptr -> phi[ 2] = 1.31492318;
   ptr -> phi[ 3] = 1.84691237;
   ptr -> phi[ 4] = 2.39296300;
   ptr -> phi[ 5] = 2.84754759;
   ptr -> phi[ 6] = 3.36523847;
   ptr -> phi[ 7] = 3.96883633;
   ptr -> phi[ 8] = 4.44257551;
   ptr -> phi[ 9] = 4.93489202;
   ptr -> phi[10] = 5.51758690;
   ptr -> phi[11] = 6.04448261;

   ptr -> w[ 0] = 0.10003630;
   ptr -> w[ 1] = 0.08168390;
   ptr -> w[ 2] = 0.11044348;
   ptr -> w[ 3] = 0.07513798;
   ptr -> w[ 4] = 0.08317102;
   ptr -> w[ 5] = 0.06100895;
   ptr -> w[ 6] = 0.11563351;
   ptr -> w[ 7] = 0.06178559;
   ptr -> w[ 8] = 0.06037355;
   ptr -> w[ 9] = 0.10539436;
   ptr -> w[10] = 0.06071698;
   ptr -> w[11] = 0.08461438;
}


/* ------------------------------------------------------------------------ */
void pluginCreateNoisyGroundTruth3_(BPEPhases **phases, int idx, double noise, int n)
/* ------------------------------------------------------------------------ */
{  BPEPhases *ptr;
   double weightTotal=0;
   int i;

   srand(idx);

   newPhases(phases, 3+n); ptr = *phases;
   ptr -> phi[0] = 2;  ptr -> w[0] = 0.5*(1-noise);
   ptr -> phi[1] = 4;  ptr -> w[1] = 0.3*(1-noise);
   ptr -> phi[2] = 5;  ptr -> w[2] = 0.2*(1-noise);

   for (i = 0; i < n; i++)
   {  ptr -> phi[3+i] = ((2*M_PI)*rand()) / RAND_MAX;
      ptr -> w[3+i]   = (1. * rand())  / RAND_MAX;
      weightTotal += ptr -> w[3+i];
   }

   for (i = 0; i < n; i++)
   {  ptr -> w[3+i] = noise * (ptr -> w[3+i] / weightTotal);
   }
}


/* ------------------------------------------------------------------------ */
void pluginCreateNoisyGroundTruth3_10_2(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  return pluginCreateNoisyGroundTruth3_(phases, idx, 0.10, 2);
}


/* ------------------------------------------------------------------------ */
void pluginCreateNoisyGroundTruth3_10_5(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  return pluginCreateNoisyGroundTruth3_(phases, idx, 0.10, 5);
}


/* ------------------------------------------------------------------------ */
void pluginCreateNoisyGroundTruth3_10_10(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  return pluginCreateNoisyGroundTruth3_(phases, idx, 0.10, 10);
}


/* ------------------------------------------------------------------------ */
void pluginCreateNoisyGroundTruth3_10_20(BPEPhases **phases, int idx)
/* ------------------------------------------------------------------------ */
{  return pluginCreateNoisyGroundTruth3_(phases, idx, 0.10, 20);
}



/* ======================================================================== */
/* MEASUREMENT FUNCTIONS                                                    */
/* ======================================================================== */

/* ------------------------------------------------------------------------ */
void pluginMeasureRounds(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{
   /* Sample the measurements */
   sampleRounds(scheme -> phases -> nEig,
                scheme -> phases -> w,
                scheme -> phases -> phi,
                scheme -> nRounds,
                scheme -> k,
                scheme -> beta,
                scheme -> betaHat,
                scheme -> measurements);

   /* Prepare the probability distribution */
   initFourierProb(scheme -> probability,
                   scheme -> nRounds,
                   scheme -> k,
                   scheme -> beta,
                   scheme -> measurements);
}


/* ------------------------------------------------------------------------ */
void pluginMeasureDecoherence(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  TFourier dist;
   double p;
   int measurement;

   /* Sanity check */
   if (scheme -> nRounds != 1)
   {  printf("Configuration error: plugingMeasureDecoherence only applies to single-round experiments\n");
      exit(-1);
   }

   /* Determine the probability to obtain a correct measurement */
   p = exp((-1. * scheme->k[0]) / (scheme ->kErr));

   /* Sample the measurements */
   sampleRounds(scheme -> phases -> nEig,
                scheme -> phases -> w,
                scheme -> phases -> phi,
                scheme -> nRounds,
                scheme -> k,
                scheme -> beta,
                scheme -> betaHat,
                scheme -> measurements);

   /* Measure the noiseless value with probability (1 - (scheme1->noise1)) */
   /* and measure zero or one with probability (scheme->noise1 / 2.).      */
   measurement = scheme -> measurements[0];
   if (sampleUniform() < (1-p)/2.)
   {  measurement = 1 - measurement;
      scheme -> measurements[0] = measurement;
   }

   /* Initialize the temporary Fourier representation */

   /* Prepare the probability distribution for the original measurement */
   initFourierProb(scheme -> probability, 1, scheme -> k, scheme -> beta, &measurement);
   scaleFourier(scheme -> probability, (1.+p)/2.); /* p + (1-p)/2 */

   /* Prepare the probability distribution for the flipped measurement */
   measurement = 1 - measurement;
   initFourier(&dist);
   initFourierProb(&dist, 1, scheme -> k, scheme -> beta, &measurement);
   scaleFourier(&dist, (1-p)/2.);
   addFourier(scheme -> probability, &dist);
   freeFourier(&dist);
}


/* ------------------------------------------------------------------------ */
void pluginMeasureFourier(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  double p,beta = 0;
   int i,n;

   /* Sample the measurements - without noise */
   sampleFourierRounds(scheme -> phases -> nEig,
                       scheme -> phases -> w,
                       scheme -> phases -> phi,
                       scheme -> nRounds,
                       scheme -> betaHat,
                       scheme -> measurements);

   if ((scheme -> noise1 == 0) && (scheme -> noise2 == 0))
   {  /* Set the beta values */
      beta = 0;
      for (i = 0; i < scheme -> nRounds; i++)
      {  scheme -> beta[i] = beta;
         beta -= scheme -> measurements[i] * M_PI;
         beta /= 2;
      }

      /* Compute the probability */
      initFourierProb(scheme -> probability,
                      scheme -> nRounds,
                      scheme -> k,
                      scheme -> beta,
                      scheme -> measurements);
   }
   else
   {  TFourier dist;
      double prob[4];

      initFourier(&dist);
      prob[0] = 1 - scheme -> noise1; /* 0 -> 0 */
      prob[1] =     scheme -> noise1; /* 0 -> 1 */
      prob[2] =     scheme -> noise2; /* 1 -> 0 */
      prob[3] = 1 - scheme -> noise2; /* 1 -> 1 */

      /* Add noise to the measurements */
      n = scheme -> nRounds;
      for (i = 0; i < n; i++)
      {  if (scheme -> measurements[i] == 0)
         {  if (sampleUniform() < scheme -> noise1)
               scheme -> measurements[i] = 1;
         }
         else
         {  if (sampleUniform() < scheme -> noise2)
               scheme -> measurements[i] = 0;
         }
      }

      /* Set the intial Fourier probability to zero */
      initFourierZero(scheme -> probability);

      /* Try all possible states */
      for (i = 0; i < n; i++)
      {  scheme -> beta[i]  = 0;
         scheme -> state[i] = 0;
      }

      while (1)
      {  /* Compute the probability */
         initFourierProb(&dist,
                      scheme -> nRounds,
                      scheme -> k,
                      scheme -> beta,
                      scheme -> state);

         /* Compute the probability of the measurement given the state */
         p = 1;
         for (i = 0; i < n; i++)
         {  p *= prob[2 * scheme -> state[i] + scheme -> measurements[i]];
         }
         scaleFourier(&dist, p);

         /* Add the contribution */
         addFourier(scheme -> probability, &dist);

         /* Move to the next state */
         i = 0;
         while (i < n)
         {  if ((++(scheme -> state[i])) == 1)
            {  break;
            }

            scheme -> state[i] = 0; i ++;
         }
         if (i == n) break;

         /* Compute beta values */
         beta = 0;
         for (i = 0; i < n; i++)
         {  scheme -> beta[i] = beta;
            beta -= scheme -> state[i] * M_PI;
            beta /= 2;
         }
      }

      freeFourier(&dist);
   }
}


/* ------------------------------------------------------------------------ */
void pluginMeasureNoisy(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  double p;
   int i,n;

   /* Sample the measurements - without noise */
   sampleRounds(scheme -> phases -> nEig,
                scheme -> phases -> w,
                scheme -> phases -> phi,
                scheme -> nRounds,
                scheme -> k,
                scheme -> beta,
                scheme -> betaHat,
                scheme -> measurements);

   if ((scheme -> noise1 == 0) && (scheme -> noise2 == 0))
   {  
      /* Compute the probability */
      initFourierProb(scheme -> probability,
                      scheme -> nRounds,
                      scheme -> k,
                      scheme -> beta,
                      scheme -> measurements);
   }
   else
   {  TFourier dist;
      double prob[4];

      initFourier(&dist);
      prob[0] = 1 - scheme -> noise1; /* 0 -> 0 */
      prob[1] =     scheme -> noise1; /* 0 -> 1 */
      prob[2] =     scheme -> noise2; /* 1 -> 0 */
      prob[3] = 1 - scheme -> noise2; /* 1 -> 1 */

      /* Add noise to the measurements */
      n = scheme -> nRounds;
      for (i = 0; i < n; i++)
      {  if (scheme -> measurements[i] == 0)
         {  if (sampleUniform() < scheme -> noise1)
               scheme -> measurements[i] = 1;
         }
         else
         {  if (sampleUniform() < scheme -> noise2)
               scheme -> measurements[i] = 0;
         }
      }

      /* Set the intial Fourier probability to zero */
      initFourierZero(scheme -> probability);

      /* Try all possible states */
      for (i = 0; i < n; i++)
      {  scheme -> state[i] = 0;
      }

      while (1)
      {  /* Compute the probability */
         initFourierProb(&dist,
                         scheme -> nRounds,
                         scheme -> k,
                         scheme -> beta,
                         scheme -> state);


         /* Compute the probability of the measurement given the state */
         p = 1;
         for (i = 0; i < n; i++)
         {  p *= prob[2 * scheme -> state[i] + scheme -> measurements[i]];
         }
         scaleFourier(&dist, p);

         /* Add the contribution */
         addFourier(scheme -> probability, &dist);

         /* Move to the next state */
         i = 0;
         while (i < n)
         {  if ((++(scheme -> state[i])) == 1)
            {  break;
            }

            scheme -> state[i] = 0; i ++;
         }
         if (i == n) break;
      }

      freeFourier(&dist);
   }
}



/* ======================================================================== */
/* SELECTION OF K VALUES                                                    */
/* ======================================================================== */

/* ------------------------------------------------------------------------ */
void pluginRandomBeta(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{
   sampleAngles(scheme -> beta, scheme -> nRounds);
}


/* ------------------------------------------------------------------------ */
void pluginAdaptiveK(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  double p;
   double kSum=0;
   int    i,k,n;

   n = density -> nEig;

   for (i = 0; i <= n-1; i++)
   {  kSum += density -> weights[i] * (1.25 / density->sigma[i]);
   }

   k = ceil(kSum);
   if (k > scheme -> kMax) k = scheme -> kMax;
   scheme -> k[0] = k;
}


/* ------------------------------------------------------------------------ */
void pluginAdaptiveCyclicK(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  double p;
   double kSum = 0;
   int    i,k,kMax,n;

   n = density -> nEig;

   for (i = 0; i <= n-1; i++)
   {  kSum += density -> weights[i] * (1.25 / density->sigma[i]);
   }

   kMax = ceil(kSum);
   if (kMax > scheme -> kMax) kMax = scheme -> kMax;

   k = scheme -> k[0] + 1;
   if (k > kMax) k = 1;
   scheme -> k[0] = k;
}


/* ------------------------------------------------------------------------ */
void pluginCyclicFixedK(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  int k;

   k = scheme -> k[0] + 1;
   if (k > scheme -> kMax) k = 1;
   scheme -> k[0] = k;
}


/* ------------------------------------------------------------------------ */
void pluginOLDCyclicK(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  int i,n,k,kMax;
   double wMax, sigma;

   wMax = 2;
   n = density -> nEig;
   for (i = 0; i < n; i++)
   {  if ((density -> weights[i] < wMax) && (density -> weights[i] > 0))
      {  wMax = density -> weights[i];
         sigma = density -> sigma[i];
      }
   }
   kMax = ceil(1.25 / sigma);
   if (kMax > scheme -> kMax) kMax = scheme -> kMax;


   k = scheme -> k[0] + 1;
   if (k > kMax) k = 1;
   scheme -> k[0] = k;
}


/* ------------------------------------------------------------------------ */
void pluginAdaptiveSampleK(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  double p;
   int i,n,k;

   p = sampleUniform();
   n = density -> nEig;

   for (i = 0; i <= n-1; i++)
   {  p -= density -> weights[i];
      if ((p < 0) || (i == n))
      {  k = ceil(1.25 / density->sigma[i]);
         if (k > scheme -> kMax) k = scheme -> kMax;
         scheme -> k[0] = k;
         break;
      }
   }
}


/* ------------------------------------------------------------------------ */
void pluginAlternateK2(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  scheme -> k[0] = 1 + (iteration % 2);
}


/* ------------------------------------------------------------------------ */
void pluginAlternateK3(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  scheme -> k[0] = 1 + (iteration % 3);
}


/* ------------------------------------------------------------------------ */
void pluginAlternateK5(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  scheme -> k[0] = 1 + (iteration % 5);
}


/* ------------------------------------------------------------------------ */
void pluginAlternateK10(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  scheme -> k[0] = 1 + (iteration % 10);
}


/* ------------------------------------------------------------------------ */
void pluginAlternateK20(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  scheme -> k[0] = 1 + (iteration % 20);
}


/* ------------------------------------------------------------------------ */
void pluginAlternateK50(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{ scheme -> k[0] = 1 + (iteration % 50);
}


/* ------------------------------------------------------------------------ */
void pluginAlternateK80(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  scheme -> k[0] = 1 + (iteration % 80);
}


/* ------------------------------------------------------------------------ */
void pluginAlternateK100(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  scheme -> k[0] = 1 + (iteration % 100);
}


/* ------------------------------------------------------------------------ */
void pluginAlternateK256(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  scheme -> k[0] = 1 + (iteration % 256);
}


/* ------------------------------------------------------------------------ */
void pluginAlternateK512(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  scheme -> k[0] = 1 + (iteration % 512);
}


/* ------------------------------------------------------------------------ */
void pluginAlternateK1024(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  scheme -> k[0] = 1 + (iteration % 1024);
}


/* ------------------------------------------------------------------------ */
void pluginAlternateK4096(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration)
/* ------------------------------------------------------------------------ */
{  scheme -> k[0] = 1 + (iteration % 4096);
}



/* ======================================================================== */
/* MEASUREMENT SCHEMES                                                      */
/* ======================================================================== */

/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeA(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 3);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> k[1] = 2;
   ptr -> k[2] = 5;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB2(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK2;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB3(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK3;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB5(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK5;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB10(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK10;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB20(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK20;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB50(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK50;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB256(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK256;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB512(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK512;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB1024(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK1024;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB4096(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK4096;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeC(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 50;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAdaptiveK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeC4096(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 4096;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAdaptiveK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeD(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 50;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAdaptiveCyclicK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeD4096(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 4096;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAdaptiveCyclicK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB10_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK10;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB20_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK20;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB50_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK50;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB80_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK80;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB100_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK100;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB10_Noisy_25(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kErr = 25;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK10;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB20_Noisy_25(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kErr = 25;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK20;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeB20_Noisy_50(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kErr = 50;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAlternateK20;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeC10_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 10;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAdaptiveK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeC20_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 20;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAdaptiveK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeC50_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 50;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAdaptiveK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeD10_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 10;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAdaptiveCyclicK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeD20_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 20;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAdaptiveCyclicK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeD50_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 50;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAdaptiveCyclicK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeFixed2(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 2;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeFixed3(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 3;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeX(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 50;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginOLDCyclicK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeXInf(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 4096;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginOLDCyclicK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeY(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 50;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAdaptiveK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeYInf(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 4096;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAdaptiveK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureRounds;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0;
   ptr -> noise2 = 0;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeX_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 50;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginOLDCyclicK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeX10_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 10;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginOLDCyclicK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeX20_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 20;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginOLDCyclicK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeX_Noisy30(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 30;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginOLDCyclicK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeX_Noisy40(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 40;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginOLDCyclicK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeX_Noisy50(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 50;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginOLDCyclicK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeX_Noisy60(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 60;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginOLDCyclicK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeY_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;

   newMeasurementScheme(scheme, 1);
   ptr = *scheme;

   /* Initialize the fields */
   ptr -> k[0] = 1;
   ptr -> kMax = 50;
   ptr -> kErr = 100;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = pluginAdaptiveSampleK;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureDecoherence;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE_Noisy_2(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.02;
   ptr -> noise2 = 0.02;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE_Noisy_5(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.05;
   ptr -> noise2 = 0.05;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE_Noisy_8(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.08;
   ptr -> noise2 = 0.08;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE_Noisy_10(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.10;
   ptr -> noise2 = 0.10;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE_Noisy_15(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.15;
   ptr -> noise2 = 0.15;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE_Noisy_20(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.2;
   ptr -> noise2 = 0.2;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE_Noisy_25(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.25;
   ptr -> noise2 = 0.25;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE_Noisy_26(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.26;
   ptr -> noise2 = 0.26;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE_Noisy_30(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.30;
   ptr -> noise2 = 0.30;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE_Noisy_35(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.35;
   ptr -> noise2 = 0.35;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE_Noisy_40(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.40;
   ptr -> noise2 = 0.40;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE_Noisy_45(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.45;
   ptr -> noise2 = 0.45;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE_Noisy_48(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.48;
   ptr -> noise2 = 0.48;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeE_Noisy_50(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k=1;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Set the k values */
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.5;
   ptr -> noise2 = 0.5;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = NULL;
   ptr -> fptrMeasure  = pluginMeasureFourier;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeF_Noisy(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Initialize the fields */
   /* Set the k values */
   k = 1;
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.02;
   ptr -> noise2 = 0.02;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureNoisy;
   ptr -> fptrFreeData = NULL;
}


/* ------------------------------------------------------------------------ */
void pluginMeasurementSchemeF_Noisy0001(BPEMeasurementScheme **scheme)
/* ------------------------------------------------------------------------ */
{  BPEMeasurementScheme *ptr;
   int i,k;

   newMeasurementScheme(scheme, 5);
   ptr = *scheme;

   /* Initialize the fields */
   /* Set the k values */
   k = 1;
   for (i = 0; i < 5; i++)
   {  ptr -> k[4-i] = k;
      k *= 2;
   }

   ptr -> noise1 = 0.0000001;
   ptr -> noise2 = 0.0000001;

   /* Set the measurement scheme */
   ptr -> fptrSetK     = NULL;
   ptr -> fptrSetBeta  = pluginRandomBeta;
   ptr -> fptrMeasure  = pluginMeasureNoisy;
   ptr -> fptrFreeData = NULL;
}



/* ======================================================================== */
/* INITIAL DENSITIES                                                        */
/* ======================================================================== */

/* ------------------------------------------------------------------------ */
void pluginBPEInitDensityB3_InitDensity(BPEInitDensity *self, BPEDensity *density)
/* ------------------------------------------------------------------------ */
{  int i;
   int s;

   density -> mu[0] = 2;
   density -> mu[1] = 4;
   density -> mu[2] = 5;

   for (i = 0; i < 3; i++)
   {  s = ((2 * (int)(2 * sampleUniform())) - 1);
      density -> mu[i]     += ((BPEInitDensityB3 *)self) -> initBias * s;
      density -> sigma[i]   = ((BPEInitDensityB3 *)self) -> initSigma;
      density -> fourier[i] = NULL;
   }
}


/* ------------------------------------------------------------------------ */
BPEInitDensity *pluginBPEInitDensityB3(double sigma, double bias)
/* ------------------------------------------------------------------------ */
{  BPEInitDensityB3 *plugin;

   plugin = malloc(sizeof(BPEInitDensityB3));
   plugin -> base.fptrInitDensity = pluginBPEInitDensityB3_InitDensity;
   plugin -> base.fptrDelete = NULL;

   plugin -> initSigma = sigma;
   plugin -> initBias  = bias;

   return (BPEInitDensity *)plugin;
}


/* ------------------------------------------------------------------------ */
void pluginInitWeightsB3(BPEConfig *config, double *weights)
/* ------------------------------------------------------------------------ */
{
   weights[0] = 0.5;
   weights[1] = 0.3;
   weights[2] = 0.2;
}


/* ------------------------------------------------------------------------ */
void pluginBPEInitDensityA_InitDensity(BPEInitDensity *self, BPEDensity *density)
/* ------------------------------------------------------------------------ */
{  int i;
   int s;

   density -> mu[0] = 5;

   s = ((2 * (int)(2 * sampleUniform())) - 1);
   density -> mu[0]     += ((BPEInitDensityA *)self) -> initBias * s;
   density -> sigma[0]   = ((BPEInitDensityA *)self) -> initSigma;
   density -> fourier[0] = NULL;

}


/* ------------------------------------------------------------------------ */
BPEInitDensity *pluginBPEInitDensityA(double sigma, double bias)
/* ------------------------------------------------------------------------ */
{  BPEInitDensityA *plugin;

   plugin = malloc(sizeof(BPEInitDensityA));
   plugin -> base.fptrInitDensity = pluginBPEInitDensityA_InitDensity;
   plugin -> base.fptrDelete = NULL;

   plugin -> initSigma = sigma;
   plugin -> initBias  = bias;

   return (BPEInitDensity *)plugin;
}


/* ------------------------------------------------------------------------ */
void pluginInitWeightsA(BPEConfig *config, double *weights)
/* ------------------------------------------------------------------------ */
{
   weights[0] = 1;
}
