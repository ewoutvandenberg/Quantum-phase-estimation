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

#include "bpe.h"
#include "bpe_plugins.h"
#include <stdlib.h>  /* malloc */


/* ------------------------------------------------------------------------ */
int main(int argc, const char *argv[])
/* ------------------------------------------------------------------------ */
{  BPEConfig             config;
   BPEMeasurementScheme *scheme;
   BPEPhases            *phases;
   double                sigma;
   int                   index, trial, biasIdx, i;
   void                 (*pluginScheme)(BPEMeasurementScheme **);


   /* Parse parameters */
   index = (argc == 2) ? atoi(argv[1]) : -1;

   /* Default configuration parameters */
   config.nIter        = 100*1000;
   config.nEig         = 1;
   config.mode         = BPE_MODE_NORMAL_ONLY;
   config.displaySolve = 0;
   config.updateDouble = 512;
   config.curvilinear  = 1;
   config.pluginInitDensity = NULL;
   config.fptrInitWeights   = pluginInitWeightsA;

   srand(0);

   switch (index)
   {  case 1 : pluginScheme = pluginMeasurementSchemeB;
               sigma        = fourierCriticalSigma(200, 1e-4);
               break;

      case 2 : pluginScheme = pluginMeasurementSchemeB2;
               sigma        = fourierCriticalSigma(200, 1e-4);
               break;

      case 3 : pluginScheme = pluginMeasurementSchemeB5;
               sigma        = fourierCriticalSigma(200, 1e-4);
               break;

      case 4 : pluginScheme = pluginMeasurementSchemeB10;
               sigma        = fourierCriticalSigma(200, 1e-4);
               break;

      case 5 : pluginScheme = pluginMeasurementSchemeB20;
               sigma        = fourierCriticalSigma(200, 1e-4);
               break;

      case 6 : pluginScheme = pluginMeasurementSchemeB10;
               sigma        = fourierCriticalSigma(1000, 1e-4);
               break;

      case 7 : pluginScheme = pluginMeasurementSchemeB10;
               sigma        = fourierCriticalSigma(5000, 1e-4);
               break;

      default :
              printf("Usage: %s <index>\n", argv[0]);
             exit(1);
   }

   /* For different values run */
   for (biasIdx = 0; biasIdx <= 200; biasIdx ++)
   {  double beta = 0.01 * biasIdx; /*(biasIdx*biasIdx) * 0.01;*/

      config.pluginInitDensity = pluginBPEInitDensityA(sigma, beta);

      /* Run different trials */
      for (trial = 0; trial < 100; trial++)
      {
         /* Initialize the random number generator per trial */
         /* to ensure the starting states are identical for  */
         /* each of the different configurations.            */

         /* Create the measurement scheme and ground truth */
         pluginScheme(&scheme);
         pluginCreateGroundTruth1(&(scheme -> phases), trial);

         /* Print start of new experiment */
         phases = scheme -> phases;
         printf("# %d", phases -> nEig);
         for (i = 0; i < phases -> nEig; i++) printf(" %.12e", phases->phi[i]);
         for (i = 0; i < phases -> nEig; i++) printf(" %.12e", phases->w[i]);
         printf("\n");

         /* Bayesian phase estimation */
         BayesianPhaseEstimation(scheme, &config);

         /* Delete the measurment scheme */
         deleteMeasurementScheme(scheme);
      }

      /* Free density initialization plugin */
      free(config.pluginInitDensity);
   }

   return 0;
}
