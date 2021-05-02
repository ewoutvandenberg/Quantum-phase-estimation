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

typedef struct
{  int                    nTrials;
   BPEConfig              config;
   void (*fptrCreateGroundTruth)(BPEPhases **, int);
   void (*fptrCreateMeasurementScheme)(BPEMeasurementScheme **);
} BPEExperiment;



/* ------------------------------------------------------------------------ */
int generateExperiment(BPEExperiment *experiment, int index)
/* ------------------------------------------------------------------------ */
{  BPEConfig *config;

   /* Initialize the experiement */
   experiment -> nTrials                     = 20;
   experiment -> fptrCreateGroundTruth       = NULL;
   experiment -> fptrCreateMeasurementScheme = NULL;

   /* Default configuration parameters */
   config = &(experiment -> config);
   config -> nIter        = 1000*1000;
   config -> nEig         = -1;
   config -> kMax         = 200;
   config -> epsilon      = 1e-4;
   config -> mode         = BPE_MODE_MIXED;
   config -> displaySolve = 0;
   config -> updateDouble = 512;
   config -> curvilinear  = 1;
   config -> pluginInitDensity = NULL;
   config -> fptrInitWeights   = NULL;

   switch (index)
   {
      /* ------------------------------------------------------------ */
      /* Plot A - Convergence of Normal and hybrid priors (nEig=2)    */
      /* ------------------------------------------------------------ */

      case 1 : /* Normal mode only */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth2; /* 2 */

               config -> mode  = BPE_MODE_NORMAL_ONLY;
               config -> nEig  = 2;
               config -> nIter = 10*1000*1000;
               break;

      case 2 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth2; /* 2 */

               config -> mode  = BPE_MODE_MIXED;
               config -> nEig  = 2;
               config -> nIter = 10*1000*1000;
               break;


      /* ------------------------------------------------------------ */
      /* Plot B - Convergence of Fourier prior with nEig=2            */
      /* ------------------------------------------------------------ */

      case 3 : /* Fourier mode only - kMax = 200 */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth2; /* 2 */

               config -> mode  = BPE_MODE_FOURIER_ONLY;
               config -> nEig  = 2;
               config -> kMax  = 200;
               break;

      case 4 : /* Fourier mode only - kMax = 1000 */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth2; /* 2 */

               config -> mode  = BPE_MODE_FOURIER_ONLY;
               config -> nEig  = 2;
               config -> kMax  = 1000;
               break;

      case 5 : /* Fourier mode only - kMax = 5000 */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth2; /* 2 */

               config -> mode  = BPE_MODE_FOURIER_ONLY;
               config -> nEig  = 2;
               config -> kMax  = 5000;
               break;


      /* ------------------------------------------------------------ */
      /* Plot C - Convergence of Normal and hybrid priors (nEig=3)    */
      /* ------------------------------------------------------------ */

      case 6 : /* Normal mode only */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 7 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;


      /* ------------------------------------------------------------ */
      /* Plot D - Convergence of Normal and hybrid priors (nEig=3)    */
      /*          using five distributions.                           */
      /* ------------------------------------------------------------ */

      case 8 : /* Normal mode only */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_NORMAL_ONLY;
               config -> nEig  = 5;
               config -> nIter = 10*1000*1000;
               break;

      case 9 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_MIXED;
               config -> nEig  = 5;
               config -> nIter = 10*1000*1000;
               break;


      /* ------------------------------------------------------------ */
      /* Experiment: Effect of kMax in test problem B                 */
      /* ------------------------------------------------------------ */

      case 10 : /* Fourier mode only */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_FOURIER_ONLY;
               config -> nEig  = 3;
               config -> nIter = 100*1000;
               config -> kMax  = 200;
               break;

      case 11 : /* Fourier mode only */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_FOURIER_ONLY;
               config -> nEig  = 3;
               config -> nIter = 100*1000;
               config -> kMax  = 1000;
               break;

      case 12 : /* Fourier mode only */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_FOURIER_ONLY;
               config -> nEig  = 3;
               config -> nIter = 100*1000;
               config -> kMax  = 5000;
               break;

      case 14 : /* Fourier mode only */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_MIXED;
               config -> nEig  = 3;
               config -> nIter = 100*1000;
               config -> kMax  = 200;
               break;

      case 15 : /* Fourier mode only */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_MIXED;
               config -> nEig  = 3;
               config -> nIter = 100*1000;
               config -> kMax  = 1000;
               break;

      case 16 : /* Fourier mode only */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_MIXED;
               config -> nEig  = 3;
               config -> nIter = 100*1000;
               config -> kMax  = 5000;
               break;

      case 17 : /* Fourier mode only - kMax = 200 (similar to case 2) */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth2; /* 2 */

               config -> mode  = BPE_MODE_MIXED;
               config -> nEig  = 2;
               config -> kMax  = 200;
               break;

      case 18 : /* Fourier mode only - kMax = 1000 */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth2; /* 2 */

               config -> mode  = BPE_MODE_MIXED;
               config -> nEig  = 2;
               config -> kMax  = 1000;
               break;

      case 19 : /* Fourier mode only - kMax = 5000 */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth2; /* 2 */

               config -> mode  = BPE_MODE_MIXED;
               config -> nEig  = 2;
               config -> kMax  = 5000;
               break;


      /* ------------------------------------------------------------ */
      /* Experiments: different weight solve methods                  */
      /* ------------------------------------------------------------ */

      case 20 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_NORMAL_ONLY;
               config -> nEig  = 5;
               config -> nIter = 10*1000*1000;
               config -> updateDouble = 2;
               config -> curvilinear = 1;
               config -> displaySolve = 1;
               break;

      case 21 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_NORMAL_ONLY;
               config -> nEig  = 5;
               config -> nIter = 10*1000*1000;
               config -> updateDouble = 64;
               config -> curvilinear = 1;
               config -> displaySolve = 1;
               break;

      case 22 : /* Normal mode (see also 8) */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_NORMAL_ONLY;
               config -> nEig  = 5;
               config -> nIter = 10*1000*1000;
               config -> updateDouble = 512;
               config -> curvilinear = 1;
               config -> displaySolve = 1;
               break;

      case 23 : /* Normal mode  */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_NORMAL_ONLY;
               config -> nEig  = 5;
               config -> nIter = 10*1000*1000;
               config -> updateDouble = 512;
               config -> curvilinear = 0; /* Backtracking line search only */
               config -> displaySolve = 1;
               break;

      case 24 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_NORMAL_ONLY;
               config -> nEig  = 5;
               config -> nIter = 10*1000*1000;
               config -> updateDouble = 100*1000;
               config -> curvilinear = 1;
               config -> displaySolve = 1;
               break;

      case 25 : /* Mixed mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_MIXED;
               config -> nEig  = 5;
               config -> nIter = 10*1000*1000;
               config -> updateDouble = 2;
               config -> curvilinear = 1;
               config -> displaySolve = 1;
               break;

      case 26 : /* Mixed mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_MIXED;
               config -> nEig  = 5;
               config -> nIter = 10*1000*1000;
               config -> updateDouble = 64;
               config -> curvilinear = 1;
               config -> displaySolve = 1;
               break;

      case 27 : /* Mixed mode (see also 9) */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_MIXED;
               config -> nEig  = 5;
               config -> nIter = 10*1000*1000;
               config -> updateDouble = 512;
               config -> curvilinear = 1;
               config -> displaySolve = 1;
               break;

      case 28 : /* Mixed mode  */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_MIXED;
               config -> nEig  = 5;
               config -> nIter = 10*1000*1000;
               config -> updateDouble = 512;
               config -> curvilinear = 0; /* Backtracking line search only */
               config -> displaySolve = 1;
               break;

      case 29 : /* Mixed mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeA;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3; /* 3 */

               config -> mode  = BPE_MODE_MIXED;
               config -> nEig  = 5;
               config -> nIter = 10*1000*1000;
               config -> updateDouble = 100*1000;
               config -> curvilinear = 1;
               config -> displaySolve = 1;
               break;


      /* ------------------------------------------------------------ */
      /* Experiment G1: cyclic K, single eigenphase                   */
      /* ------------------------------------------------------------ */

      case 30 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 31 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 32 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB2;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 33 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB2;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 34 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB3;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 35 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB3;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 36 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB5;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 37 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB5;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 38 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB10;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 39 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;


      /* ------------------------------------------------------------ */
      /* Experiment G4: cyclic K, three eigenphases                   */
      /* ------------------------------------------------------------ */

      case 40 : /* Normal mode - fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 41 : /* Mixed mode - fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 42 : /* Normal mode - fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB2;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 43 : /* Mixed mode - initially converges then diverges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB2;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 44 : /* Normal mode - fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB3;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 55 : /* Mixed mode - initially converges then diverges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB3;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 46 : /* Normal mode - fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB5;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 47 : /* Mixed mode - initially converges then diverges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB5;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 48 : /* Mixed mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB10;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 49 : /* Mixed mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 50 : /* Mixed mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB5;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 1000;
               break;

      case 51 : /* Mixed mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB5;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 5000;
               break;

      case 52 : /* Mixed mode - 5 eigenphase - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB10;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 5;
               break;

      case 53 : /* Mixed mode - 5 eigenphases - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 5;
               break;

/*    case 54 : // Mixed mode
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB5;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 5000;
               config -> nIter = 7*1000*1000;
               break;
*/


      /* ------------------------------------------------------------ */
      /* Experiments G4: dynamic update of k, single eigenphase       */
      /* ------------------------------------------------------------ */

      case 60 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 61 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 62 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeC;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 63 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeC;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 64 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeD;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 65 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeD;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 66 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB50;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 67 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB50;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 68 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeC4096;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 69 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeC4096;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 70 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeD4096;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 71 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeD4096;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 72 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB1024;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 73 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB4096;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;


      /* ------------------------------------------------------------ */
      /* Experiments G5: dynamic update of k, three eigenphases (3/3) */
      /* ------------------------------------------------------------ */

      case 80 : /* Fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 81 : /* Fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 82 : /* Purely adaptvie - Fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeC;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 83 : /* Purely adaptive - Fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeC;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 84 : /* Adaptive cyclic - Fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeD;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 85 : /* Adaptive cyclic - Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeD;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 86 : /* Fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB50;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 87 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB50;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;


      /* ------------------------------------------------------------ */
      /* Experiments G5-B: dynamic update of k, tree eigenphase (3/5) */
      /* ------------------------------------------------------------ */

      case 90 : /* Fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 5;
               break;

      case 91 : /* Pure adaptive - Fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeC;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 5;
               break;

      case 92 : /* Adaptive cyclic - Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeD;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 5;
               break;

      case 93 : /* Converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB50;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 5;
               break;



      /* ============================================================ */
      /*                 PHASE ESTIMATION WITH NOISE                  */
      /* ============================================================ */

      /* ------------------------------------------------------------ */
      /* Experiments N1: decoherence, single eigenphase               */
      /* ------------------------------------------------------------ */

      case 100 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB50_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 101 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB50_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 102 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeC50_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 103 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeC50_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 104 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeD50_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 105 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeD50_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 106 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 107 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;


      /* ------------------------------------------------------------ */
      /* Experiments N2: decoherence, three eigenphases               */
      /* ------------------------------------------------------------ */

      case 110 : /* Cyclic - Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB10_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 111 : /* Cyclic - Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 112 : /* Cyclic - Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB50_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 113 : /* Purely adaptive - Fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeC50_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 114 : /* Adaptive cyclic - Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeD50_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 115 : /* Purely adaptive - Fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeC20_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 116 : /* Adaptive cyclic - Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeD20_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 117 : /* Purely adaptive - Fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeC10_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 118 : /* Adaptive cyclic - Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeD10_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 119 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB50_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 1000;
               break;

      case 120 : /* Partial convergence */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB80_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 121 : /* Partial convergence */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB100_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 122 : /* Succeeds (kMax=1000) */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB80_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 1000;
               break;

      case 123 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB100_Noisy;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 1000;
               break;


      /* ------------------------------------------------------------ */
      /* Experiments I0: Fourier measurements, exact eigenphase       */
      /* ------------------------------------------------------------ */

      case 130 : /* Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth4;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 131 : /* Fourier & Normal mode */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth4;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;


      case 132 : /* - */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_2;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth4;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 133 : /* - */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_2;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth4;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;


      /* ------------------------------------------------------------ */
      /* Experiments I1: Fourier measurements, single eigenphase      */
      /* ------------------------------------------------------------ */

      case 140 : /* Noiseless - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 141 : /* Noiseless - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 142 : /* 2% noise - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_2;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 143 : /* 2% noise - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_2;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 144 : /* 5% noise - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_5;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 145 : /* 5% noise - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_5;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 146 : /* 20% noise - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 147 : /* 20% noise - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 148 : /* 25% noise - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_25;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 149 : /* 25% noise - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_25;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;


      case 150 : /* 26% noise - only just converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_26;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 151 : /* 26% noise - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_26;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 152 : /* 30% noise - fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_30;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 153 : /* 30% noise - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_30;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 154 : /* 40% noise - fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_40;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 155 : /* 40% noise - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_40;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 156 : /* 45% noise - fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_45;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 157 : /* 45 noise - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_45;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;

      case 158 : /* 48% noise - fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_48;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 1;
               break;

      case 159 : /* 48% noise - converges */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_48;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth1;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1;
               break;


      /* ------------------------------------------------------------ */
      /* Experiments I2: Fourier measurements, three eigenphases      */
      /* ------------------------------------------------------------ */

      case 160 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 161 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 162 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_2;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 163 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_2;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 164 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_5;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 165 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_5;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 166 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_8;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 167 : /* Partially fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_8;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 168 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_8;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 1000;
               break;

      case 169 : /* Partially fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_10;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 170 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_10;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 1000;
               break;

      case 171 : /* Fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_15;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 172 : /* Partially fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_15;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 173 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_15;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 1000;
               break;

      case 174 : /* Fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_NORMAL_ONLY;
               config -> nEig = 3;
               break;

      case 175 : /* Partially fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               break;

      case 176 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 1000;
               break;

      case 177 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_25;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 1000;
               break;

      case 178 : /* Partially fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_30;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 1000;
               break;

      case 179 : /* Succeeds */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_30;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 2000;
               break;

      case 180 : /* Partially fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_35;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 1000;
               break;

      case 181 : /* Partially fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_35;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 2000;
               break;

      case 182 : /* Largely fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_40;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 1000;
               break;

      case 183 : /* Fails */
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeE_Noisy_45;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruth3;

               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3;
               config -> kMax = 1000;
               break;


      /* ------------------------------------------------------------ */
      /* Experiments with many phases                                 */
      /* ------------------------------------------------------------ */
      case 201 :
      case 201+20 :
      case 201+40 :
      case 201+60 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthMultiple1;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 1 + (index - 201) / 20;
               break;

      case 202 :
      case 202+20 :
      case 202+40 :
      case 202+60 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthMultiple2;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 2 + (index - 202) / 20;
               break;

      case 203 :
      case 203+20 :
      case 203+40 :
      case 203+60 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthMultiple3;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 3 + (index - 203) / 20;
               break;

      case 204 :
      case 204+20 :
      case 204+40 :
      case 204+60 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthMultiple4;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 4 + (index - 204) / 20;
               break;

      case 205 :
      case 205+20 :
      case 205+40 :
      case 205+60 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthMultiple5;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 5 + (index - 205) / 20;
               break;

      case 206 :
      case 206+20 :
      case 206+40 :
      case 206+60 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthMultiple6;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 6 + (index - 206) / 20;
               break;

      case 207 :
      case 207+20 :
      case 207+40 :
      case 207+60 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthMultiple7;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 7 + (index - 207) / 20;
               break;

      case 208 :
      case 208+20 :
      case 208+40 :
      case 208+60 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthMultiple8;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 8 + (index - 208) / 20;
               break;

      case 209 :
      case 209+20 :
      case 209+40 :
      case 209+60 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthMultiple9;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 9 + (index - 209) / 20;
               break;

      case 210 :
      case 210+20 :
      case 210+40 :
      case 210+60 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthMultiple10;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 10 + (index - 210) / 20;
               break;

      case 211 :
      case 211+20 :
      case 211+40 :
      case 211+60 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthMultiple11;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 11 + (index - 211) / 20;
               break;

      case 212 :
      case 212+20 :
      case 212+40 :
      case 212+60 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthMultiple12;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 12 + (index - 212) / 20;
               break;



      /* ------------------------------------------------------------ */
      /* Fixed problems with many phases                              */
      /* ------------------------------------------------------------ */
      case 301 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthFixed12;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 14;
               break;

      case 302 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthFixed7;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 10;
               break;

      case 303 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateGroundTruthFixed10;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 12;
               break;


      /* ------------------------------------------------------------ */
      /* Problems with spurious eigenstates                           */
      /* ------------------------------------------------------------ */
      case 310 :
      case 311 :
      case 312 :
      case 313 :
      case 314 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateNoisyGroundTruth3_10_2;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 4+(index-310);
               break;

      case 315 :
      case 316 :
      case 317 :
      case 318 :
      case 319 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateNoisyGroundTruth3_10_5;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 4+(index-315);
               break;

      case 320 :
      case 321 :
      case 322 :
      case 323 :
      case 324 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateNoisyGroundTruth3_10_10;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 4+(index-320);
               break;

      case 325 :
      case 326 :
      case 327 :
      case 328 :
      case 329 :
               experiment -> fptrCreateMeasurementScheme = pluginMeasurementSchemeB20;
               experiment -> fptrCreateGroundTruth       = pluginCreateNoisyGroundTruth3_10_20;
               experiment -> nTrials                     = 100;
               config -> mode = BPE_MODE_MIXED;
               config -> nEig = 4+(index-325);
               break;

      default :
               return -1;
   }

   return 0;
}


/* ------------------------------------------------------------------------ */
int main(int argc, const char *argv[])
/* ------------------------------------------------------------------------ */
{  BPEExperiment         experiment;
   BPEMeasurementScheme *scheme;
   BPEPhases            *phases;
   int                   index, trial, i;

   #ifdef DEBUG_RANDOM
   initRandomBuffer();
   #endif

   /* Parse parameters */
   index = (argc == 2) ? atoi(argv[1]) : -1;

   /* Generate the experiment */
   if (generateExperiment(&experiment, index) == -1)
   {  printf("Usage: %s <index>\n", argv[0]);
      exit(1);
   }

   /* Run different trials */
   for (trial = 0; trial < experiment.nTrials; trial++)
   {
      /* Initialize the random number generator per trial */
      /* to ensure the starting states are identical for  */
      /* each of the different configurations.            */
      srand(trial);

      /* Create the measurement scheme and ground truth */
      experiment.fptrCreateMeasurementScheme(&scheme);
      experiment.fptrCreateGroundTruth(&(scheme -> phases), trial);

      /* Configure the measurement scheme */

      /* Print start of new experiment */
      phases = scheme -> phases;
      printf("# %d", phases -> nEig);
      for (i = 0; i < phases -> nEig; i++) printf(" %.12e", phases->phi[i]);
      for (i = 0; i < phases -> nEig; i++) printf(" %.12e", phases->w[i]);
      printf("\n");

      /* Bayesian phase estimation */
      BayesianPhaseEstimation(scheme, &(experiment.config));

      /* Delete the measurment scheme */
      deleteMeasurementScheme(scheme);
   }

   /* Free BPEInitDensity in config, if needed */
   if (experiment.config.pluginInitDensity)
   {  free(experiment.config.pluginInitDensity);
   }

   #ifdef DEBUG_RANDOM
   freeRandombuffer();
   #endif

   return 0;
}
