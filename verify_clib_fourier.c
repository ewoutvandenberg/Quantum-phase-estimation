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

int main(int argc, const char *argv[])
{  TFourier *dist;
   int index;

   /* Combine unit tests */
   index = (argc >= 2) ? atoi(argv[1]) : 0;

   dist = newFourier();

   switch (index)
   {
      case 0 :
      {  initFourierZero(dist);
         printFourierCoef(dist);
         break;
      }

      case 1 :
      {  initFourierOne(dist);
         printFourierCoef(dist);
         break;
      }

      case 2 :
      {  initFourierUniform(dist);
         printFourierCoef(dist);
         break;
      }

      case 3 :
      {  int    k1[]    = {0,2,4};
         double beta1[] = {0.1,2.3,-1.7}; 
         int    k2[]    = {1,3};
         double beta2[] = {0.2,0.9};
         initFourierCosine(dist,3,k1,beta1,NULL);
         printFourierCoef(dist);
         initFourierSine(dist,2,k2,beta2,NULL);
         printFourierCoef(dist);
         break;
      }

      case 4 :
      {  int    k[]       = {0,1,3};
         double beta[]    = {0.1,2.3,-1.7};
         double weights[] = {1,2,3};
         initFourierCosine(dist,3,k,beta,weights);
         printFourierCoef(dist);
         initFourierSine(dist,3,k,beta,weights);
         printFourierCoef(dist);
         break;
      }

      case 5 :
      {  int k[] = {2};
         double beta[] = {0.3};
         initFourierProb(dist, 1, k, beta, NULL);
         printFourierCoef(dist);
         break;
      }

      case 6 :
      {  int k[] = {1,1};
         double beta[] = {0.3,2.7};
         initFourierProb(dist, 2, k, beta, NULL);
         printFourierCoef(dist);
         break;
      }

      case 7 :
      {  int k[] = {1,2,1,5};
         double beta[] = {0.3,2.7,-9.2,1.2};
         initFourierProb(dist, 4, k, beta, NULL);
         printFourierCoef(dist);
         break;
      }

      case 8 :
      {  int k[] = {1,2,5};
         double beta[] = {0.3,2.7,-9.2};
         initFourierProb(dist, 3, k, beta, NULL);
         printFourierCoef(dist);
         break;
      }

      case 9 :
      {  initFourierNormal(dist, 10, 0.4, 0.8);
         printFourierCoef(dist);
         break;
      }

      case 10 :
      {  double delta = 2*M_PI / 100.;
         int i;
         initFourierNormal(dist, 10, 0.4, 0.8);
         for (i = 0; i < 100; i++)
         {  printf("%.12e ", evaluateFourier(dist, i*delta));
         }
         printf("\n");
         break;
      }

      case 11 :
      {  int k[] = {1,2,1,5};
         double beta[] = {0.3,2.7,-9.2,1.2};
         initFourierProb(dist, 4, k, beta, NULL);
         printf("%.12e\n", getFourierMean(dist));
         printf("%.12e\n", getFourierVar(dist));
         printf("%.12e\n", getFourierStd(dist));
         break;
      }

      case 12 :
      {  printf("%.12e\n", fourierNormalError(7,0.2));
         break;
      }

      case 13 :
      {  printf("%.12e\n", fourierCriticalSigma(7,1e-5));
         break;
      }

      case 14 :
      {  TFourier *dist2;
   
         dist2 = newFourier();
         initFourierOne(dist);
         initFourierOne(dist2);
         multiplyFourier(dist, dist2);
         printFourierCoef(dist);

         deleteFourier(dist2);
         break;
      }

      case 15 :
      {  int k1[] = {1};
         double beta1[] = {0.3};
         TFourier *dist2;

         dist2 = newFourier();
         initFourierProb(dist, 1, k1, beta1, NULL);
         initFourierOne(dist2);
         multiplyFourier(dist, dist2);
         printFourierCoef(dist);

         deleteFourier(dist2);
         break;
      }

      case 16 :
      {  int k1[] = {1};
         int k2[] = {1};
         double beta1[] = {0.3};
         double beta2[] = {2.5};
         TFourier *dist2;
   
         dist2 = newFourier();
         initFourierProb(dist,  1, k1, beta1, NULL);
         initFourierProb(dist2, 1, k2, beta2, NULL);
         multiplyFourier(dist, dist2);
         printFourierCoef(dist);

         deleteFourier(dist2);
         break;
      }


      case 17 :
      {  int k1[] = {1,5,8};
         double beta1[] = {0.3,0.2,1.8};
         TFourier *dist2;
   
         dist2 = newFourier();
         initFourierProb(dist, 3, k1, beta1, NULL);
         initFourierNormal(dist2,14, 0.4, 0.8);
         multiplyFourier(dist, dist2);
         printFourierCoef(dist);

         deleteFourier(dist2);
         break;
      }
   }

   deleteFourier(dist);

   return 0;
}
