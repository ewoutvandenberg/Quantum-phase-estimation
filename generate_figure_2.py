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

import matplotlib.pyplot as plt
import scipy.special
import numpy as np
import bpe
import generic


plt.rc('font', size=16)

phi   = np.linspace(0,2*np.pi,10001)
sigma = 0.1

for (kMax,label) in zip([12,15,20,50],['2c','2d','2a','2b']) :
   mu = np.pi
   coef = bpe.DensityFourier()
   coef.initNormal(mu,sigma,kMax)
   p1 = coef.evaluate(phi)

   # Wrap-around function
   coef2 = bpe.DensityNormal(mu,sigma)
   p2 = coef2.evaluate(phi)

   # Compute the error   
   err = p2-p1

   errBound = coef.normalError(sigma,kMax)
   
   #f = np.exp(-(((phi-mu)/sigma)**2)/2) / (sigma*np.sqrt(2*np.pi))
   plt.plot(phi,err,color=[0.8,0.8,0.8],linestyle='-')
   plt.plot(phi,err,color=[0.2,0.2,1.0],linestyle=':')
   
   plt.plot(phi,p2,color=[0.9,0.5,0.2],linestyle='-')
   plt.plot(phi,p1,'b-')
   plt.plot([0,2*np.pi],[errBound,errBound],color=[0.5,0.5,0.5],linestyle='--')
   plt.plot([0,2*np.pi],[-errBound,-errBound],color=[0.5,0.5,0.5],linestyle='--')
   plt.xlabel('$\phi$')
   #plt.tight_layout()
   plt.subplots_adjust(bottom=0.15)

   generic.exportFigure("FigTruncatedFourier%d"%(kMax), label=label)
