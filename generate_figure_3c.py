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

import numpy as np
import scipy.special
import matplotlib.pyplot as plt


import generic


def evalError1(sigma,kMax) :
   x = 0.
   termMax=max(kMax,100000)
   for k in range(termMax,kMax,-1) :
      x = x + np.exp(-((sigma*k)**2)/2)

   x = x / np.pi
   x = x  + scipy.special.erfc(termMax*sigma / np.sqrt(2)) / (sigma * np.sqrt(2*np.pi))
   return x

def evalError2(sigma,kMax) :
   return scipy.special.erfc(kMax*sigma / np.sqrt(2)) / (sigma * np.sqrt(2*np.pi))


sigma = np.logspace(0,-6,101)
error1 = np.zeros(sigma.size)
error2 = np.zeros(sigma.size)
fontsize = 14

legendKey = []
legendStr = []
for [j,kMax] in enumerate([5000,1000,200,50,10]):
   for idx in range(sigma.size) :
      error1[idx] = evalError1(sigma[idx], kMax)
      error2[idx] = evalError2(sigma[idx], kMax)

   plt.plot(sigma,error1,'k:',alpha=0.7)
   plt.plot(sigma,error2,'-',color='C%d'%j)

   idx = int(np.min(np.argwhere(error2 > 1e-10)))
   plt.text(1.2*sigma[idx],1.4e-10,chr(ord('A')+j),fontweight='medium', fontsize=fontsize)
   legendKey.append(chr(ord('A')+j))
   legendStr.append('Max. coef = %d' % kMax)


plt.xscale('log')
plt.yscale('log')
plt.xlim([3e-7,1.2])
plt.ylim([1e-10,1e6])
plt.xlabel('$\sigma$', fontsize=fontsize)
plt.ylabel('Fourier truncation error', fontsize=fontsize)
#plt.legend(loc='lower left', fontsize=fontsize)
plt.legend(legendKey, legendStr,
           handler_map={str: generic.StringObjectHandler()}, loc='lower left', fontsize=fontsize)
generic.setTickFontsize(fontsize)
plt.tight_layout()
generic.exportFigure("FigFourierError", label='3c')

