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
from   matplotlib import colors
import numpy as np
import bpe

from generic import *


# ----------------------------------------------------------------------
# Plot G1 - Cyclic exponent, single eigenphase
# ----------------------------------------------------------------------

# Load the results
instancesA = loadInstances('./cache/experiment_hybrid_30.dat') # Normal, kMax = 1
instancesB = loadInstances('./cache/experiment_hybrid_36.dat') # Normal, kMax = 5
instancesC = loadInstances('./cache/experiment_hybrid_31.dat') # Mixed, kMax = 1
instancesD = loadInstances('./cache/experiment_hybrid_37.dat') # Mixed, kMax = 5
instancesE = loadInstances('./cache/experiment_hybrid_39.dat') # Mixed, kMax = 20

# Display the averate switch iterations along with standard deviation
for instances in [instancesC,instancesD,instancesE] :
   data  = [problem.switchIter for problem in instances]
   width = max([len(s) for s in data])
   switches = np.zeros((len(data),width),dtype=int)
   for (i,d) in enumerate(data) :
      switches[i,:len(d)] = d

   print("Mean switch: %s (%s)" % (", ".join(["%d"%f for f in np.mean(switches,axis=0)]),
                                   ", ".join(["%.2f"%f for f in np.std(switches,axis=0)])))


# Plot
plotCollatedMu(instancesA,'-',alpha=0.4,color='C1')
plotCollatedMu(instancesB,'-',alpha=1.0,color='C1')
plotCollatedMu(instancesC,'-',alpha=0.3)
plotCollatedMu(instancesD,'-',alpha=0.6)
plotCollatedMu(instancesE,'-',alpha=1.0)


for (alpha,instances) in zip([0.3,0.6,1.0],[instancesC,instancesD,instancesE]) :
   problem = instances[0]
   error = np.zeros((len(instances),problem.phiMu.shape[0],problem.phiStar.size),dtype=np.double)
    
   for i in range(len(instances)) :
      problem = instances[i]
      d = np.abs(problem.collatedPhi - problem.phiStar)
      idx = d > np.pi
      d[idx] = 2*np.pi - d[idx]
      error[i,:,:] = d


   medianError = np.nanmedian(error,axis=0)
   switch = np.mean([problem.switchIter for problem in instances])
   idx1   = np.argwhere(problem.phiIter <= switch); idx1 = idx1[-1]
   idx2   = np.argwhere(problem.phiIter >  switch); idx2 = idx2[0]
   value  = ((problem.phiIter[idx2] - switch) * medianError[idx1] +
             (switch - problem.phiIter[idx1]) * medianError[idx2]) / (problem.phiIter[idx2] - problem.phiIter[idx1])


   c = colors.to_rgba('C0')
   c = tuple([(alpha * c[i] + (1-alpha) * 1) for i in range(3)])
   c = [0.75 * c[i] for i in range(3)]
   plt.plot(switch, value, color=[1,1,1], marker='.',markersize=15)
   plt.plot(switch, value, color=c, marker='.',markersize=10)


plt.text(1.24e2,1e-1,'A',fontweight='medium', fontsize=fontsize)
plt.text(1.15e6,5.0e-1,'B',fontweight='medium', fontsize=fontsize)
plt.text(1.15e6,7.0e-4,'C',fontweight='medium', fontsize=fontsize)
plt.text(1.15e6,1.3e-4,'D',fontweight='medium', fontsize=fontsize)
plt.text(1.15e6,5.0e-5,'E',fontweight='medium', fontsize=fontsize)
plt.legend(['A','B','C','D','E'],
           ['Normal ($c_{\max}$ = 1)','Normal ($c_{\max}$ = 5)',
            'Mixed ($c_{\max}$ = 1)','Mixed ($c_{\max}$ = 5)','Mixed ($c_{\max}$ = 20)'],
           handler_map={str: StringObjectHandler()}, loc='lower left', fontsize=12)
plt.ylim([1e-5,5e0])
plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Median phase error', fontsize=fontsize)
setTickFontsize(fontsize)
exportFigure("FigExperimentHybrid_G1", label='4a')


# ----------------------------------------------------------------------
# Plot G2 - Phase transitions
# ----------------------------------------------------------------------

epsilon = 1e-4; kMax = 200
sigmaCritical = bpe.DensityFourier.normalCriticalSigma(None, epsilon, kMax=kMax)
print(sigmaCritical)


for i in [1,2,3,4,6,7] :
   filename = './cache/experiment_transition_%d.dat' % i
   instances = loadInstances(filename)
   groups = []
   for i in range(len(instances) // 100) :
      groups.append(instances[i*100:(i+1)*100])
  
   results = np.zeros(len(groups))
   for i in range(len(groups)) :
      instances = groups[i]
      for j in range(len(instances)) :
         error = np.abs(instances[j].phiMu[-1] - instances[j].phiStar[0])
         error = np.minimum(error, 2*np.pi - error)
         if (error >= 1e-1) :
            results[i] += 1

   plt.plot(np.arange(len(results)) * 0.01, 100 * results / 100)
  

plt.text(1.42, 64.4,'A',fontweight='medium', fontsize=fontsize)
plt.text(0.98, 64.4,'B',fontweight='medium', fontsize=fontsize)
plt.text(0.70, 72.0,'C',fontweight='medium', fontsize=fontsize)
plt.text(0.49, 81.0,'D',fontweight='medium', fontsize=fontsize) # 0.49, 82.0
plt.text(0.32, 82.5,'E',fontweight='medium', fontsize=fontsize) # 0.31, 83.0
plt.text(0.09, 80.0,'F',fontweight='medium', fontsize=fontsize)
plt.legend(['A','B','C','D','E','F'],
           ['$\sigma_{\epsilon}(200)$, $c_{\max}$ = 1','$\sigma_{\epsilon}(200)$, $c_{\max}$ = 2','$\sigma_{\epsilon}(200)$, $c_{\max}$ = 5',
            '$\sigma_{\epsilon}(200)$, $c_{\max}$ = 10','$\sigma_{\epsilon}(1000)$, $c_{\max}$ = 10','$\sigma_{\epsilon}(5000)$, $c_{\max}$ = 10'],
           handler_map={str: StringObjectHandler()}, loc='lower right', fontsize=12)
plt.xlabel('Bias', fontsize=fontsize)
plt.ylabel('Failure probability (%)', fontsize=fontsize)
setTickFontsize(fontsize)

exportFigure("FigExperimentHybrid_G2", label='4b')
