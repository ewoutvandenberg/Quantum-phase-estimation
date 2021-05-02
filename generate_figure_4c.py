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
# Plot G3 - Cyclic exponent, three eigenphases
# ----------------------------------------------------------------------

def displaySwitchPoints(instance) :
   for i in range(instance.phiMu.shape[1]) :
      v = np.copy(instance.phiMu[:,i])
      plt.plot(instance.phiIter, v, '-',color=('C%d' % i))

      # Plot reference phase values
      for x in instance.phiStar :
         plt.plot([1,1e6],[x,x],'k--')

      # Plot switch points
      for (idx,iter) in zip(instance.switchIdx, instance.switchIter) :
         idx2 = np.argwhere(instance.phiIter <= iter)
         idx2 = idx2[-1]
         plt.plot([iter,iter],[instance.phiMu[idx2,idx]-0.5, instance.phiMu[idx2,idx]+0.5],'k--')

   plt.xscale('log')
   plt.show()
    
def getSwitchPointBias(instances) :
   dist = np.zeros((len(instances),instances[0].phiMu.shape[1]),dtype=np.double)
   for (i,instance) in enumerate(instances) :
      # Plot switch points
      for (j,(idx,iter)) in enumerate(zip(instance.switchIdx, instance.switchIter)) :
        idx2 = np.argwhere(instance.phiIter <= iter)
        idx2 = idx2[-1]
        d = np.abs(instance.phiMu[idx2,idx] - instance.phiStar)
        idx2 = (d > np.pi)
        d[idx2] = 2*np.pi - d[idx2]
        dist[i,j] =np.min(d)
   return dist


# Show mean deviations at the transition point
#kv  = [ 1,  2,  3,  5, 10, 20,  5,  5]
#idx = [41, 43, 45, 47, 48, 49, 50, 51]
#for (idx,k) in zip(idx,kv) :
#   instances = loadInstances('./cache/experiment_hybrid_%d.dat' % idx)
#   dist = getSwitchPointBias(instances)
#   print("%2d %s" % (k, str(np.median(dist,axis=0))))

instancesA = loadInstances('./cache/experiment_hybrid_47.dat')
instancesB = loadInstances('./cache/experiment_hybrid_48.dat')
instancesC = loadInstances('./cache/experiment_hybrid_49.dat')
instancesD = loadInstances('./cache/experiment_hybrid_51.dat')

for i in range(len(instancesA)) : instancesA[i].collateMaxWeight(matchIter=1000)
for i in range(len(instancesB)) : instancesB[i].collateMaxWeight(matchIter=1000)
for i in range(len(instancesC)) : instancesC[i].collateMaxWeight(matchIter=1000)
for i in range(len(instancesD)) : instancesD[i].collateMaxWeight(matchIter=1000)


# Display the averate switch iterations along with standard deviation
for instances in [instancesB,instancesC,instancesD] :
   data  = [problem.switchIter for problem in instances]
   width = max([len(s) for s in data])
   switches = np.zeros((len(data),width),dtype=int)
   for (i,d) in enumerate(data) :
      switches[i,:len(d)] = d

   print("Mean switch: %s (%s)" % (", ".join(["%d"%f for f in np.mean(switches,axis=0)]),
                                   ", ".join(["%.2f"%f for f in np.std(switches,axis=0)])))

h1 = plotCollatedMu(instancesD,'-',alpha=0.40,color='C1')
h2 = plotCollatedMu(instancesA,'-',alpha=0.80,color='C1')
h3 = plotCollatedMu(instancesB,'-',alpha=0.40,color='C2')
h4 = plotCollatedMu(instancesC,'-',alpha=1.00,color='C0')


plt.text(1.147e6,1.12e0,'A',fontweight='medium', fontsize=fontsize)
plt.text(1.147e6,3.4e-4,'B',fontweight='medium', fontsize=fontsize)
plt.text(1.147e6,1.7e-4,'C',fontweight='medium', fontsize=fontsize)
plt.text(1.147e6,2.7e-3,'D',fontweight='medium', fontsize=fontsize)
plt.legend(['A','B','C','D'],
           ['Mixed ($c_{\max}$ = 5)', 'Mixed ($c_{\max}$ = 10)','Mixed ($c_{\max}$ = 20)','Mixed ($c_{\max}$ = 5, $n_{\max}$=5,000)'],
           handler_map={str: StringObjectHandler()}, loc='lower left', fontsize=12)



plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Median phase error', fontsize=fontsize)
plt.ylim([3e-5,4e0])
setTickFontsize(fontsize)
exportFigure("FigExperimentHybrid_G3", label='4c')
