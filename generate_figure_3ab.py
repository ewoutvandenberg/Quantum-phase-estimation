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


# Load the results
instancesA = loadInstances('./cache/experiment_hybrid_3.dat')
instancesB = loadInstances('./cache/experiment_hybrid_4.dat')
instancesC = loadInstances('./cache/experiment_hybrid_5.dat')
instancesE = loadInstances('./cache/experiment_hybrid_17.dat')  # Mixed (kMax = 200)
instancesF = loadInstances('./cache/experiment_hybrid_18.dat')  # Mixed (kMax = 1000)
instancesG = loadInstances('./cache/experiment_hybrid_19.dat')  # Mixed (kMax = 5000)

# Collate results
for i in range(len(instancesA)) : instancesA[i].collate(matchIter=1e2)
for i in range(len(instancesB)) : instancesB[i].collate(matchIter=1e2)
for i in range(len(instancesC)) : instancesC[i].collate(matchIter=1e2)
for i in range(len(instancesE)) : instancesE[i].collate(matchIter=1e2)
for i in range(len(instancesF)) : instancesF[i].collate(matchIter=1e2)
for i in range(len(instancesG)) : instancesG[i].collate(matchIter=1e2)

# Display the averate switch iterations along with standard deviation
for instances in [instancesG] :
   data  = [problem.switchIter for problem in instances]
   width = max([len(s) for s in data])
   switches = np.zeros((len(data),width),dtype=int)
   for (i,d) in enumerate(data) :
      switches[i,:len(d)] = d

   print("Mean switch: %s (%s)" % (", ".join(["%d"%f for f in np.mean(switches,axis=0)]),
                                   ", ".join(["%.2f"%f for f in np.std(switches,axis=0)])))


# Display the average runtimes along with standard deviation
for instances in [instancesA, instancesB, instancesC] :
   runtime = np.asarray([problem.runtime for problem in instances])
   print("Median runtime: %.2f (%.2f)" % (np.mean(runtime),np.std(runtime)))

# Generate figure 3a
plotCollatedMu(instancesE,'-',alpha=0.3,color='C4')
plotCollatedMu(instancesA,'-',alpha=0.3)
plotCollatedMu(instancesB,'-',alpha=0.65)
plotCollatedMu(instancesC,'-',alpha=1.0)

plt.text(1.5e3,1e-1,'A',fontweight='medium', fontsize=fontsize)
plt.text(2.6e4,2.5e-2,'B',fontweight='medium', fontsize=fontsize)
plt.text(1.5e5,1.085e-2,'C',fontweight='medium', fontsize=fontsize)
plt.text(6.0e5,6.0e-4,'D',fontweight='medium', fontsize=fontsize)

plt.legend(['A','B','C','D'], ['Max. coef = 200','Max. coef = 1000', 'Max. coef = 5000','Mixed approach (200)'],
           handler_map={str: StringObjectHandler()}, loc='lower left', fontsize=12)
plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Median phase error', fontsize=fontsize)
setTickFontsize(fontsize)
exportFigure("FigExperimentHybrid_B1",label="3a")




# Generate figure 3b
plotCollatedSigma(instancesG,'-',alpha=0.3,color='C4')
plotCollatedSigma(instancesA,'-',alpha=0.3)
plotCollatedSigma(instancesB,'-',alpha=0.65)
plotCollatedSigma(instancesC,'-',alpha=1.0)

epsilon = 1e-4
for kMax in [200,1000,5000] :
   sigmaCritical = bpe.DensityFourier.normalCriticalSigma(None, epsilon, kMax=kMax)
   plt.plot([1,1e6],[sigmaCritical,sigmaCritical],'k-',alpha=0.2)
   plt.text(1.2,1.2*sigmaCritical,'$\sigma(%d)$'%kMax,fontweight='medium', fontsize=11,alpha=0.5)

plt.text(1.5e3,8e-2,'A',fontweight='medium', fontsize=fontsize) # 1.5e3, 1e-1
plt.text(2.3e4,5e-2,'B',fontweight='medium', fontsize=fontsize)
plt.text(1.5e5,3.55e-2,'C',fontweight='medium', fontsize=fontsize) #1.4e5, 1.6e-3
plt.text(6e5,3.2e-4,'D',fontweight='medium', fontsize=fontsize)
plt.legend(['A','B','C','D'], ['Max. coef = 200','Max. coef = 1000', 'Max. coef = 5000','Mixed approach (5000)'],
           handler_map={str: StringObjectHandler()}, loc='upper right', fontsize=12)

plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Median standard deviation', fontsize=fontsize)
plt.xlim([1,1e6])
setTickFontsize(fontsize)
exportFigure("FigExperimentHybrid_B2",label="3b")
