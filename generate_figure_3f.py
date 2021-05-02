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
instances8 = loadInstances('./cache/experiment_hybrid_8.dat')
instances9 = loadInstances('./cache/experiment_hybrid_9.dat')

# Collate the results
for i in range(len(instances8)) : instances8[i].collate(matchIter=2e5)
for i in range(len(instances9)) : instances9[i].collate(matchIter=2e5)


# Normal only
instancesA = loadInstances('cache/experiment_hybrid_20.dat') # Curvilinear,  2
instancesB = loadInstances('cache/experiment_hybrid_21.dat') # Curvilinear,  64
instancesC = loadInstances('cache/experiment_hybrid_22.dat') # Curvilinear,  512 (same as 8)
instancesD = loadInstances('cache/experiment_hybrid_23.dat') # Backtracking, 512
instancesE = loadInstances('cache/experiment_hybrid_24.dat') # Curvilinear,  100,000

for i in range(len(instancesA))  : instancesA[i].collate(matchIter=2e5)
for i in range(len(instancesB))  : instancesB[i].collate(matchIter=2e5)
for i in range(len(instancesC))  : instancesC[i].collate(matchIter=2e5)
for i in range(len(instancesD))  : instancesD[i].collate(matchIter=9e5)
for i in range(len(instancesE))  : instancesE[i].collate(matchIter=2e5)

# Display the average runtimes along with standard deviation
for instances in [instancesA, instancesB, instancesC, instancesD, instancesE] :
   runtime = np.asarray([problem.runtime for problem in instances])
   print("Mean runtime (%s): %.2f (%.2f)" % (instances[0].filename, np.mean(runtime),np.std(runtime)))

plotCollatedWeights(instancesA,'--',alpha=0.70,style='median',color='C1')
plotCollatedWeights(instancesB,'-',alpha=0.30,style='median',color='C2')
plotCollatedWeights(instancesD,'-',alpha=0.50,style='median',color='C2')
plotCollatedWeights(instancesE,'-',alpha=0.70,style='median',color='C2')


# Mixed
instancesA = loadInstances('cache/experiment_hybrid_25.dat') # Curvilinear,  2
instancesB = loadInstances('cache/experiment_hybrid_26.dat') # Curvilinear,  64
instancesC = loadInstances('cache/experiment_hybrid_27.dat') # Curvilinear,  512 (same as 9)
instancesD = loadInstances('cache/experiment_hybrid_28.dat') # Backtracking, 512
instancesE = loadInstances('cache/experiment_hybrid_29.dat') # Curvilinear,  100,000

for i in range(len(instancesA))  : instancesA[i].collate(matchIter=2e5)
for i in range(len(instancesB))  : instancesB[i].collate(matchIter=2e5)
for i in range(len(instancesC))  : instancesC[i].collate(matchIter=2e5)
for i in range(len(instancesD))  : instancesD[i].collate(matchIter=2e5)
for i in range(len(instancesE))  : instancesE[i].collate(matchIter=2e5)

# Display the average runtimes along with standard deviation
for instances in [instancesA, instancesB, instancesC, instancesD, instancesE] :
   runtime = np.asarray([problem.runtime for problem in instances])
   print("Mean runtime (%s): %.2f (%.2f)" % (instances[0].filename, np.mean(runtime),np.std(runtime)))

plotCollatedWeights(instancesA,'--',alpha=0.90,style='median',color='C0')
plotCollatedWeights(instancesB,'-',alpha=0.30,style='median',color='C0')
plotCollatedWeights(instancesD,'-',alpha=0.50,style='median',color='C0')
plotCollatedWeights(instancesE,'-',alpha=0.70,style='median',color='C0')

plotCollatedWeights(instances8,'-',alpha=1.0,color='C2',style='median')
plotCollatedWeights(instances9,'-',alpha=1.0,color='C0',style='median')

plt.text(3.7e2,6.5e-3,'A',fontweight='medium', fontsize=fontsize)
plt.text(2.5e3,3e-2,'B',fontweight='medium', fontsize=fontsize)
plt.text(1.3e6,2.7e-2,'C',fontweight='medium', fontsize=fontsize)
plt.legend(['A','B','C'], ['Mixed approach','Normal prior only','Normal prior only (T=2)'],
           handler_map={str: StringObjectHandler()}, loc='lower left', fontsize=fontsize)
plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Median absolute weight error', fontsize=fontsize)
plt.xlim([10,2.2e6])
setTickFontsize(fontsize)
exportFigure("FigExperimentHybrid_D3", label='3f')
