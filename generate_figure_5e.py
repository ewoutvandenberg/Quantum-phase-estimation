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
# Plot I2 - Fourier-based sampling, three eigenphases
# ----------------------------------------------------------------------

# Curve A
instancesA = loadInstances('./cache/experiment_hybrid_160.dat') # Noiseless, normal, succeeds
instancesB = loadInstances('./cache/experiment_hybrid_161.dat') # Noiseless, mixed, succeeds
plotCollatedMu(instancesA,'-',alpha=0.3,color='C1')
plotCollatedMu(instancesB,'-',alpha=1.0,color='C1')

#instancesA = loadInstances('./cache/experiment_hybrid_162.dat') # 2% noise, normal, succeeds
#instancesB = loadInstances('./cache/experiment_hybrid_163.dat') # 2% noise, mixed, succeeds
#plotCollatedMu(instancesA,'-',alpha=0.3,color='C1')
#plotCollatedMu(instancesB,'-',alpha=1.0,color='C1')

#instancesA = loadInstances('./cache/experiment_hybrid_164.dat') # 5% noise, normal, succeeds
#instancesB = loadInstances('./cache/experiment_hybrid_165.dat') # 5% noise, mixed, succeeds
#plotCollatedMu(instancesA,'-',alpha=0.3,color='C1')
#plotCollatedMu(instancesB,'-',alpha=1.0,color='C1')

# Curve B
instancesA = loadInstances('./cache/experiment_hybrid_166.dat') # 8% noise, normal, succeeds
#instancesB = loadInstances('./cache/experiment_hybrid_167.dat') # 8% noise, mixed, partially fails
instancesC = loadInstances('./cache/experiment_hybrid_168.dat') # 8% noise, mixed, succeeds (kMax=2000)
plotCollatedMu(instancesA,'-',alpha=0.4,color='C0')
#plotCollatedMu(instancesB,'-',alpha=0.7,color='C0')
plotCollatedMu(instancesC,'-',alpha=1.0,color='C0')

# Normal only starts to deteriorate, need to use larger kMax
#instancesA = loadInstances('./cache/experiment_hybrid_169.dat') # 10% noise, normal, partially fails
#instancesC = loadInstances('./cache/experiment_hybrid_170.dat') # 10% noise, mixed, succeeds (kMax=2000)
#plotCollatedMu(instancesA,'-',alpha=0.4,color='C9')
#plotCollatedMu(instancesC,'-',alpha=1.0,color='C9')

#instancesA = loadInstances('./cache/experiment_hybrid_171.dat') # 15% noise, normal, fails
#instancesB = loadInstances('./cache/experiment_hybrid_172.dat') # 15% noise, mixed, partially fails
#instancesC = loadInstances('./cache/experiment_hybrid_173.dat') # 15% noise, mixed, succeeds (kMax=2000)
#plotCollatedMu(instancesA,'-',alpha=0.4,color='C3')
#plotCollatedMu(instancesC,'-',alpha=1.0,color='C3')

#instancesA = loadInstances('./cache/experiment_hybrid_174.dat') # 20% noise, normal, fails
#instancesB = loadInstances('./cache/experiment_hybrid_175.dat') # 20% noise, mixed, partially fails
#instancesC = loadInstances('./cache/experiment_hybrid_176.dat') # 20% noise, mixed, succeeds (kMax=2000)
#plotCollatedMu(instancesA,'-',alpha=0.4,color='C3')
#plotCollatedMu(instancesB,':',alpha=0.4,color='C3')
#plotCollatedMu(instancesC,'-',alpha=1.0,color='C3')

#instancesC = loadInstances('./cache/experiment_hybrid_177.dat') # 25% noise, mixed, succeeds (kMax=1000)
#plotCollatedMu(instancesC,'-',alpha=1.0,color='C5')

instancesC1 = loadInstances('./cache/experiment_hybrid_178.dat') # 30% noise, mixed, partially succeeds (kMax=1000)
plotCollatedMu(instancesC1,'--',alpha=0.6,color='C2')
instancesC2 = loadInstances('./cache/experiment_hybrid_179.dat') # 30% noise, mixed, succeeds (kMax=2000)
plotCollatedMu(instancesC2,'-',alpha=1.0,color='C2')


# AT WHICH ITERATION DO THE DISTRIBUTIONS SWITCH?
#instancesC = loadInstances('./cache/experiment_hybrid_180.dat') # 35% noise, mixed, partially succeeds (kMax=1000)
#plotCollatedMu(instancesC,'-',alpha=0.5,color='C7')
#instancesC = loadInstances('./cache/experiment_hybrid_181.dat') # 30% noise, mixed, partially succeeds (kMax=2000)
#plotCollatedMu(instancesC,'-',alpha=1.0,color='C7')


for instances in [instancesC1,instancesC2] :
   print("--------- Switch iteration ---------")
   for i in range(3) :
      s = []
      for problem in instances :
         if (len(problem.switchIter) > i) :
            s.append(problem.switchIter[i])
      print(np.median(s))
   print("")


# Bracket
x1 = 1.05e6; x2 = 1.12e6
y1 = 2.1e-3; y2 = 6.5e-4
plt.plot([x1,x2,x2,x1],[y1,y1,y2,y2],'k-',linewidth=0.75)

plt.text(1.195e6,6.0e-5,'A',fontweight='medium', fontsize=fontsize)
plt.text(1.195e6,1.1e-4,'B',fontweight='medium', fontsize=fontsize)
plt.text(1.220e6,9.2e-4,'D',fontweight='medium', fontsize=fontsize)
plt.text(1.195e6,5.5e-2,'C',fontweight='medium', fontsize=fontsize)

plt.legend(['A','B','C','D'],
           ['Noiseless (normal / mixed)',
            '8% noise (normal / mixed)',
            '30% noise, mixed ($n_{\max} = 1000$)',
            '30% noise, mixed ($n_{\max} = 2000$)'],
           handler_map={str: StringObjectHandler()}, loc='upper right', fontsize=11)

xl = plt.xlim()
plt.xlim([1e1,1.9e6])
plt.ylim([3.5e-5, 9])
plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Median phase error', fontsize=fontsize)
setTickFontsize(fontsize)
exportFigure("FigExperimentHybrid_I2", label='5e')
