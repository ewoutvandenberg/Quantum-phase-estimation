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
# N2 - Decoherence - Three eigenstates
# ----------------------------------------------------------------------

#instancesA = loadInstances('./cache/experiment_hybrid_87.dat')
#for i in range(len(instancesA)) : instancesA[i].collateMaxWeight(matchIter=1000)
#plotCollatedMu(instancesA,'-',alpha=1,color='C1')

# Purely adaptive -- fail
#instancesA = loadInstances('./cache/experiment_hybrid_113.dat')
#instancesA = loadInstances('./cache/experiment_hybrid_115.dat')
#instancesA = loadInstances('./cache/experiment_hybrid_117.dat')

# Adaptive cyclic -- succeed
#instancesA = loadInstances('./cache/experiment_hybrid_114.dat')
#plotCollatedMu(instancesA,'-',alpha=0.4,color=(1,0,0))
#instancesA = loadInstances('./cache/experiment_hybrid_116.dat')
#plotCollatedMu(instancesA,'-',alpha=0.7,color=(1,0,0))
#instancesA = loadInstances('./cache/experiment_hybrid_118.dat')
#plotCollatedMu(instancesA,'-',alpha=1.0,color=(1,0,0))

# Cyclic with cmax values 50, 80, and 100:
#instancesA = loadInstances('./cache/experiment_hybrid_119.dat')
#instancesA = loadInstances('./cache/experiment_hybrid_120.dat')
#instancesA = loadInstances('./cache/experiment_hybrid_121.dat')
#instancesA = loadInstances('./cache/experiment_hybrid_122.dat') # Same as #120 but with kMax=1000
#instancesA = loadInstances('./cache/experiment_hybrid_123.dat') # Same as #121 but with kMax=1000

# Cyclic with cmax=10
#instancesA = loadInstances('./cache/experiment_hybrid_110.dat')
#plotCollatedMu(instancesA,'-',alpha=0.3,color='C0')



# Cyclic with cmax values 20 and 50
instancesA = loadInstances('./cache/experiment_hybrid_111.dat')
h1=plotCollatedMu(instancesA,'-',alpha=0.3,color='C0')
instancesA = loadInstances('./cache/experiment_hybrid_112.dat')
h2=plotCollatedMu(instancesA,'-',alpha=0.6,color='C1')

# Cyclic with cmax=100, and kmax values 200 and 1000
instancesA = loadInstances('./cache/experiment_hybrid_121.dat')
h3=plotCollatedMu(instancesA,'--',alpha=1.0,color='C2')
instancesA = loadInstances('./cache/experiment_hybrid_123.dat')
h4=plotCollatedMu(instancesA,'-',alpha=1.0,color='C2')


plt.legend([h1[0],h2[0],h3[0],h4[0]],
           ['Cyclic, $c_{\max}$ = 20',
            'Cyclic, $c_{\max}$ = 50',
            'Cyclic, $c_{\max}$ = 100, $n_{\max}$ = 200',
            'Cyclic, $c_{\max}$ = 100, $n_{\max}$ = 1000'],loc='lower left', fontsize=12)

plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Median phase error', fontsize=fontsize)
plt.xlim([1e1,1e6])
plt.ylim([3e-6,4])
setTickFontsize(fontsize)
exportFigure("FigExperimentHybrid_N2", label='5b')
