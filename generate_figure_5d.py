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
# Plot I1 - Fourier-based sampling, single phase (phi = 5.0)
# ----------------------------------------------------------------------

instancesA = loadInstances('./cache/experiment_hybrid_150.dat') # 26%
#instancesB = loadInstances('./cache/experiment_hybrid_151.dat')
plotCollatedMu(instancesA,'-',alpha=0.4,color='C1')
#plotCollatedMu(instancesB,'-',alpha=1.0,color='C1')

instancesA = loadInstances('./cache/experiment_hybrid_140.dat')
instancesB = loadInstances('./cache/experiment_hybrid_141.dat')
plotCollatedMu(instancesA,'-',alpha=0.3,color='C0') # Normal only
plotCollatedMu(instancesB,'-',alpha=1.0,color='C0') # Mixed

#instancesA = loadInstances('./cache/experiment_hybrid_142.dat') # 2% read-out error
#instancesB = loadInstances('./cache/experiment_hybrid_143.dat')
#plotCollatedMu(instancesA,'-',alpha=0.3,color='C3')
#plotCollatedMu(instancesB,'-',alpha=1.0,color='C3')

#instancesA = loadInstances('./cache/experiment_hybrid_144.dat') # 5%
#instancesB = loadInstances('./cache/experiment_hybrid_145.dat')
#plotCollatedMu(instancesA,'-',alpha=0.3,color='C3')
#plotCollatedMu(instancesB,'-',alpha=1.0,color='C3')

instancesA = loadInstances('./cache/experiment_hybrid_146.dat') # 20%
instancesB = loadInstances('./cache/experiment_hybrid_147.dat')
plotCollatedMu(instancesA,'-',alpha=0.3,color='C2') 
plotCollatedMu(instancesB,'-',alpha=1.0,color='C2')

#instancesA = loadInstances('./cache/experiment_hybrid_148.dat') # 25%
#instancesB = loadInstances('./cache/experiment_hybrid_149.dat')
#plotCollatedMu(instancesA,'-',alpha=0.3,color='C9') # Succeeds
#plotCollatedMu(instancesB,'-',alpha=1.0,color='C9') # Succeeds

#instancesA = loadInstances('./cache/experiment_hybrid_152.dat') # 30%
#instancesB = loadInstances('./cache/experiment_hybrid_153.dat')
#plotCollatedMu(instancesA,'-',alpha=0.3,color='C8') # Fails
#plotCollatedMu(instancesB,'-',alpha=1.0,color='C8')

#instancesA = loadInstances('./cache/experiment_hybrid_154.dat') # 40%
instancesB = loadInstances('./cache/experiment_hybrid_155.dat')
#plotCollatedMu(instancesA,':',alpha=0.3,color='C0') # Fails 
plotCollatedMu(instancesB,'-',alpha=1.0,color='C8')

#instancesA = loadInstances('./cache/experiment_hybrid_156.dat') # 45%
instancesB = loadInstances('./cache/experiment_hybrid_157.dat')
#plotCollatedMu(instancesA,':',alpha=0.3,color='C9') # Fails
plotCollatedMu(instancesB,'-',alpha=1.0,color='C4')

#instancesA = loadInstances('./cache/experiment_hybrid_158.dat') # 48%
instancesB = loadInstances('./cache/experiment_hybrid_159.dat')
#plotCollatedMu(instancesA,':',alpha=0.3,color='C9') # Fails
plotCollatedMu(instancesB,'-',alpha=1.0,color='C9')

plt.text(2.0e1,3.0e-3,'A',fontweight='medium', fontsize=fontsize)
plt.text(1.1e2,1.13e-2,'B',fontweight='medium', fontsize=fontsize)
plt.text(6.1e2,2.2e-2,'C',fontweight='medium', fontsize=fontsize)
plt.text(2.5e3,3.0e-2,'D',fontweight='medium', fontsize=fontsize)
plt.text(1.4e4,3.8e-2,'E',fontweight='medium', fontsize=fontsize)
plt.text(4.8e5,6.9e-2,'F',fontweight='medium', fontsize=fontsize)

plt.legend(['A','B','C','D','E','F'],
           ['Noiseless',
            '20% noise',
            '40% noise',
            '45% noise',
            '48% noise',
            '26% noise (normal-only)'],
           handler_map={str: StringObjectHandler()}, loc='lower left', ncol=2, fontsize=12)
           
plt.ylim([1e-6,5e0])
plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Median phase error', fontsize=fontsize)
setTickFontsize(fontsize)
exportFigure("FigExperimentHybrid_I1", label='5d')
