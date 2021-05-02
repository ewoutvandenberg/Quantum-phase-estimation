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
# G5 - Adaptive schemes - Three eigenstates
# ----------------------------------------------------------------------

"""
# Normal only - does not work
#instances80 = loadInstances('./cache/experiment_hybrid_80.dat')
#instances82 = loadInstances('./cache/experiment_hybrid_82.dat')
#instances84 = loadInstances('./cache/experiment_hybrid_84.dat')
#instances86 = loadInstances('./cache/experiment_hybrid_86.dat')

#for i in range(len(instances80)) : instances80[i].collateMaxWeight(matchIter=1000)
#for i in range(len(instances82)) : instances82[i].collateMaxWeight(matchIter=1000)
#for i in range(len(instances84)) : instances84[i].collateMaxWeight(matchIter=1000)
#for i in range(len(instances86)) : instances86[i].collateMaxWeight(matchIter=1000)

#plotCollatedMu(instances80,'-',alpha=1.0,color='C1') # Does not converge
#plotCollatedMu(instances82,'-',alpha=1.0,color='C1') # Does not converge
#plotCollatedMu(instances84,'-',alpha=1.0,color='C1') # Does not converge
#plotCollatedMu(instances86,'-',alpha=1.0,color='C1') # Does not converge

# Non-adaptive (k=1) - does not work
#instances81 = loadInstances('./cache/experiment_hybrid_81.dat')
#for i in range(len(instances81)) : instances81[i].collateMaxWeight(matchIter=1000)
#plotCollatedMu(instances81,'-',alpha=1.0,color='C1') # Does not converge

# Pure adaptive - does not work
#instances83 = loadInstances('./cache/experiment_hybrid_83.dat')
#for i in range(len(instances83)) : instances83[i].collateMaxWeight(matchIter=1000)
#plotCollatedMu(instances83,'-',alpha=1.0,color='C0') # Does not converge
"""

# Adaptive with cycling, nEig = 5
instancesA = loadInstances('./cache/experiment_hybrid_92.dat')
plotCollatedMu(instancesA,'-',alpha=0.3,color='C1') # Converges

# Cyclic (cmax=50), nEig = 5
instancesA = loadInstances('./cache/experiment_hybrid_93.dat')
plotCollatedMu(instancesA,'-',alpha=0.3,color='C0') # Converges

# Adaptive with cycling, nEig = 3
instances85 = loadInstances('./cache/experiment_hybrid_85.dat')
plotCollatedMu(instances85,'-',alpha=0.6,color='C1') # Converges

# Cyclic (cmax=50), nEig = 3
instances87 = loadInstances('./cache/experiment_hybrid_87.dat')
plotCollatedMu(instances87,'-',alpha=1.0,color='C0') # Converges


# Plot results
plt.text(6.5e1,7.8e-3,'A',fontweight='medium', fontsize=fontsize)
plt.text(1.6e1,1.2e-1,'B',fontweight='medium', fontsize=fontsize)
plt.legend(['A','B'],
           ['Cyclic ($c_{\max}$ = 50)',
            'Adaptive cyclic ($c_{\max}$ = 50)'],
           handler_map={str: StringObjectHandler()}, loc='lower left', fontsize=12)


plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Median phase error', fontsize=fontsize)
setTickFontsize(fontsize)
plt.ylim([5.2e-5,3])
exportFigure("FigExperimentHybrid_G5", label='4e')
