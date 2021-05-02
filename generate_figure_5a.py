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
# N1 - Decoherence - Single eigenstate
# ----------------------------------------------------------------------

instancesA = loadInstances('./cache/experiment_hybrid_102.dat')
plotCollatedMu(instancesA,'-',alpha=0.3,color='C0')
instancesA = loadInstances('./cache/experiment_hybrid_103.dat')
plotCollatedMu(instancesA,'-',alpha=1.0,color='C0')

instancesA = loadInstances('./cache/experiment_hybrid_104.dat')
plotCollatedMu(instancesA,'-',alpha=0.3,color='C0')
instancesA = loadInstances('./cache/experiment_hybrid_105.dat')
plotCollatedMu(instancesA,'-',alpha=1.0,color='C0')

instancesA = loadInstances('./cache/experiment_hybrid_100.dat') # Fails
plotCollatedMu(instancesA,'-',alpha=0.5,color='C1')
instancesA = loadInstances('./cache/experiment_hybrid_101.dat')
plotCollatedMu(instancesA,'-',alpha=1.0,color='C0')

instancesA = loadInstances('./cache/experiment_hybrid_106.dat')
plotCollatedMu(instancesA,'-',alpha=0.3,color='C2')
instancesA = loadInstances('./cache/experiment_hybrid_107.dat')
plotCollatedMu(instancesA,'-',alpha=1.0,color='C2')

plt.text(1.147e6,5.0e-1,'A',fontweight='medium', fontsize=fontsize)
plt.text(1.147e6,8.3e-4,'B',fontweight='medium', fontsize=fontsize)
plt.text(1.147e6,2.8e-5,'C',fontweight='medium', fontsize=fontsize)

plt.legend(['A','B','C'],
           ['Cyclic ($c_{\max}$ = 50), normal-only',
            'Fixed ($k = 1$)',
            'Adaptive/cyclic ($c_{\max}$ = 50)'],
           handler_map={str: StringObjectHandler()}, loc='lower left', fontsize=12)


plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Median phase error', fontsize=fontsize)
setTickFontsize(fontsize)
plt.ylim([4e-6,3])

exportFigure("FigExperimentHybrid_N1", label='5a')
