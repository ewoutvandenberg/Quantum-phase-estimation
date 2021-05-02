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
# G4 - Adaptive schemes - Single eigenstate
# ----------------------------------------------------------------------

# Non-adaptive (k=1)
instances60 = loadInstances('./cache/experiment_hybrid_60.dat')
instances61 = loadInstances('./cache/experiment_hybrid_61.dat')
plotCollatedMu(instances60,'-',alpha=0.5,color='C2') # Normal only - converges
plotCollatedMu(instances61,'-',alpha=1.0,color='C2') # Mixed - converges

# Pure adaptive
instances62 = loadInstances('./cache/experiment_hybrid_62.dat')
instances63 = loadInstances('./cache/experiment_hybrid_63.dat')
plotCollatedMu(instances62,'-',alpha=0.5,color='C1') # Normal only - converges
plotCollatedMu(instances63,'-',alpha=1.0,color='C1') # Mixed - converges

# Adaptive with cycling
instances64 = loadInstances('./cache/experiment_hybrid_64.dat')
instances65 = loadInstances('./cache/experiment_hybrid_65.dat')
plotCollatedMu(instances64,'-',alpha=0.5,color='C1') # Normal only - converges
plotCollatedMu(instances65,'-',alpha=1.0,color='C1') # Mixed - converges

# Cyclic (cmax=50)
#instances66 = loadInstances('./cache/experiment_hybrid_66.dat')
#plotCollatedMu(instances66,'-',alpha=0.5,color='C1') # Normal only - fails
instances67 = loadInstances('./cache/experiment_hybrid_67.dat')
plotCollatedMu(instances67,'-',alpha=1.0,color='C1') # Mixed - converges


# Pure adaptive (kmax=4096)
instances68 = loadInstances('./cache/experiment_hybrid_68.dat')
instances69 = loadInstances('./cache/experiment_hybrid_69.dat')
plotCollatedMu(instances68,'-',alpha=0.5,color='C3') # Normal only - converges
plotCollatedMu(instances69,'-',alpha=1.0,color='C3') # Mixed - converges

# Adaptive with cycling (kmax=4096)
instances70 = loadInstances('./cache/experiment_hybrid_70.dat')
instances71 = loadInstances('./cache/experiment_hybrid_71.dat')
plotCollatedMu(instances70,'-',alpha=0.5,color='C9') # Normal only - converges
plotCollatedMu(instances71,'-',alpha=1.0,color='C9') # Mixed - converges

# Cyclic (kmax=4096)
#instances72 = loadInstances('./cache/experiment_hybrid_72.dat')
#plotCollatedMu(instances72,':',alpha=0.5,color='C11') # Normal only - fails
instances73 = loadInstances('./cache/experiment_hybrid_73.dat')
plotCollatedMu(instances73,'-',alpha=1.0,color='C4') # Mixed - converges


plt.text(1.147e6,6.0e-4,'A',fontweight='medium', fontsize=fontsize)
plt.text(1.147e6,1.4e-5,'B',fontweight='medium', fontsize=fontsize)
plt.text(1.147e6,2.4e-6,'C',fontweight='medium', fontsize=fontsize)
plt.text(1.147e6,4.5e-7,'D',fontweight='medium', fontsize=fontsize)
plt.legend(['A','B','C','D'],
           ['Fixed (k = 1)',
            'Adaptive/cyclic ($k_{\max}$ = 50)',
            'Adaptive/pure cyclic ($k_{\max}$ = 4096)',
            'Pure adaptive ($k_{\max}$ = 4096)'],
           handler_map={str: StringObjectHandler()}, loc='upper right', fontsize=10)


plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Median phase error', fontsize=fontsize)
setTickFontsize(fontsize)
plt.ylim([2e-7,4e1])
exportFigure("FigExperimentHybrid_G4", label='4d')
