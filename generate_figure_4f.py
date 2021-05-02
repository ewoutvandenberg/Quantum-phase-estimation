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
# G6 - Adaptive schemes - Adaptive k values
# ----------------------------------------------------------------------

# Three eigenvalues (3/5)
instancesA = loadInstances('./cache/experiment_hybrid_92.dat')
plotKVals(instancesA,linestyle='-',color='C0',alpha=0.7,kMaxValue=50)

# Three eigenvalues (3/3)
instancesA = loadInstances('./cache/experiment_hybrid_85.dat')
plotKVals(instancesA,linestyle='-',color='C0',kMaxValue=50)

# Single eigenvalue (1/1)
instancesA = loadInstances('./cache/experiment_hybrid_65.dat')
plotKVals(instancesA,linestyle='-',color='C1',kMaxValue=50)

plt.xscale('log')
plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Median $c_{\max}$', fontsize=fontsize)
setTickFontsize(fontsize)

plt.text(10.2, 15.5,'A',fontweight='medium', fontsize=fontsize)
plt.text(140,17,'B',fontweight='medium', fontsize=fontsize)

plt.legend(['A','B'],
           ['Single eigenstate',
            'Three eigenstates (nEig = 3,5)'],
           handler_map={str: StringObjectHandler()}, loc='upper right', fontsize=12)

plt.ylim([0,63])
exportFigure("FigExperimentHybrid_G6", label='4f')
