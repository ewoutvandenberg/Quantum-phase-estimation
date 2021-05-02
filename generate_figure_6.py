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

from generic import *

# Figure 6a
data = loadInstances('./cache/experiment_hybrid_314.dat')
plotInstanceFilter(data[15],3)
plt.xlim([1,1e6])
plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Estimated phase', fontsize=fontsize)
setTickFontsize(fontsize)
plt.ylim([0,6.5])
plt.xlim([1,1e6])
exportFigure("FigExperimentHybrid_SpuriousA", label='6a')


# Figure 6b
filterPhaseData(data[15],3,True)
plt.xlim([1,1e6])
plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Estimated phase', fontsize=fontsize)
setTickFontsize(fontsize)
plt.ylim([0,6.5])
plt.xlim([1,1e6])
exportFigure("FigExperimentHybrid_SpuriousB", label='6b')


# Figure 6c
data = loadInstances('./cache/experiment_hybrid_2%02d.dat' % (9+3*20))
plotInstanceFilter(data[1], 9)
plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Estimated phase', fontsize=fontsize)
setTickFontsize(fontsize)
plt.ylim([0,6.3])
plt.xlim([1,1e6])
exportFigure("FigExperimentHybrid_Multiple", label='6c')
