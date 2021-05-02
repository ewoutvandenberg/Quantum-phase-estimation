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
# Plot D - Estimate using five eigenvalues
# ----------------------------------------------------------------------

# Load the results
instances8 = loadInstances('./cache/experiment_hybrid_8.dat')
instances9 = loadInstances('./cache/experiment_hybrid_9.dat')

# Collate the results
for i in range(len(instances8)) : instances8[i].collate(matchIter=2e5)
for i in range(len(instances9)) : instances9[i].collate(matchIter=2e5)

# Display the average runtimes along with standard deviation
for instances in [instances8, instances9] :
   runtime = np.asarray([problem.runtime for problem in instances])
   print("Mean runtime: %.2f (%.2f)" % (np.mean(runtime),np.std(runtime)))

# Display the averate switch iterations along with standard deviation
for instances in [instances8, instances9] :
   data  = [problem.switchIter for problem in instances]
   width = max([len(s) for s in data])
   switches = np.zeros((len(data),width),dtype=int)
   for (i,d) in enumerate(data) :
      switches[i,:len(d)] = d

   print("Mean switch: %s (%s)" % (", ".join(["%d"%f for f in np.mean(switches,axis=0)]),
                                   ", ".join(["%.2f"%f for f in np.std(switches,axis=0)])))


# Generate figure 3d
plotCollatedMu(instances8,'-',alpha=0.4,color='C2',style='median')
plotCollatedMu(instances9,'-',alpha=1.0,color='C0',style='median')

plt.text(7.8e1,1.4e-2,'A',fontweight='medium', fontsize=fontsize)
plt.text(4.5e3,2e-1,'B',fontweight='medium', fontsize=fontsize)
plt.legend(['A','B'], ['Mixed approach','Normal prior only'],
           handler_map={str: StringObjectHandler()}, loc='lower left', fontsize=fontsize)
plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Median phase error', fontsize=fontsize)
setTickFontsize(fontsize)
exportFigure("FigExperimentHybrid_D1", label='3d')


# Generate figure 3e
plotCollatedSigma(instances8,'-',alpha=0.4,color='C2')
plotCollatedSigma(instances9,'-',alpha=1.0,color='C0')
plt.text(6e1,1.7e-2,'A',fontweight='medium', fontsize=fontsize)
plt.text(4.4e3,1.9,'B',fontweight='medium', fontsize=fontsize)
plt.legend(['A','B'], ['Mixed approach','Normal prior only'],
           handler_map={str: StringObjectHandler()}, loc='lower left', fontsize=fontsize)
plt.xlabel('Iteration', fontsize=fontsize)
plt.ylabel('Median standard deviation', fontsize=fontsize)
setTickFontsize(fontsize)
exportFigure("FigExperimentHybrid_D2", label='3e')
