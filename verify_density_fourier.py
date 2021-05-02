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

import numpy as np
import bpe

phi = np.linspace(0,2*np.pi,10000,endpoint=False)

# Function for printing
def isSmall(value) :
   if (np.abs(value) < 1e-14) :
      return "OK"
   else :
      return "Failed (%s)" % (np.abs(value))

# Initialization to zero
density = bpe.DensityFourier()
density.initZero()
p1 = 0 * phi
p2 = density.evaluate(phi)
print("Initialization to zero           : %s" % isSmall(np.max(np.abs(p1-p2))))

# Initialization to one
density.initOne()
p1 = np.ones(phi.size,dtype=np.double)
p2 = density.evaluate(phi)
print("Initialization to one            : %s" % isSmall(np.max(np.abs(p1-p2))))

# Initialization to uniform distribution
density.initUniform()
p1 = np.ones(phi.size,dtype=np.double) / (2*np.pi)
p2 = density.evaluate(phi)
print("Initialization to uniform distr. : %s" % isSmall(np.max(np.abs(p1-p2))))

# Initialization to cosine
kVals  = [  0,   2,  3]
bVals  = [0.1, 0.0,2.4]
wVals  = [0.4,-0.1,1.5]
density.initCosine(kVals,bVals,wVals)
p1 = np.zeros(phi.size,dtype=np.double)
for (k,beta,weight) in zip(kVals,bVals,wVals) :
   p1 += weight * np.cos(k*phi + beta)
p2 = density.evaluate(phi)
print("Initialization to cosine         : %s" % isSmall(np.max(np.abs(p1-p2))))

# Initialization to sine
kVals  = [  0,   1,  4]
bVals  = [0.3, 1.0,1.2]
wVals  = [0.8, 0.4,0.5]
density.initSine(kVals,bVals,wVals)
p1 = np.zeros(phi.size,dtype=np.double)
for (k,beta,weight) in zip(kVals,bVals,wVals) :
   p1 += weight * np.sin(k*phi + beta)
p2 = density.evaluate(phi)
print("Initialization to sine           : %s" % isSmall(np.max(np.abs(p1-p2))))

# Initialization based on probability - single term
kVals = [  4]
bVals = [0.5]
density.initProbability(kVals, bVals)
p1 = np.ones(phi.size,dtype=np.double)
for (k,beta) in zip(kVals,bVals) :
   p1 *= (1. + np.cos(k*phi + beta)) / 2.
p2 = density.evaluate(phi)
print("Initialization to probability (1): %s" % isSmall(np.max(np.abs(p1-p2))))

# Initialization based on probability - multiple terms
kVals = [  4,    2,    1]
bVals = [0.5, -0.2, -1.2]
density.initProbability(kVals, bVals)
p1 = np.ones(phi.size,dtype=np.double)
for (k,beta) in zip(kVals,bVals) :
   p1 *= (1. + np.cos(k*phi + beta)) / 2.
p2 = density.evaluate(phi)
print("Initialization to probability (2): %s" % isSmall(np.max(np.abs(p1-p2))))

# Initialization based on normal distribution - with warp-around
mu = 0.4
sigma = 0.3
density.initNormal(mu,sigma,200)
p1 = np.zeros(phi.size,dtype=np.double)
for offset in range(-6,7) :
   p1 += np.exp(-(((phi-mu + offset*2*np.pi)/sigma)**2) / 2) / (sigma * np.sqrt(2*np.pi))
p2 = density.evaluate(phi)
print("Initialization to normal         : %s" % isSmall(np.max(np.abs(p1-p2))))


# Multiplication of distributions
mu1 = 0.4; sigma1 = 0.6
mu2 = 2.3; sigma2 = 0.9

density1 = bpe.DensityFourier()
density2 = bpe.DensityFourier()
density1.initNormal(mu1,sigma1,100)
density2.initNormal(mu2,sigma2,100)
density = density1 * density2
p1a = np.zeros(phi.size,dtype=np.double)
p1b = np.zeros(phi.size,dtype=np.double)
for offset in range(-6,7) :
   p1a += np.exp(-(((phi-mu1 + offset*2*np.pi)/sigma1)**2) / 2) / (sigma1 * np.sqrt(2*np.pi))
   p1b += np.exp(-(((phi-mu2 + offset*2*np.pi)/sigma2)**2) / 2) / (sigma2 * np.sqrt(2*np.pi))
p1 = p1a * p1b
p2 = density.evaluate(phi)
print("Multiplication of distributions  : %s" % isSmall(np.max(np.abs(p1-p2))))

density = density1 + density2
p1 = p1a + p1b
p2 = density.evaluate(phi)
print("Addition of distributions        : %s" % isSmall(np.max(np.abs(p1-p2))))

density = density1 - density2
p1 = p1a - p1b
p2 = density.evaluate(phi)
print("Subtraction of distributions     : %s" % isSmall(np.max(np.abs(p1-p2))))

density = -3. * density1 + density2 * 2.
p1 = -3 * p1a + p1b * 2.
p2 = density.evaluate(phi)
print("Scaling of distributions         : %s" % isSmall(np.max(np.abs(p1-p2))))


#import matplotlib.pyplot as plt

