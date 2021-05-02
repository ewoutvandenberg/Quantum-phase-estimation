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
import generic
import bpe

generic.makeExperiment('verify_clib_fourier')


def unittest(index, ref) :
   v = generic.runExperiment('verify_clib_fourier', param=[str(index)], verbose=False, showinfo=False)
   v = v.splitlines()
   v = [s for s in v if s.strip()] # Filter out empty lines

   if (len(v) != len(ref)) :
      print("Test %2d: Failed . . . mismatch in output (%d) and reference (%d)" % (index, len(v), len(ref)))
      return False

   success = True
   for [i,s] in enumerate(v) :
      r1 = np.asarray([float(x) for x in s.split()])

      # Copy reference vector and zero-pad vectors to ensure same length
      r2 = np.copy(np.reshape(ref[i],ref[i].size))
      if (r2.size < r1.size) :
        r2.resize(r1.size,refcheck=False)
      if (r1.size < r2.size) :
        r1.resize(r2.size,refcheck=False)

      # Compute the vectors
      if (np.max(np.abs(r1-r2)) > 1e-7) :
         print("         Line %d" % (i+1))
         print("         %s" % str(r1))
         print("         %s" % str(r2))
         success = False

   print("Test %2d: %s" % (index, "Success" if success else "Failed"))

   return success



dist = bpe.DensityFourier()


dist.initZero()
unittest(0,[dist.coefCos,dist.coefSin])

dist.initOne()
unittest(1,[dist.coefCos,dist.coefSin])

dist.initUniform()
unittest(2,[dist.coefCos,dist.coefSin])

dist.initCosine([0,2,4],[0.1,2.3,-1.7])
ref1 = [dist.coefCos,dist.coefSin]
dist.initSine([1,3],[0.2,0.9])
ref2 = [dist.coefCos,dist.coefSin]
unittest(3,ref1+ref2)

dist.initCosine([0,1,3],[0.1,2.3,-1.7],[1,2,3])
ref1 = [dist.coefCos,dist.coefSin]
dist.initSine([0,1,3],[0.1,2.3,-1.7],[1,2,3])
ref2 = [dist.coefCos,dist.coefSin]
unittest(4,ref1+ref2)

dist.initProbability([2],[0.3])
unittest(5,[dist.coefCos,dist.coefSin])

dist.initProbability([1,1],[0.3,2.7])
unittest(6,[dist.coefCos,dist.coefSin])

dist.initProbability([1,2,1,5],[0.3,2.7,-9.2,1.2])
unittest(7,[dist.coefCos,dist.coefSin])

dist.initProbability([1,2,5],[0.3,2.7,-9.2])
unittest(8,[dist.coefCos,dist.coefSin])

dist.initNormal(0.4,0.8,10)
unittest(9,[dist.coefCos,dist.coefSin])

dist.initNormal(0.4,0.8,10)
ref = dist.evaluate(np.arange(100) * (2*np.pi/100.))
unittest(10,[ref])

dist.initProbability([1,2,1,5],[0.3,2.7,-9.2,1.2])
unittest(11,[dist.mean(), dist.variance(), dist.std()])

unittest(12,[np.asarray(dist.normalError(sigma=0.2,kMax=7))])

unittest(13,[np.asarray(dist.normalCriticalSigma(epsilon=1e-5,kMax=7))])


dist1 = bpe.DensityFourier()
dist2 = bpe.DensityFourier()
dist1.initOne()
dist2.initOne()
dist3 = dist1 * dist2
unittest(14,[dist3.coefCos,dist3.coefSin])
dist = dist1 # Must be the same
if ((np.max(np.abs(dist.coefCos - dist3.coefCos)) > 1e-8) or
    (np.max(np.abs(dist.coefSin - dist3.coefSin)) > 1e-8)) :
    print("Test 14: INTERNAL CHECK FAILED");


dist1 = bpe.DensityFourier()
dist2 = bpe.DensityFourier()
dist1.initProbability([1],[0.3])
dist2.initOne()
dist3 = dist1 * dist2
unittest(15,[dist3.coefCos,dist3.coefSin])
dist = dist1 # Must be the same
if ((np.max(np.abs(dist.coefCos - dist3.coefCos)) > 1e-8) or
    (np.max(np.abs(dist.coefSin - dist3.coefSin)) > 1e-8)) :
    print("Test 15: INTERNAL CHECK FAILED");


dist1 = bpe.DensityFourier()
dist2 = bpe.DensityFourier()
dist1.initProbability([1],[0.3])
dist2.initProbability([1],[2.5])
dist3 = dist1 * dist2
unittest(16,[dist3.coefCos,dist3.coefSin])
dist.initProbability([1,1],[0.3,2.5]) # Must be the same
if ((np.max(np.abs(dist.coefCos - dist3.coefCos)) > 1e-8) or
    (np.max(np.abs(dist.coefSin - dist3.coefSin)) > 1e-8)) :
    print("Test 16: INTERNAL CHECK FAILED");


dist1 = bpe.DensityFourier()
dist2 = bpe.DensityFourier()
dist1.initProbability([1,5,8],[0.3,0.2,1.8])
dist2.initNormal(0.4,0.8,14);
dist3 = dist1 * dist2
unittest(17,[dist3.coefCos,dist3.coefSin])

