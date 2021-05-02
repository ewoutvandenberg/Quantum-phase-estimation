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

import bpe
import matplotlib.pyplot as plt
import numpy as np

angle = np.linspace(0,2*np.pi,1001,dtype=np.float)
values= np.zeros(angle.size,angle.dtype)

for i in range(angle.size) :
   # Probability of measuring 1 = (1 - cos(angle)) / 2
   values[i] = 1 - 2 * (np.sum(bpe.measure(angle[i],1000)) / 1000.0)

plt.plot(angle,values,'b.',angle,np.cos(angle),'r-')
plt.show()
