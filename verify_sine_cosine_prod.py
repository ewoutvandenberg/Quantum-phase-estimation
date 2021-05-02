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
import matplotlib.pyplot as plt

phi = np.linspace(0,2*np.pi,10001)

(fig, axs) = plt.subplots(2, 2)

theta = 0.2
axs[0,0].plot(phi,np.cos(theta-phi),'b-',phi,np.cos(theta)*np.cos(phi) + np.sin(theta)*np.sin(phi),'r:')
axs[0,1].plot(phi,np.cos(theta+phi),'b-',phi,np.cos(theta)*np.cos(phi) - np.sin(theta)*np.sin(phi),'r:')
axs[1,0].plot(phi,np.sin(theta+phi),'b-',phi,np.sin(theta)*np.cos(phi) + np.cos(theta)*np.sin(phi),'r:')
axs[1,1].plot(phi,np.sin(theta-phi),'b-',phi,np.sin(theta)*np.cos(phi) - np.cos(theta)*np.sin(phi),'r:')
plt.show()

