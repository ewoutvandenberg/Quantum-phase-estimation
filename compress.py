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

import sys

def writeBuffer(fpDst, buffer, options, instance) :
   if (buffer) :
      if (instance in options['full']) :
         for line in buffer:
             fpDst.write('%s' % line)
      else :
         # Count the number of lines starting with '*'
         nAsterisk = 0; nWeight = 0
         for line in buffer :
            if (line[0] == '*') :
               nAsterisk += 1
            elif (line[0] == 'w') :
               nWeight += 1

         idxAsterisk = 0; idxWeight = 0
         for line in buffer :
            if (line[0] == '*') :
               idxAsterisk += 1
               if (idxAsterisk < nAsterisk-26) :
                   continue
            elif (line[0] == 'w') :
               idxWeight += 1
               if (idxWeight < nWeight) :
                  continue
            fpDst.write('%s' % line)
            
             
def compressFile(filename, options) :
   buffer = []
   instance = 0

   # First read the entire source file
   with open(filename,'r') as fp :
      lines = fp.readlines()
      
   with open(filename,'w') as fp :
      for line in lines :
            buffer.append(line)
            if (line[0] == '!') :
               writeBuffer(fp, buffer, options, instance)
               buffer = []; instance += 1
      writeBuffer(fp, buffer, options, instance)


# Filter data for the 200 range problems
for i in range(4) :
   for j in range(1,13) :
      index = 200+i*20+j
      options = {}
      if (index == 269) :
         options['full'] = [1]
      else  :
         options['full'] = []
      print("Compressing %d . . ." % index)
      compressFile('./cache/experiment_hybrid_%d.dat' % index, options)


# Filter data for the 300 range problems
for index in range(310,329+1) :
      options = {}
      if (index == 314) :
         options['full'] = [15]
      else  :
         options['full'] = []
      print("Compressing %d . . ." % index)
      compressFile('./cache/experiment_hybrid_%d.dat' % index, options)
