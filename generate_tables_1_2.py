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

# Make sure the 'tables' directory exists
ensureDirExists('tables')


# Table 1
filename = 'tables/TableSpurious.tex'
print("Writing table 1: %s . . ." % filename)
with open(filename,'w') as fp :
   fp.write("\\begin{tabular}{|lllll|rrrr|}\n")
   fp.write("\\hline\n")
   fp.write("spurious & additional & \multicolumn{3}{c}{\\# failed} & \\multicolumn{2}{c}{phase error} & \\multicolumn{2}{c|}{weight error}\\\\\n")
   fp.write("phases & parameters & \multicolumn{3}{c}{instances} & maximum & median & maximum & median\\\\\n")
   fp.write("\\hline\n")
   for (i,nSpurious) in enumerate([2,5,10,20]) :
      for j in range(5) :
         try :
            data = loadInstances('./cache/experiment_hybrid_%d.dat' % (310 + i*5 + j))
            nFail = [0,0,0]
            phiErr = []
            weightErr = []
            for inst in data :
               for [tIdx,threshold] in enumerate([1,3,5]) :
                  result = filterPhaseData(inst,3,threshold=threshold)
                  if (result is None) :
                     nFail[tIdx] += 1
                  elif (threshold == 3) :
                     phiErr.extend(result[0])
                     weightErr.extend(result[1])
         except :
            nFail = 0
            data = []
            phiErr = [0]
            weightErr =[0]
         
         fp.write("%d & %d & %s & %s & %s & %s & %s\\\\\n" %
               (nSpurious,j+1,"&".join([str(f) for f in nFail]),
                "--" if (nFail == len(data)) else ("%.5f" % (np.max(phiErr))),
                "--" if (nFail == len(data)) else ("%.5f" % (np.median(phiErr))),
                "--" if (nFail == len(data)) else ("%.5f" % np.max(weightErr)),
                "--" if (nFail == len(data)) else ("%.5f" % np.median(weightErr))))
      fp.write("\\hline\n")
   fp.write("\\end{tabular}")


# Table 2
def countRecovery(data, n, threshold) :
   phiExact = 0
   for inst in data :
      result = filterPhaseData(inst, n, flagPlot=False, threshold=threshold)
      if (result is not None) :
         (phiErr, weightErr) = result
         if (all(np.asarray(phiErr) < 0.005)) :
            phiExact += 1
   return phiExact

filename = 'tables/TableMultiple.tex'
print("Writing table 2: %s . . ." % filename)
with open(filename,'w') as fp :
   # Header
   fp.write('\\begin{tabular}{|l|l|%s|}\n' % ('l'*12))
   fp.write('\\hline\n')
   fp.write('phase&&\\multicolumn{12}{|c|}{$n$}\\\\\n')
   fp.write('parameters & $\\tau$')
   for i in range(1,13) :
      fp.write('& %2d' % i)
   fp.write("\\\\\n")
   fp.write('\\hline\n')

   for threshold in [0.01,1,3,5] :
      # Body
      for j in range(0,4) :
         fp.write('$n$+%d & %s' % (j, "" if (j > 0) else str(threshold)))
         for i in range(1,13) :
            #try :
            if (True) :
               data = loadInstances('./cache/experiment_hybrid_2%02d.dat' % (i+j*20))
               v = "%3d" % countRecovery(data, i, threshold)
            #except :
            if (False) :
               v = '---'
            fp.write("& %s" % v)
         fp.write("\\\\\n")
      fp.write('\\hline\n')

   # Footer
   fp.write('\\end{tabular}')


