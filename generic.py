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
import matplotlib.text as mtext
from   matplotlib import colors
import numpy as np
import bpe

import subprocess
import os



# ============================================================================
# Directories, makefile, and running an experiment
# ============================================================================
cachedir = "cache/"

def ensureDirExists(directory) :
   try :
      os.makedirs(directory)
   except FileExistsError :
      pass

def ensureCacheDir() :
   ensureDirExists(cachedir)

def cacheFilename(filename) :
   return os.path.join(cachedir, filename)


def makeExperiment(name) :
   print("Running make %s . . ." % name)
   result = subprocess.run(["make",name], capture_output=True)
   if (result.returncode != 0) :
      print("=== Error running make ===\n%s" % (result.stdout.decode("utf-8")))
      print("%s" % (result.stderr.decode("utf-8")))
      raise Exception("Error running make")

def runExperiment(name, param=[], verbose=False, showinfo=True) :
   if (showinfo) :
      print("Running %s %s" % (name, " ".join([p for p in param+[". . ."]])))
   result = subprocess.run([name] + param, capture_output=True)
   if (result.returncode != 0) :
      print("=== Error running experiment %s ===\n%s" % (name, result.stdout.decode("utf-8")))
      raise Exception("Error running experiment")
   if (verbose):
      print("%s" % result.stdout.decode("utf-8"))
   else :
      return result.stdout.decode("utf-8")



# ============================================================================
# Parsing of results
# ============================================================================

def parseDistribution(s) :
   distribution = bpe.distribution();
   distribution.parse(s);
   return distribution


# ============================================================================
# Plotting
# ============================================================================

def setTickFontsize(fontsize) :
    for tick in plt.gca().yaxis.get_major_ticks():
      tick.label.set_fontsize(fontsize)
    for tick in plt.gca().xaxis.get_major_ticks():
      tick.label.set_fontsize(fontsize)

class StringObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        h = mtext.Text(x=x0+width/2,y=y0+height,text=orig_handle,fontweight='semibold',fontsize=fontsize,
                       verticalalignment='top',horizontalalignment='center')
        handlebox.add_artist(h)
        return h

class ColorStringObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height

        s = orig_handle
        idx = s.find('@')
        if (idx == -1) :
           h = mtext.Text(x=x0+width/2,y=y0+height,text=orig_handle,fontweight='semibold',fontsize=fontsize,
                          verticalalignment='top',horizontalalignment='center')
        else :
           h = mtext.Text(x=x0+width/2,y=y0+height,text=s[:idx],fontweight='semibold',fontsize=fontsize,
                          verticalalignment='top',horizontalalignment='center',color=s[idx+1:])
        handlebox.add_artist(h)
        return h

def exportFigure(prefix, closeFig=True,label='') :
   ensureDirExists('fig/')
   prefix    = 'Figure-%s' % label
   uncropped = "fig/%s-uncropped.pdf" % prefix
   cropped   = "fig/%s.pdf" % prefix

   # Export the figure
   print("")
   print("Generating figure %s" % (label))
   print("Uncropped version: %s" % (uncropped))
   plt.savefig(uncropped, transparent=True)

   # Close the figure if needed
   if (closeFig) :
      plt.close()

   print("Cropping disabled; see the exportFigure function in generic.py")
   #try :
   #   # Crop the file
   #   if (os.system("pdfcrop %s %s" % (uncropped, cropped)) == 0) :
   #      # Delete the uncropped figure
   #      os.remove(uncropped)
   #except :
   #   pass


def plotDistribution(distribution) :
   density = distribution.mass / (distribution.gridpoints[1:] - distribution.gridpoints[:-1])

   b = 0; p = 0
   for i in range(len(distribution.mass)) :
      pPrev = p
      a = b;
      b = distribution.gridpoints[i+1]
      p = distribution.mass[i] / (b - a)
      plt.plot([a,b],[p,p],'b-')
      if (i > 0) :
         plt.plot([a,a],[pPrev, p],'b:')





# ============================================================================
# Problem instance class
# ============================================================================

class ProblemInstance :
   def __init__(self, filename, s) :
      s = s.split()
      n = int(s[1])
      phiStar = np.asarray([float(s[i]) for i in range(2,2+n)])
      weightStar = np.asarray([float(s[i]) for i in range(2+n,2+2*n)]) 

      self.filename   = filename
      self.phiStar    = phiStar
      self.wStar      = weightStar
      self.weightIter = []
      self.weightVals = []
      self.phiIter    = []
      self.phiMu      = []
      self.phiSigma   = []
      self.runtime    = 0
      self.switchIter = []
      self.switchIdx  = []

   def addRuntime(self, s):
      s = s.split()
      self.runtime = float(s[1])
   
   def addSwitch(self, s):
      s = s.split()
      self.switchIter.append(int(s[1]))
      self.switchIdx.append(int(s[2]))

   def addWeight(self, s) :
      s = s.split()
      self.weightIter.append(int(s[1]))
      self.weightVals.append([float(s[i]) for i in range(2,len(s))])

   def addPhi(self, s) :
      s = s.split()
      self.phiIter.append(int(s[1]))
      self.phiMu.append([float(s[i]) for i in range(2,len(s),2)])
      self.phiSigma.append([float(s[i]) for i in range(3,len(s),2)])

   def convert(self) :
      self.weightIter = np.asarray(self.weightIter)
      self.weightVals = np.asarray(self.weightVals)
      self.phiIter    = np.asarray(self.phiIter)
      self.phiMu      = np.asarray(self.phiMu)
      self.phiSigma   = np.asarray(self.phiSigma)

   def match(self, matchIter=-1) :
       self.phiMatch = np.zeros(self.phiMu.shape[1],dtype=np.double)
       for i in range(self.phiMu.shape[1]) :
          dMin = 10
          for j in range(self.phiStar.size) :
             d = np.abs(self.phiStar[j] - self.phiMu[matchIter,i])
             d = np.minimum(d, 2*np.pi-d)
             if (d < dMin) :
                dMin = d
                self.phiMatch[i] = self.phiStar[j]


   def collatePerIter(self, matchIter=None) :
      # Note: we ignore matchIter, collatedWeights are not set!

      m = self.phiStar.size
      n = self.phiMu.shape[0]
      l = self.phiMu.shape[1]

      # Map weight values to phi
      self.phiWeights = np.zeros(self.phiMu.shape, dtype=np.double)
      phiIdx = 0;

      for (weightIdx,weightIter) in enumerate(self.weightIter) :
         while ((phiIdx < n) and (self.phiIter[phiIdx] <= weightIter)) :
            self.phiWeights[phiIdx,:] = self.weightVals[weightIdx,:]
            phiIdx += 1

      self.collatedWeights = None
      self.collatedPhi     = np.zeros((n,self.phiStar.size),dtype=np.double)
      self.collatedSigma   = np.ones((n,self.phiStar.size),dtype=np.double) * 20

      dist = np.zeros((m,l),dtype=np.double)
      for i in range(n) :
         # Compute the distances
         for j in range(m) :
            dist[j,:] = np.abs(self.phiMu[i,:] - self.phiStar[j])
         idx = dist > np.pi;
         dist[idx] = 2*np.pi - dist[idx]
         minDist = np.min(dist,axis=0)

         for j in range(m) :
            # Allect all elements closest to phiStar[j]
            phi = 0
            w   = 0
            s   = 0
            for k in range(l) :
               if (dist[j,k] == minDist[k]) :
                  w   += self.phiWeights[i,k]
                  phi += self.phiWeights[i,k] * np.exp(1j * self.phiMu[i,k])
                  s   += self.phiWeights[i,k] * self.phiSigma[i,k]
            if (w > 0) :
               phi = np.angle(phi)
               self.collatedPhi[i,j] = phi if (phi > 0) else 2*np.pi + phi
               self.collatedSigma[i,j] = s / w


   def collate(self, matchIter=-1) :
       self.phiMatch = np.zeros(self.phiMu.shape[1],dtype=np.double)
       index1 = np.zeros(self.phiMu.shape[1],dtype=int)
       index2 = np.zeros(self.phiMu.shape[1],dtype=int)

       if (matchIter != -1) :
          matchIdx = int(np.max(np.argwhere(self.phiIter <= matchIter)))
       else :
          matchIdx = -1
     
       for i in range(self.phiMu.shape[1]) :
          dMin = 10
          for j in range(self.phiStar.size) :
             d = np.abs(self.phiStar[j] - self.phiMu[matchIdx,i])
             d = np.minimum(d, 2*np.pi-d)
             if (d < dMin) :
                dMin = d
                index1[i] = j

       idx = 0
       for j in range(self.phiStar.size) :
          for i in range(self.phiMu.shape[1]) :
             if (index1[i] == j) :
               index2[idx] = i
               self.phiMatch[idx] = self.phiStar[j]
               idx += 1

       # Map weight values to phi mapping
       self.phiWeights = np.zeros(self.phiMu.shape, dtype=np.double)

       n = self.weightVals.shape[0]
       for i in range(n) :
          if (i+1 < n) :
            self.phiWeights[self.weightIter[i]:self.weightIter[i+1],:] = self.weightVals[i,:]
          else :
            self.phiWeights[self.weightIter[i]:,:] = self.weightVals[i,:]
          

       self.collatedWeights = np.zeros((self.phiMu.shape[0],self.phiStar.size),dtype=np.double)
       self.collatedPhi     = np.zeros((self.phiMu.shape[0],self.phiStar.size),dtype=np.double)
       self.collatedSigma   = np.zeros((self.phiMu.shape[0],self.phiStar.size),dtype=np.double)
       for i in range(self.phiMu.shape[1]) :
          idx = index1[i]
          p = np.copy(self.phiMu[:,i])
          d = np.abs(self.collatedPhi[:,idx] - p)
          p[d > np.pi] -= 2*np.pi
          w = self.collatedWeights[:,idx] + self.phiWeights[:,i]
          d = (self.collatedWeights[:,idx] * self.collatedPhi[:,idx] + self.phiWeights[:,i] * p)
          s = (self.collatedWeights[:,idx] * self.collatedSigma[:,idx] + self.phiWeights[:,i] * self.phiSigma[:,i])
          wIdx = (w > 0)
          d[wIdx] /= w[wIdx]
          s[wIdx] /= w[wIdx]
          d[d < 0] += 2*np.pi
          self.collatedWeights[:,idx] = w
          self.collatedPhi[:,idx] = d
          self.collatedSigma[:,idx] = s

       if (self.weightVals.size > 0) :
          self.collatedWeights = np.zeros((self.weightVals.shape[0],self.phiStar.size),dtype=np.double)
          for i in range(self.weightVals.shape[1]) :
             idx = index1[i]
             self.collatedWeights[:,idx] += self.weightVals[:,i]


   def collateMaxWeight(self, matchIter=-1) :
       self.phiMatch = np.zeros(self.phiMu.shape[1],dtype=np.double)
       index1 = np.zeros(self.phiMu.shape[1],dtype=int)
       index2 = np.zeros(self.phiMu.shape[1],dtype=int)

       if (matchIter != -1) :
          matchIdx = int(np.max(np.argwhere(self.phiIter <= matchIter)))
       else :
          matchIdx = -1

       # index1[i] = j: estimate i maps to phi state j
       for i in range(self.phiMu.shape[1]) :
          dMin = 10
          for j in range(self.phiStar.size) :
             d = np.abs(self.phiStar[j] - self.phiMu[matchIdx,i])
             d = np.minimum(d, 2*np.pi-d)
             if (d < dMin) :
                dMin = d
                index1[i] = j

       idx = 0
       for j in range(self.phiStar.size) :
          for i in range(self.phiMu.shape[1]) :
             if (index1[i] == j) :
               index2[idx] = i
               self.phiMatch[idx] = self.phiStar[j]
               idx += 1

       # Map weight values to phi mapping
       self.phiWeights = np.zeros(self.phiMu.shape, dtype=np.double)

       n = self.weightVals.shape[0]
       for i in range(n) :
          if (i+1 < n) :
            self.phiWeights[self.weightIter[i]:self.weightIter[i+1],:] = self.weightVals[i,:]
          else :
            self.phiWeights[self.weightIter[i]:,:] = self.weightVals[i,:]
          

       self.collatedWeights = np.zeros((self.phiMu.shape[0],self.phiStar.size),dtype=np.double)
       self.collatedPhi     = np.zeros((self.phiMu.shape[0],self.phiStar.size),dtype=np.double)
       self.collatedSigma   = np.zeros((self.phiMu.shape[0],self.phiStar.size),dtype=np.double)

       for j in range(self.phiStar.size) :
          idx = -1; w = 0
          for i in range(self.phiMu.shape[1]) :
             if ((index1[i] == j) and (self.phiWeights[matchIdx,i] > w)) :
                self.collatedPhi[:,j] = self.phiMu[:,i]
                self.collatedSigma[:,j] = self.phiSigma[:,i]
                self.collatedWeights[:,j] = self.phiWeights[:,i]
                w = self.phiWeights[matchIdx,i]
               


def loadInstances(filename,max=-1) :
   instances = []; problem = None
   with open(filename,'r') as f :
      for s in f :
         if (not s) : break

         if (s[0] == '!') :
            problem.addRuntime(s)
            problem.convert()
            problem.collate(matchIter=-1)
            instances.append(problem)
            problem = None
            if (len(instances) == max) :
               break

         if (s[0] == 'w') :
            problem.addWeight(s)
            
         if (s[0] == '*') :
            problem.addPhi(s)

         if (s[0] == '#') :
            problem = ProblemInstance(filename,s)

         if (s[0] == 's') :
            problem.addSwitch(s)
            
   return instances


def plotFile(instances,linestyle,alpha) :
   problem = instances[0]
  
   error = np.zeros((len(instances),problem.phiMu.shape[0],problem.phiMu.shape[1]),dtype=np.double)
   for i in range(len(instances)) :
      problem = instances[i]
      error[i,:,:] = np.abs(problem.phiMu - problem.phiMatch)
   medianError = np.nanmedian(error,axis=0)
   for i in range(medianError.shape[1]) :
      c = colors.to_rgba('C%d'%i)
      c = tuple([(alpha * c[i] + (1-alpha) * 1) for i in range(3)])
      h = plt.plot(problem.phiIter,medianError[:,i],linestyle=linestyle,color=c)

   plt.title(problem.filename)

 
def plotSigma(instances,linestyle,alpha) :
   problem = instances[0]
   sigma = np.zeros((len(instances),problem.phiMu.shape[0],problem.phiMu.shape[1]),dtype=np.double)
   for i in range(len(instances)) :
      problem = instances[i]
      sigma[i,:,:] = problem.phiSigma

   medianSigma = np.nanmedian(sigma,axis=0)
   for i in range(medianSigma.shape[1]) :
      c = colors.to_rgba('C%d'%i)
      c = tuple([(alpha * c[i] + (1-alpha) * 1) for i in range(3)])
      h = plt.plot(problem.phiIter,medianSigma[:,i],linestyle=linestyle,color=c)

   plt.title(problem.filename)


def plotCollatedMu(instances,linestyle,alpha=None,style='median',color=None) :
   problem = instances[0]
   error = np.zeros((len(instances),problem.phiMu.shape[0],problem.phiStar.size),dtype=np.double)
   handles = []
   
   for i in range(len(instances)) :
      problem = instances[i]
      d = np.abs(problem.collatedPhi - problem.phiStar)
      idx = d > np.pi
      d[idx] = 2*np.pi - d[idx]
      error[i,:,:] = d


   if (style == 'median') :
      medianError = np.nanmedian(error,axis=0)
      for i in range(medianError.shape[1]) :
         c = colors.to_rgba('C%d'%i) if (color is None) else colors.to_rgba(color)
         c = tuple([(alpha * c[i] + (1-alpha) * 1) for i in range(3)])
         h = plt.plot(problem.phiIter,medianError[:,i],linestyle=linestyle,color=c)
         handles.append(h)

   elif (style == 'mean') :
      meanError = np.nanmean(error,axis=0)
      for i in range(meanError.shape[1]) :
         c = colors.to_rgba('C%d'%i) if (color is None) else colors.to_rgba(color)
         c = tuple([(alpha * c[i] + (1-alpha) * 1) for i in range(3)])
         h = plt.plot(problem.phiIter,meanError[:,i],linestyle=linestyle,color=c)
         handles.append(h)

   elif (style == 'all') :
      for i in range(error.shape[2]) :
         for j in range(error.shape[0]) :
            c = colors.to_rgba('C%d'%i) if (color is None) else colors.to_rgba(color)
            c = tuple([(alpha * c[i] + (1-alpha) * 1) for i in range(3)])
            h = plt.plot(problem.phiIter,error[j,:,i],linestyle=linestyle,color=c)
         handles.append(h)

   plt.xscale('log')
   plt.yscale('log')
   return h


def plotCollatedSigma(instances,linestyle,alpha,style='median',color=None) :
   problem = instances[0]
   sigma = np.zeros((len(instances),problem.phiSigma.shape[0],problem.phiStar.size),dtype=np.double)
   for i in range(len(instances)) :
      problem = instances[i]
      s = np.copy(problem.collatedSigma)
      s[s==0] = 20
      sigma[i,:,:] = s


   if (style == 'median') :
      medianSigma = np.nanmedian(sigma,axis=0)
    
      handles = []
      for i in range(medianSigma.shape[1]) :
         c = colors.to_rgba('C%d'%i) if (color is None) else colors.to_rgba(color)
         c = tuple([(alpha * c[i] + (1-alpha) * 1) for i in range(3)])
         h = plt.plot(problem.phiIter,medianSigma[:,i],linestyle=linestyle,color=c)
         handles.append(h)
         
   elif (style == 'all') :
      handles = []
      for k in range(sigma.shape[0]) :
       for i in range(sigma.shape[2]) :
         c = colors.to_rgba('C%d'%i) if (color is None) else colors.to_rgba(color)
         c = tuple([(alpha * c[i] + (1-alpha) * 1) for i in range(3)])
         h = plt.plot(problem.phiIter,sigma[k,:,i],linestyle=linestyle,color=c)
         handles.append(h)

   plt.xscale('log')
   plt.yscale('log')
   return handles


def plotCollatedWeights(instances,linestyle,alpha,style,color=None) :
   problem = instances[0]
   error = np.zeros((len(instances),problem.collatedWeights.shape[0],problem.wStar.size),dtype=np.double)
   for i in range(len(instances)) :
      problem = instances[i]
      error[i,:,:] = np.abs(problem.collatedWeights - problem.wStar)

   if (style == 'median') :
      medianError = np.median(error,axis=0)
      for i in range(medianError.shape[1]) :
         c = colors.to_rgba('C%d'%i) if (color is None) else colors.to_rgba(color)
         c = tuple([(alpha * c[i] + (1-alpha) * 1) for i in range(3)])
         h = plt.plot(problem.weightIter, medianError[:,i],linestyle=linestyle,color=c)

   elif (style == 'mean') :
      meanError = np.mean(error,axis=0)
      for i in range(meanError.shape[1]) :
         c = colors.to_rgba('C%d'%i) if (color is None) else colors.to_rgba(color)
         c = tuple([(alpha * c[i] + (1-alpha) * 1) for i in range(3)])
         h = plt.plot(problem.weightIter, meanError[:,i],linestyle=linestyle,color=c)

   elif (style == 'all') :
      for i in range(error.shape[2]) :
         for j in range(error.shape[0]) :
            c = colors.to_rgba('C%d'%i) if (color is None) else colors.to_rgba(color)
            c = tuple([(alpha * c[i] + (1-alpha) * 1) for i in range(3)])
            h = plt.plot(problem.weightIter, error[j,:,i],linestyle=linestyle,color=c)

   plt.xscale('log')
   plt.yscale('log')


def plotCollatedKMax(instances,linestyle,alpha=1.0,color=None,kMaxValue=50) :
   problem = instances[0]
   kMax    = np.zeros((len(instances),problem.phiSigma.shape[0]),dtype=np.double)
   for i in range(len(instances)) :
      problem = instances[i]
      n = problem.collatedSigma.shape[0]
      for j in range(n) :
         # Find the corresponding weight index
         weightIdx = np.argwhere(problem.weightIter <= j)
         weightIdx = weightIdx[-1]
         w = np.copy(problem.collatedWeights[weightIdx,:]); w[w==0] = 2
         idx = np.argmin(w)
         s = problem.collatedSigma[j,idx]
         kMax[i,j] = np.minimum(np.ceil(1.25 / s), kMaxValue)

   medianKMax = np.median(kMax,axis=0)
   c = colors.to_rgba('C0') if (color is None) else colors.to_rgba(color)
   c = tuple([(alpha * c[i] + (1-alpha) * 1) for i in range(3)])
   handles = plt.plot(problem.phiIter,medianKMax,linestyle=linestyle,color=c)

   plt.xscale('log')
   return handles


def plotKVals(instances,linestyle='-',alpha=1.0,color=None,kMaxValue=50) :
   problem = instances[0]
   kVals   = np.zeros((len(instances),problem.phiIter.shape[0]),dtype=np.double)
   for i in range(len(instances)) :
      problem = instances[i]

      n = problem.phiSigma.shape[0]
      for j in range(n) :
         # Find the corresponding weight index
         weightIdx = np.argwhere(problem.weightIter <= problem.phiIter[j])
         weightIdx = weightIdx[-1]
         kVals[i,j] = np.sum(problem.weightVals[weightIdx,:]*np.ceil(1.25 / problem.phiSigma[j,:]))

      kVals = np.minimum(np.ceil(kVals), kMaxValue)

   medianKVals = np.median(kVals,axis=0)
   c = colors.to_rgba('C0') if (color is None) else colors.to_rgba(color)
   c = tuple([(alpha * c[i] + (1-alpha) * 1) for i in range(3)])
   handles = plt.plot(problem.phiIter,medianKVals,linestyle=linestyle,color=c)

   plt.xscale('log')
   return handles


# ============================================================================
# Remove oscillating and low-weight estimates, combine nearby pairs
# ============================================================================

def angleDiff(a,b) :
   # Both a and b in [0,2*pi)
   v = np.abs(b-a)
   return np.minimum(v, 2*np.pi - v)

def filterPhaseData(inst, n, flagPlot=False, threshold=5) :
   # Determine the 'valid' phases
   nPhi = inst.phiMu.shape[1]
   wMax = np.max(inst.weightVals[-1,:])
   valid = np.ones(nPhi,dtype=np.bool)
   for i in range(nPhi) :
      totalVariation = np.sum(np.abs(inst.phiMu[-26:-1,i] - inst.phiMu[-25:,i]))
      if ((totalVariation > 0.5) or
          (inst.weightVals[-1,i] < 0.1*wMax)) :
         valid[i] = False

   # Combine close pairs
   if (threshold > 0) :
      threshold *= (2*np.pi / 360) # Convert to radians
   
      phiCombined = []
      weightCombined = []
      for i in range(inst.phiMu.shape[1]) :
         mu     = inst.phiMu[-1,i]
         weight = inst.weightVals[-1,i]
         if (not valid[i]) : continue;

         # Find the phase closest to the current phase
         minDiff  = 2*np.pi
         minIndex = -1
         for (j,muRef) in enumerate(phiCombined) :
            diff = angleDiff(mu,muRef[-1])
            if (diff < minDiff) :
               minDiff = diff
               minIndex = j

         if ((minIndex == -1) or (minDiff > threshold)) :
            # New or not close to any previous cluster
            phiCombined.append(inst.phiMu[:,i])
            weightCombined.append(inst.weightVals[-1,i])
         else :
            muSum = 0; weightSum = 0
            refAngle = np.minimum(inst.phiMu[:,i],phiCombined[minIndex])
            for (mu,weight) in [(inst.phiMu[:,i],inst.weightVals[-1,i]),
                                (phiCombined[minIndex], weightCombined[minIndex])] :
               mu = np.copy(mu)
               idx = ((mu - refAngle) > 2*np.pi)
               mu[idx] -= 2*np.pi

               # Add up the variables (data type changes)
               weightSum = weightSum + weight
               muSum     = muSum + weight * mu

            muSum = muSum / weightSum
            idx = (muSum < 0)
            muSum[idx] += 2*np.pi
            phiCombined[minIndex] = muSum
            weightCombined[minIndex] = weightSum
   else :
      # Do not aggregate the estimates
      phiCombined    = [[inst.phiMu[-1,idx]] for idx in range(inst.phiMu.shape[1])]
      weightCombined = inst.weightVals[-1,:]

   if (flagPlot) :
      for phi in phiCombined :
         plt.plot(inst.phiIter, phi, linewidth=2)
      for i in range(3) :
         plt.plot([1,inst.phiIter[-1]],
                  [inst.phiStar[i],inst.phiStar[i]],'k--',alpha=0.5, zorder=-10)
      plt.ylim([0,6.5])
      plt.xscale('log')


   # Final values after sorting
   phiCombined = [(phi[-1],idx) for (idx,phi) in enumerate(phiCombined)]
   phiCombined.sort()
   weightCombined = [weightCombined[idx] for (_,idx) in phiCombined]
   phiCombined = [phi for (phi,_) in phiCombined]

   # Assume the reference phi Star is already sorted
   if (len(phiCombined) != n) :
      return None
   else :
      phiStar = np.sort(inst.phiStar[:n])
      phiErr = [angleDiff(phiCombined[i],phiStar[i]) for i in range(n)]
      weightErr = [np.abs(weightCombined[i] - inst.wStar[i]) for i in range(n)]
      return (phiErr, weightErr)


def plotInstanceFilter(inst, n) :
   wMax = np.max(inst.weightVals[-1,:])
   colorIdx=0
   for i in range(inst.phiMu.shape[1]) :
      if ((np.sum(np.abs(inst.phiMu[-26:-1,i] - inst.phiMu[-25:,i])) > 1) or
          (inst.weightVals[-1,i] < 0.05*wMax)) :
         color = (0.85,0.85,0.88)
         zorder = -10
         linewidth = 1
      else :
         color = 'C%d' % colorIdx
         colorIdx += 1
         zorder = None
         linewidth = 2
      plt.plot(inst.phiIter,inst.phiMu[:,i], color=color,linewidth=linewidth,zorder=zorder)

   for i in range(n) :
      plt.plot([0,inst.phiIter[-1]],
               [inst.phiStar[i],inst.phiStar[i]],'k--',alpha=0.5, zorder=-10)
   plt.ylim([0,6.5])
   plt.xscale('log')



# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------

np.seterr(invalid='ignore')
fontsize = 14
