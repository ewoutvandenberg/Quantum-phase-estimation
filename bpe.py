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
import scipy.special


M_PI = float(np.pi)

class distribution :
   def __init__(self) :
      self.gridpoints = []
      self.mass       = []

   def parse(self, s) :
      s = s.split()
      n = int(s[0])
      self.gridpoints = np.asarray([float(v) for v in s[1:n+2]])
      self.mass       = np.asarray([float(v) for v in s[n+2:2*n+2]])
      self.normalize()

   def normalize(self) :
      self.mass /= np.sum(self.mass)

   def uniformPrior(self, n) :
      delta = 2*M_PI / n
      self.gridpoints = np.linspace(0,2*M_PI,n+1)
      self.mass = np.ones(n,dtype=np.double) / n

   def gaussianPrior(self, n, mu, sigma) :
      delta = 2*M_PI / n;
      self.gridpoints = np.linspace(0,2*M_PI,n+1)
      grid = (self.gridpoints[1:] + self.gridpoints[:-1]) / 2
      self.mass = np.zeros(n,dtype=np.double)
      k = max(1,int(np.ceil(8*sigma / (2*M_PI))))
      for i in range(-k,k+1) :
         self.mass += np.exp(-((mu+grid+i*2*M_PI)**2)/(2*(sigma**2)))
      self.normalize()

   def cosineProb(self,n,s,theta,value) :
      gamma = 1 if (value == 0) else -1
      delta = 2*M_PI / n;
      self.gridpoints = np.linspace(0,2*M_PI,n+1)
      grid = (self.gridpoints[1:] + self.gridpoints[:-1]) / 2
      self.mass = (1 + gamma*np.cos(s*grid-theta)) / 2.

   def getErrorSquared(self, ref) :
      error = 0
      rbar = (ref + M_PI) if (ref < M_PI) else (ref - M_PI)

      # We should ensure that the segments are sufficiently small

      b = 0
      for i in range(len(self.mass)) :
         a = b
         b = self.gridpoints[i+1]
         p = self.mass[i] / (b - a) # Probability density

         if (b - a > 3) :
            raise Exception("Grid cell exceeds maximum size")

         # Compute the signed angle differences
         da = a - ref
         db = b - ref
         opposite = ((a <= rbar) and (rbar < b))
         if (opposite) :
            if (ref < M_PI) :
               db -= 2*M_PI
            else :
               da += 2*M_PI
         else :
            if (da < -M_PI) :
               da += 2*M_PI
               db += 2*M_PI
            if (da > M_PI) :
               da -= 2*M_PI
               db -= 2*M_PI

         # Integral of d^2*p dx
         error += p * (db*db*db - da*da*da) / 3.
      
         if ((a <= rbar) and (rbar < b)) :
            # Add correction term
            error += p * (2 * M_PI * M_PI * M_PI) / 3.

      return error

  
   def getMean(self) :
      # Initialize the mean squared error
      mean = 0
      mseBest = self.getErrorSquared(mean)

      # Traverse the array to find the opposite cell
      idx = 0
      for idxOpposite in range(len(self.gridpoints)-1) :
         if (self.gridpoints[idxOpposite+1] > M_PI) :
            break

      mu = 0; mubar = M_PI; mode = 0
      segments = []
      while True :
         # Distance to end of current and opposite interval
         d0 = self.gridpoints[idx+1] - mu
         d1 = self.gridpoints[idxOpposite+1] - mubar
         if (d0 < d1) :
            delta = d0; mode = 0
         else :
            delta = d1; mode = 1

         # Initialize the polynomial coefficients
         polycoef0 = 0
         polycoef1 = 0
         polycoef2 = 0

         b = 0
         for i in range(len(self.gridpoints)-1) :
            a = b
            b = self.gridpoints[i+1]
            p = self.mass[i] / (b - a) # Probability density

            # Determine the signed distance
            da = a - mu; db = b - mu
            if (i != idxOpposite) :
               if (idx < idxOpposite) :
                  if (i > idxOpposite) :
                     da = da - 2*M_PI
                     db = db - 2*M_PI
               elif (idx > idxOpposite) :
                  if (i < idxOpposite) :
                     da = 2*M_PI + da
                     db = 2*M_PI + db
            else : # (i == idxOpposite)
               if (idx < idxOpposite) :
                  db -= 2*M_PI
               else :
                  da += 2*M_PI

            # The integral of [d - delta]^2 from d=da to d=db
            # The distance da should always be less than db
            dap = p * da
            dbp = p * db
            polycoef2 += (dbp - dap)
            dap *= da
            dbp *= db
            polycoef1 -= (dbp - dap)
            dap *= da
            dbp *= db
            polycoef0 += (dbp - dap) / 3.

            if (i == idxOpposite) :
               polycoef0 += p * (2 * M_PI * M_PI * M_PI / 3.)

         # Check for local maxima
         if (polycoef2 > 0) :
            root = -polycoef1 / (2 * polycoef2)
            if ((root >= 0) and (root <= delta)) :
               mse = polycoef0 + polycoef1 * root + polycoef2 * (root*root)
               if (mse < mseBest) :
                  mseBest = mse
                  mean = mu + root

         # Add the segment
         segments.append((mu, delta, polycoef0, polycoef1, polycoef2))

         # Update the segment
         if (mode == 0) :
            idx += 1
            if (idx >= len(self.gridpoints)-1) :
               break

            mu = self.gridpoints[idx]
            mubar = (mu + M_PI) if (mu < M_PI) else mu - M_PI
         else :
            idxOpposite = (idxOpposite + 1) % (len(self.gridpoints)-1)
            mubar = self.gridpoints[idxOpposite]
            mu = (mubar + M_PI) if (mubar <= M_PI) else mubar - M_PI

      return (mean, mseBest, segments)
      



def angleDiff(alpha, beta) :
   if ((alpha < 0) or (alpha >= 2*M_PI)) :
      alpha -= np.floor(alpha / (2*M_PI)) * (2*M_PI);
   if ((beta < 0) or (beta >= 2*M_PI)) :
      beta  -= np.floor(beta  / (2*M_PI)) * (2*M_PI);

   # Compute and normalize the difference
   diff = (beta - alpha) if (alpha < beta) else (alpha - beta);
   if (diff > M_PI) :
      diff = 2*M_PI - diff;

   return diff;




def p0(phi) :
   return (1 + np.cos(phi)) / 2

def p1(phi) :
   return (1 - np.cos(phi)) / 2

def measure(phi,n,k=1,theta=0) :
   p = (1 + np.cos(k*phi + theta)) / 2.
   return (np.random.random(n) >= p).astype(np.int8)

def pmeas(meas,phi,k=1,theta=0) :
   p0 = (1 + np.cos(k*phi + theta)) / 2.
   return ((1-2*m)*p0 + m) # ((1-m)*p0 + m*(1-p0))




# Sample multiple rounds, as in OBR2019TTa, see verify_OBR2019TTa_sample_rounds.py
def sampleRounds(A,phi,k,beta) :
   # A_j = |a_j|^2
   A = np.asarray(A,dtype=np.double)
   k = np.asarray(k,dtype=int)
   beta = np.asarray(beta,dtype=np.double)
   m = np.zeros(k.size,dtype=int)

   for i in range(k.size) :
      pZero = np.sum(A * (np.cos(k[i]*phi/2 + beta[i]/2)**2))
      m[i] = 1 * (np.random.random(1) > pZero)

      # Update A and normalize
      A = A * np.cos(k[i]*phi/2 + (beta[i] - m[i]*np.pi)/2) ** 2
      A /= np.sum(A)

   return m


# Rho from OBR2019TTA, Appendix A, see verify_OBR2019TTa_rho.py
def rhoScalar(ell,m,K) :
   K2 = K // 2
   p1 = np.arange(ell//2 + 1)
   s1 = scipy.special.binom(m,2*p1) * scipy.special.binom(K2-m,ell-2*p1) / scipy.special.binom(K2,ell)
   return 2*np.sum(s1) - 1

def rhoMatrix(K) :
   K2 = K // 2
   rho = np.zeros((K2+1,K2+1),dtype=np.double)
   for ell in range(K2+1) :
      for m in range(K2+1) :
         rho[ell,m] = rhoScalar(ell,m,K)
   return rho


# Chi from OBR2019TTA, Appendix A, see verify_OBR2019TTa_chi.py
def chiMatrix(K,k) :
   K2 = K // 2
   chi = np.zeros((K2+1,K2+1),dtype=np.complex)
   rho = rhoMatrix(K)
   
   for ell in range(k+1) :
      scale = scipy.special.binom(k,ell) * np.power(-1j,k-ell)
      chi += np.dot(scale * rho[ell,:].reshape((K2+1,1)), rho[k-ell,:].reshape((1,K2+1)))
   return chi


# -------------------------------------------------------------------------
# Fourier representation
# -------------------------------------------------------------------------


def normalize_density(density) :
   coefSin, coefCos = density
   scale = 1. / (coefCos.get(0,0) * 2 * np.pi)
   coefSin = {k : v*scale for (k,v) in coefSin.items()}
   coefCos = {k : v*scale for (k,v) in coefCos.items()}
   return (coefSin, coefCos)


def evaluate_density(density, phi) :
   coefSin, coefCos = density
   p = 0 * phi
   for (k,v) in coefSin.items() :
      p += np.sin(k*phi) * v
   for (k,v) in coefCos.items() :
      p += np.cos(k*phi) * v
   return p

def uniform_density() :
   return ({},{0: 1./(2*np.pi)})

def truncate_density(density, kMax) :
   coefSin = {k:v for (k,v) in density[0].items() if (k <= kMax)}
   coefCos = {k:v for (k,v) in density[1].items() if (k <= kMax)}
   return (coefSin,coefCos)

def update_density(density, kVec, betaVec, measVec, flagNormalize=True) :
   coefSin, coefCos = density
   for i in range(kVec.size) :
      k = int(kVec[i])
       
      newCoefSin = {}
      newCoefCos = {}

      # Coefficients
      cc = np.cos(betaVec[i] - measVec[i]*np.pi) / 4.0
      cs = np.sin(betaVec[i] - measVec[i]*np.pi) / 4.0

      # Initialize the coefficients
      for d in [coefSin,coefCos] :
         for j in d :
            idx = max(k-j,j-k)
            newCoefSin[j]   = 0
            newCoefSin[j+k] = 0
            newCoefSin[idx] = 0
            newCoefCos[j]   = 0
            newCoefCos[j+k] = 0
            newCoefCos[idx] = 0
         
      for (j,v) in coefCos.items() :
         idx  = max(k-j,j-k)
         sign = 1 if (k >= j) else -1
         newCoefCos[j]   += v / 2.0
         newCoefCos[k+j] += cc*v
         newCoefCos[idx] += cc*v
         newCoefSin[k+j] -= cs*v
         newCoefSin[idx] -= cs*v*sign
         
      for (j,v) in coefSin.items() :
         idx = max(k-j,j-k)
         sign = 1 if (k >= j) else -1
         newCoefSin[j]   += v / 2.0
         newCoefSin[k+j] += cc*v
         newCoefSin[idx] -= cc*v*sign
         newCoefCos[k+j] += cs*v
         newCoefCos[idx] -= cs*v

      del newCoefSin[0]
      coefSin = newCoefSin
      coefCos = newCoefCos

   # Normalize if needed
   if (flagNormalize) :
      return normalize_density((coefSin, coefCos))
   else :
      return (coefSin, coefCos)

# Fourier representation for (wrapped) normal distribution
def fourier_coef_normal(mu,sigma,kMax,wrapped=True) :
   coefSin = {}
   coefCos = {}

   if (wrapped) :
      for k in range(kMax+1) :
         #t = np.exp(1j*mu*k - (sigma*k)**2 / 2)
         #coefSin[k] = np.imag(t)
         #coefCos[k] = np.real(t)
         scale = np.exp(-(sigma*k)**2 / 2)
         coefSin[k] = np.sin(mu*k) * scale
         coefCos[k] = np.cos(mu*k) * scale
   else :
      sigma = np.sqrt(2) * sigma
      for k in range(kMax+1) :
         c = np.complex((2*mu + 1j*(sigma**2)*k) / (2*sigma))
         t = (1/2.) * np.exp(-(mu/sigma)**2 + c**2) * (scipy.special.erf(2*np.pi/sigma - c) + scipy.special.erf(c))
         coefSin[k] = np.imag(t)
         coefCos[k] = np.real(t)

   # Normalize the sine and cosine functions such that square integral is one
   coefCos[0] /= 2
   for k in coefCos : coefCos[k] /= np.pi
   for k in coefSin : coefSin[k] /= np.pi
      
   return (coefSin, coefCos)



# -------------------------------------------------------------------------
# Optimization for determining the coefficients A
# -------------------------------------------------------------------------

def project_simplex(x,r=1) :
   # Sort elements in increasing order - we expect that for most
   # projections all elements will be in the support, so we start with
   # smaller values of tau. We add an offset to simplify the code as
   # tau will never be larger than the largest element in xs. The
   # projection of x and x+const are identical.
   if (r <= 0) :
      return 0. * x
   
   xs= np.sort(x + np.min(x) + r)
   s = np.sum(xs)
   k = 0
   n = xs.size
   while (True) :
      # Check if tau = x[k] is too large
      if (s - n*xs[k] < r) : break
      s -= xs[k]
      n -= 1
      k += 1

   tau = (s - r) / n
   return np.maximum(x+(np.min(x)+r-tau),0)



# -------------------------------------------------------------------------
# Distributions
# -------------------------------------------------------------------------

class DensityFourier(object) :
   # Maintain the probability density representation:
   #
   #    density(phi) = sum_{k=0}^kMax (coefCos[k] * cos(k*phi) + coefSin[k] * sin(k*phi))
   #
   # Even though coefSin[0] is irrelevant we still keep it for simpler code.

   def __init__(self, coefCos=0, coefSin=0) :
      self.coefCos = np.asarray(coefCos,dtype=np.double)
      self.coefSin = np.asarray(coefSin,dtype=np.double)

   def initZero(self, kMax=0) :
      self.coefCos = np.zeros(kMax+1,dtype=np.double)
      self.coefSin = np.zeros(kMax+1,dtype=np.double)

   def initOne(self) :
      self.coefCos = np.ones(1,dtype=np.double)
      self.coefSin = np.zeros(1,dtype=np.double)

   def initUniform(self) :
      self.coefCos = np.ones(1,dtype=np.double) / (2*np.pi)
      self.coefSin = np.zeros(1,dtype=np.double)

   def initCosine(self, k, beta, weights=None) :
      # Initialize to sum of cos(k*phi+beta)
      # Note: k values are assumed to be nonnegative
      k       = np.asarray(k, dtype=int)
      beta    = np.asarray(beta, dtype=np.double)
      weights = np.asarray(weights, dtype=np.double) if (weights is not None) else np.ones(k.size)
      cb      = np.cos(beta)
      sb      = np.sin(beta)
      self.initZero(kMax=np.max(k))
      for idx in range(k.size) :
         self.coefCos[k[idx]] += weights[idx] * cb[idx]
         self.coefSin[k[idx]] -= weights[idx] * sb[idx]

      if (self.coefSin.size > 0) :
         self.coefSin[0] = 0

   def initSine(self, k, beta, weights=None) :
      # Initialize to sum of sin(k*phi+beta)
      # Note: k values are assumed to be nonnegative
      k       = np.asarray(k, dtype=int)
      beta    = np.asarray(beta, dtype=np.double)
      weights = np.asarray(weights, dtype=np.double) if (weights is not None) else np.ones(k.size)
      cb      = np.cos(beta)
      sb      = np.sin(beta)
      self.initZero(kMax=np.max(k))
      for idx in range(k.size) :
         self.coefCos[k[idx]] += weights[idx] * sb[idx]
         self.coefSin[k[idx]] += weights[idx] * cb[idx]

      if (self.coefSin.size > 0) :
         self.coefSin[0] = 0

   def initProbability(self, k, beta) :
      # Initialize to the product of cos(k*phi/2 + beta/2)^2
      kVals = np.asarray(k, dtype=int)
      beta  = np.asarray(beta, dtype=np.double)

      if (kVals.size == 0) :
         self.initOne()
         return ;

      if (kVals.size == 1) :
         k = int(kVals[0])
         self.coefCos = np.zeros(k+1, dtype=np.double)
         self.coefSin = np.zeros(k+1, dtype=np.double)

         # We add values to index k, just in case k is zero
         self.coefCos[0]  = 0.5
         self.coefCos[k] += np.cos(beta[0])/2
         self.coefSin[k] -= np.sin(beta[0])/2
         return
      
      nterms  = 1
      coefCos = np.ones(1,dtype=np.double)
      coefSin = np.zeros(1,dtype=np.double)
      
      for idx in range(kVals.size) :
         # Multiply by (1 + cos(k[idx]*phi + beta[idx])) / 2
         k = int(kVals[idx])
         newCoefCos = np.zeros(k + nterms, dtype=np.double)
         newCoefSin = np.zeros(k + nterms, dtype=np.double)

         # Scaling factors
         cosBeta = np.cos(beta[idx]) / 4.
         sinBeta = np.sin(beta[idx]) / 4.
         
         # Terms cos(j*phi) and sin(j*phi)
         newCoefCos[:nterms] = coefCos / 2.
         newCoefSin[:nterms] = coefSin / 2.

         # Term cos((k+j)*phi)
         v  = newCoefCos[k:k+nterms]
         v += (cosBeta * coefCos)
         v += (sinBeta * coefSin)

         # Term sin((k+j)*phi)
         v  = newCoefSin[k:k+nterms]
         v -= (sinBeta * coefCos)
         v += (cosBeta * coefSin)

         # Terms cos((k-j)*phi) and sin((k-j)*phi)for (j < k)
         t1 = min(nterms,k)
         if (t1 > 0) :
            v  = newCoefCos[k:k-t1:-1]
            v += (cosBeta * coefCos[:t1])
            v -= (sinBeta * coefSin[:t1])
            
            v  = newCoefSin[k:k-t1:-1]
            v -= (sinBeta * coefCos[:t1])
            v -= (cosBeta * coefSin[:t1])

         # Term cos((k-j)*phi) and sin((k-j)*phi) for (j >= k)
         t2 = nterms - t1
         if (t2 > 0) :
            v  = newCoefCos[:t2]
            v += (cosBeta * coefCos[t1:])
            v -= (sinBeta * coefSin[t1:])

            v  = newCoefSin[:t2]
            v += (sinBeta * coefCos[t1:]) # Sign negated
            v += (cosBeta * coefSin[t1:]) # Sign negated
         

         # Update the distribution
         coefCos = newCoefCos
         coefSin = newCoefSin
         coefSin[0] = 0.
         nterms += k

      # Set the distribution
      self.coefCos = coefCos
      self.coefSin = coefSin
      
      if (self.coefSin.size > 0) :
         self.coefSin[0] = 0

   def initNormal(self, mu, sigma, kMax) :
      k = np.arange(kMax+1)
      scale = np.exp(-(sigma*k)**2 / 2)
      self.coefCos = np.cos(mu*k) * scale
      self.coefSin = np.sin(mu*k) * scale

      # --------------------------------------------------------------
      # Without wrap-around we have -- one problem for larger values
      # of kMax is that the exponent term goes to zero, while the erf
      # value goes towards infinity. Given the numerical instability
      # and the fact that we don't really need this functionality for
      # now, we omit this option.
      # --------------------------------------------------------------
      # sigma = np.sqrt(2) * sigma
      # k = np.arange(0,kMax+1)         
      # c = (2*mu + 1j*(sigma**2)*k) / (2*sigma)
      # t = (1/2.) * np.exp(-(mu/sigma)**2 + c**2) * (scipy.special.erf(2*np.pi/sigma - c) + scipy.special.erf(c))
      # self.coefCos = np.real(t)
      # self.coefSin = np.imag(t)

      # Normalize the sine and cosine functions such that square integral is one
      self.coefCos[0] /= 2
      self.coefCos    /= np.pi
      self.coefSin    /= np.pi
       
   def truncate(self, kMax) :
      if (self.coefCos.size > kMax+1) :
         self.coefCos = self.coefCos[:kMax+1]
      if (self.coefSin.size > kMax+1) :
         self.coefSin = self.coefSin[:kMax+1]

   def ensureSize(self, kMax) :
      # Ensure that the size is at least the given kMax
      if (self.coefCos.size < kMax+1) :
         c = np.zeros(kMax+1,dtype=np.double)
         c[:self.coefCos.size] = self.coefCos
         self.coefCos = c
      if (self.coefSin.size < kMax) :
         s = np.zeros(kMax+1,dtype=np.double)
         s[:self.coefSin.size] = self.coefSin
         self.coefSin = s

   def ensureSameSize(self) :
      self.ensureSize(max(self.coefCos.size,self.coefSin.size)-1)

   def clone(self) :
      return DensityFourier(np.copy(self.coefCos), np.copy(self.coefSin))

   def normalize(self) :
      scale = (self.coefCos[0] * 2 * np.pi)
      self.coefCos /= scale
      self.coefSin /= scale 

   def evaluate(self, phi) :
      p = 0. * phi
      for k in range(self.coefCos.size) :
         p += self.coefCos[k] * np.cos(k * phi)
      for k in range(1,self.coefSin.size) :
         p += self.coefSin[k] * np.sin(k * phi)
      return p

   def mean(self) :
      self.ensureSize(1)
      angle = np.angle(self.coefCos[1] + 1j * self.coefSin[1])
      if (angle < 0) : angle += 2 * np.pi
      return angle

   def variance(self) :
      # Holevo variance
      self.ensureSize(1)
      return (1. / ((np.pi**2) * ((self.coefCos[1]**2) + (self.coefSin[1]**2)))) - 1.

   def std(self) :
      # Holevo standard deviation
      return np.sqrt(self.variance())

   def normalError(self, sigma=None, kMax=None) :
      # Error bound on the absolute value between a N(0,sigma) distribution
      # and its Fourier representation using sinusoids up to kMax.
      if (kMax  is None) : kMax = max(self.coefCos.size,self.coefSin.size) - 1
      if (sigma is None) : sigma = self.holevoStd()
      return scipy.special.erfc(kMax*sigma/np.sqrt(2)) / (sigma * np.sqrt(2 * np.pi))

   def normalCriticalSigma(self, epsilon, kMax=None) :
      # Approximate the minimum sigma value for which normalError(sigma,kMax) <= epsilon
      if (kMax is None) : kMax = max(self.coefCos.size,self.coefSin.size) - 1

      # Interval search
      low = 0.0
      high = 1.0
      epsilon *= np.sqrt(2*np.pi)
      scale = kMax / np.sqrt(2.)

      # Double the interval
      while ((scipy.special.erfc(high*scale) / high) > epsilon) :
         low = high
         high *= 2

      # Bisection search until sigma is sufficiently accurate
      while (high-low > 1e-7) :
         sigma = (low+high) / 2.
         error = scipy.special.erfc(sigma*scale) / sigma
         if (error < epsilon) :
           high = sigma
         else :
           low = sigma

      # Return sigma for which error < epsilon
      return high
  

   def scale(self, weight) :
      self.coefCos *= weight
      self.coefSin *= weight

   def add(self, other) :
      v = self + other
      self.coefCos = v.coefCos
      self.coefSin = v.coefSin

   def __add__(self, other) :
      if (self.coefCos.size >= other.coefCos.size) :
         coefCos = np.copy(self.coefCos)
         coefCos[:other.coefCos.size] += other.coefCos
      else :
         coefCos = np.copy(other.coefCos)
         coefCos[:self.coefCos.size] += self.coefCos

      if (self.coefSin.size >= other.coefSin.size) :
         coefSin = np.copy(self.coefSin)
         coefSin[:other.coefSin.size] += other.coefSin
      else :
         coefSin = np.copy(other.coefSin)
         coefSin[:self.coefSin.size] += self.coefSin

      return DensityFourier(coefCos,coefSin)

   def __sub__(self, other) :
      if (self.coefCos.size >= other.coefCos.size) :
         coefCos = np.copy(self.coefCos)
         coefCos[:other.coefCos.size] -= other.coefCos
      else :
         coefCos = -1 * other.coefCos
         coefCos[:self.coefCos.size] += self.coefCos

      if (self.coefSin.size >= other.coefSin.size) :
         coefSin = np.copy(self.coefSin)
         coefSin[:other.coefSin.size] -= other.coefSin
      else :
         coefSin = -1 * other.coefSin
         coefSin[:self.coefSin.size] += self.coefSin

      return DensityFourier(coefCos,coefSin)

   def __mul__(self, other) :
      if (not isinstance(other,DensityFourier)) :
         scale = float(other)
         return DensityFourier(self.coefCos * scale, self.coefSin * scale)
         
      # Multiply by another Fourier-based density function - we temporarily
      # work with coefficient arrays that allow negative values for efficient
      # processing.
      n1 = max(self.coefCos.size,self.coefSin.size)
      n2 = max(other.coefCos.size,other.coefSin.size)
      k  = n1 + n2 - 2; offset = k+1

      # We add additional entry before the most negative index
      # to avoid indexing with a range #:-1:-1, which does not
      # works as intended due to the Python contention that -1
      # indexes the last entry.
      coefCos = np.zeros(2*k+2,dtype=np.double)
      coefSin = np.zeros(2*k+2,dtype=np.double)

      def __mul_helper__(dest, offset, coef1, coef2, signs) :
         # Sign variables contains signs for the j1+j2 and j1-j2 index
         # terms.

         n = coef2.size
         for j in range(coef1.size) :
            if (coef1[j] != 0) :
               dest[offset+j:offset+j+n]    += (signs[0]*coef1[j]) * coef2
               dest[offset+j:offset+j-n:-1] += (signs[1]*coef1[j]) * coef2
      
      # Product of self.cosine and other.cosine terms
      __mul_helper__(coefCos, offset, self.coefCos, other.coefCos, (1,1))

      # Product of self.cosine and other.sine terms
      __mul_helper__(coefSin, offset, self.coefCos, other.coefSin, (1,-1))
         
      # Product of self.sine and other.cosine terms
      __mul_helper__(coefSin, offset, self.coefSin, other.coefCos, (1,1))
          
      # Product of self.sine and other.sine terms
      __mul_helper__(coefCos, offset, self.coefSin, other.coefSin, (-1,1))

      # Map the negative indices back to the positive ones
      coefCos[offset+1:] += np.flip(coefCos[1:offset])
      coefSin[offset+1:] -= np.flip(coefSin[1:offset]) # Include sign swap: sin(-x) = -sin(x)

      # Zero out the zeroth sine coefficients
      coefSin[offset] = 0;

      # Apply the factor 1/2 that appears in the product-to-sum formulas
      return DensityFourier(coefCos[offset:]/2., coefSin[offset:]/2.)
      
     
   def __rmul__(self, other) :
      return self.__mul__(other)



class DensityNormal(object) :
   def __init__(self, mu=0, sigma=1) :
      self.mu    = mu
      self.sigma = sigma

   def clone(self) :
      return DensityNormal(self.mu, self.sigma)

   def normalize(self) :
      pass

   def evaluate(self, phi) :
      # Shift mean by multiples of 2*pi to be closer to interval phi
      mu = self.mu - (2*np.pi) * np.floor((self.mu - phi[0]) / (2*np.pi))
      interval = (np.max(phi) - np.min(phi))
      nWrap = int(np.minimum(100, np.maximum(1,np.ceil(interval / (5*self.sigma)))))
      v = np.zeros(phi.shape,dtype=np.double)
      for k in range(-nWrap,nWrap+1) :
         v += np.exp(-(((phi-self.mu+k*(2*np.pi))/self.sigma)**2) / 2)
      v /= (self.sigma * np.sqrt(2*np.pi))
      return v


   def mean(self) :
      return self.mu

   def variance(self) :
      return (self.sigma**2)

   def std(self) :
      return self.sigma

