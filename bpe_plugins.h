#ifndef __BPE_PLUGINS_H__
#define __BPE_PLUGINS_H__

#include "bpe.h"


/* ======================================================================== */
/* Ground truth functions                                                   */
/* ======================================================================== */
void pluginCreateGroundTruth1(BPEPhases **phases, int idx);
void pluginCreateGroundTruth2(BPEPhases **phases, int idx);
void pluginCreateGroundTruth3(BPEPhases **phases, int idx);
void pluginCreateGroundTruth4(BPEPhases **phases, int idx);

void pluginCreateGroundTruthMultiple1(BPEPhases **phases, int idx); 
void pluginCreateGroundTruthMultiple2(BPEPhases **phases, int idx); 
void pluginCreateGroundTruthMultiple3(BPEPhases **phases, int idx); 
void pluginCreateGroundTruthMultiple4(BPEPhases **phases, int idx); 
void pluginCreateGroundTruthMultiple5(BPEPhases **phases, int idx); 
void pluginCreateGroundTruthMultiple6(BPEPhases **phases, int idx); 
void pluginCreateGroundTruthMultiple7(BPEPhases **phases, int idx); 
void pluginCreateGroundTruthMultiple8(BPEPhases **phases, int idx); 
void pluginCreateGroundTruthMultiple9(BPEPhases **phases, int idx); 
void pluginCreateGroundTruthMultiple10(BPEPhases **phases, int idx);
void pluginCreateGroundTruthMultiple11(BPEPhases **phases, int idx);
void pluginCreateGroundTruthMultiple12(BPEPhases **phases, int idx);

void pluginCreateGroundTruthFixed7(BPEPhases **phases, int idx);
void pluginCreateGroundTruthFixed10(BPEPhases **phases, int idx);
void pluginCreateGroundTruthFixed12(BPEPhases **phases, int idx);

void pluginCreateNoisyGroundTruth3_10_2(BPEPhases **phases, int idx);
void pluginCreateNoisyGroundTruth3_10_5(BPEPhases **phases, int idx);
void pluginCreateNoisyGroundTruth3_10_10(BPEPhases **phases, int idx);
void pluginCreateNoisyGroundTruth3_10_20(BPEPhases **phases, int idx);


/* ======================================================================== */
/* Measurement functions                                                    */
/* ======================================================================== */

/* Plugin measurement function for standard noiseless rounds */
void pluginMeasureRounds(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);

/* Plugin measurement function for single-round measurement with decoherence noise */
void pluginMeasureDecoherence(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);

/* Plugin measurement function for Fourier */
void pluginMeasureFourier(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);

/* Plugin measurement function for noisy measurements */
void pluginMeasureNoisy(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);


/* ======================================================================== */
/* Selection of k values                                                    */
/* ======================================================================== */

/* Plugin for random beta values */
void pluginRandomBeta(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);

/* Selection of k values */
void pluginAdaptiveK(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginAdaptiveCyclicK(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginCyclicFixedK(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginOLDCyclicK(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginAdaptiveSampleK(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);

void pluginAlternateK2(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginAlternateK3(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginAlternateK5(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginAlternateK10(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginAlternateK20(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginAlternateK50(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginAlternateK80(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginAlternateK100(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginAlternateK256(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginAlternateK512(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginAlternateK1024(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);
void pluginAlternateK4096(BPEMeasurementScheme *scheme, BPEDensity *density, int iteration);


/* ======================================================================== */
/* Measurement schemes                                                      */
/* ======================================================================== */

/* Measurement scheme A: {k=[1,2,5], plugingSetRandomBeta, pluginMeasureRounds} */
void pluginMeasurementSchemeA(BPEMeasurementScheme **scheme);

/* Measurement scheme B: {Adaptive k, plugingSetRandomBeta, pluginMeasureRounds} */
void pluginMeasurementSchemeB(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB2(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB3(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB5(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB10(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB20(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB50(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB256(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB512(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB1024(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB4096(BPEMeasurementScheme **scheme);

/* Measurement scheme C: Purely adaptive, mean k value */
void pluginMeasurementSchemeC(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeC4096(BPEMeasurementScheme **scheme);

/* Measurement scheme D: Cyclic adaptive, mean k value */
void pluginMeasurementSchemeD(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeD4096(BPEMeasurementScheme **scheme);

/* Noise schemes */
void pluginMeasurementSchemeB10_Noisy(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB20_Noisy(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB50_Noisy(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB80_Noisy(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB100_Noisy(BPEMeasurementScheme **scheme);

void pluginMeasurementSchemeB10_Noisy_25(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB20_Noisy_25(BPEMeasurementScheme **scheme);

void pluginMeasurementSchemeC10_Noisy(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeC20_Noisy(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeC50_Noisy(BPEMeasurementScheme **scheme);

void pluginMeasurementSchemeD10_Noisy(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeD20_Noisy(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeD50_Noisy(BPEMeasurementScheme **scheme);

void pluginMeasurementSchemeE(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeE_Noisy_2(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeE_Noisy_5(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeE_Noisy_8(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeE_Noisy_10(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeE_Noisy_15(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeE_Noisy_20(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeE_Noisy_25(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeE_Noisy_26(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeE_Noisy_30(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeE_Noisy_35(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeE_Noisy_40(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeE_Noisy_45(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeE_Noisy_48(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeE_Noisy_50(BPEMeasurementScheme **scheme);

void pluginMeasurementSchemeFixed2(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeFixed3(BPEMeasurementScheme **scheme);

void pluginMeasurementSchemeX(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeXInf(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeY(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeYInf(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeB_Noisy(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeX_Noisy(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeX10_Noisy(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeX20_Noisy(BPEMeasurementScheme **scheme);

void pluginMeasurementSchemeX_Noisy30(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeX_Noisy40(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeX_Noisy50(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeX_Noisy60(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeY_Noisy(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeF_Noisy(BPEMeasurementScheme **scheme);
void pluginMeasurementSchemeF_Noisy0001(BPEMeasurementScheme **scheme);


/* ======================================================================== */
/* Initial densities                                                        */
/* ======================================================================== */

typedef struct
{  BPEInitDensity base;
   double initSigma;
   double initBias;
} BPEInitDensityB3;

BPEInitDensity *pluginBPEInitDensityB3(double sigma, double bias);
void pluginBPEInitDensityB3_InitDensity(BPEInitDensity *self, BPEDensity *);
void pluginInitWeightsB3(BPEConfig *config, double *weights);


typedef struct
{  BPEInitDensity base;
   double initSigma;
   double initBias;
} BPEInitDensityA;

BPEInitDensity *pluginBPEInitDensityA(double sigma, double bias);
void pluginBPEInitDensityA_InitDensity(BPEInitDensity *self, BPEDensity *);
void pluginInitWeightsA(BPEConfig *config, double *weights);

#endif
