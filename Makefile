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


%.o : %.c
	gcc -c -o $@ $<

run_experiment_hybrid : run_experiment_hybrid.o bpe.o bpe_plugins.o
	gcc -o $@ $^ -lm

run_experiment_transition : run_experiment_transition.o bpe.o bpe_plugins.o
	gcc -o $@ $^ -lm

verify_clib_fourier : verify_clib_fourier.o bpe.o
	gcc -o $@ $^ -lm

