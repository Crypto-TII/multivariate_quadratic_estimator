# *****************************************************************************
# Multivariate Quadratic (MQ) Estimator
# Copyright (C) 2021-2022 Emanuele Bellini, Rusydi H. Makarim, Javier Verbel
# Cryptography Research Centre, Technology Innovation Institute LLC
#
# This file is part of MQ Estimator
#
# MQ Estimator is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# MQ Estimator is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# MQ Estimator. If not, see <https://www.gnu.org/licenses/>.
# *****************************************************************************


from mpkc.mq_estimator import MQEstimator
from sage.functions.log import log
from mpkc.utils import ngates, nbits


#Algorithm used (theoretical)
#Run the folowing lines of codes to generate the entries
#corresponding to the column Algorithm used (theoretical)

#Row 1
n = 74; m = 148; q = 2; w = 2.81; theta = 2; h = 10;
E = MQEstimator(n=n, m=m, q=q, w=w, h=h)
time = E.crossbred.time_complexity(k=22,D=4,d=1)
time = ngates(q=q, n=time, theta=theta)
memory = E.crossbred.memory_complexity(k=22,D=4,d=1)
memory = nbits(q=q, n=memory)
print("Row 1 (Type I)")
print(float(log(time,2)), float(log(memory,2)))

#Row 2
n = 99; m = 66; q = 2; w = 2.81; theta=2; h = 0;
E = MQEstimator(n=n, m=m, q=q, w=w, h=h)
time = E.exhaustive_search.time_complexity()
time = ngates(q=q, n=time, theta=theta)
memory = E.exhaustive_search.memory_complexity()
memory = nbits(q=q, n=memory)
print("Row 2 (Type IV)")
print(float(log(time, 2)), float(log(memory, 2)))

#Row 3
n = 103; m = 69; q = 2; w = 2.81; theta=2; h = 17;
E = MQEstimator(n=n, m=m, q=q, w=w, h=h)
time = E.crossbred.time_complexity(k=15, D=4, d=2)
time = ngates(q=q, n=time, theta=theta)
memory = E.crossbred.memory_complexity(k=15, D=4, d=2)
memory = nbits(q=q, n=memory)
print("Row 3 (Type IV)")
print(float(log(time, 2)), float(log(memory, 2)))

#Row 4
n = 37; m = 74; q = 2**8; w = 2.81; theta = 2; h = 0;
E = MQEstimator(n=n, m=m, q=q, w=w, h=h)
time = E.f5.time_complexity()
time = ngates(q=q, n=time, theta=theta)
memory = E.f5.memory_complexity()
memory = nbits(q=q, n=memory)
print("Row 4 (Type II)")
print(float(log(time, 2)), float(log(memory, 2)))

#Row 5
n = 28; m = 19; q = 2**8; w = 2.81; theta = 2; h = 0;
E = MQEstimator(n=n, m=m, q=q, w=w, h=h)
time = E.hybrid_f5.time_complexity(k=1)
time = ngates(q=q, n=time, theta=theta)
memory = E.hybrid_f5 .memory_complexity(k=1)
memory = nbits(q=q, n=memory)
print("Row 5 (Type V)")
print(float(log(time, 2)), float(log(memory, 2)))

#Row 6
n = 38; m = 76; q = 31; w = 2.81; theta = 2; h = 0;
E = MQEstimator(n=n, m=m, q=q, w=w, h=h)
time = E.boolean_solve_fxl.time_complexity()
time = ngates(q=q, n=time, theta=theta)
memory = E.boolean_solve_fxl.memory_complexity()
memory = nbits(q=q, n=memory)
print("Row 6 (Type III)")
print(float(log(time, 2)), float(log(memory, 2)))

#Row 7
n = 30; m = 20; q = 31; w = 2.81; theta=2; h=0;
E = MQEstimator(n=n, m=m, q=q, w=w, h=h)
time = E.hybrid_f5.time_complexity(k=2)
time = ngates(q=q, n=time, theta=theta)
memory = E.hybrid_f5.memory_complexity(k=2)
memory = nbits(q=q, n=memory)
print("Row 7 (Type VI)")
print(float(log(time, 2)), float(log(memory, 2)))

#Table 1
#best algorithms column
w = 2.81
theta = 2

print("Type I")
E = MQEstimator(n=74, m=148, q=2, w=w, theta=theta)
print(E.table(precision=1))

print("Type IV")
E1 = MQEstimator(n=99, m=66, q=2, w=w, theta=theta)
print(E1.table(precision=1))

print("Type IV")
E2 = MQEstimator(n=103, m=69, q=2, w=w, theta=theta)
print(E2.table(precision=1))

print("Type II")
E3 = MQEstimator(n=37, m=74, q=2**8, w=w, theta=theta)
print(E3.table(precision=1))

print("Type V")
E4 = MQEstimator(n=28, m=19, q=2**8, w=w, theta=theta)
print(E4.table(precision=1))

print("Type III")
E5 = MQEstimator(n=38, m=76, q=31, w=w, theta=theta)
print(E5.table(precision=1))

print("Type VI")
E6 = MQEstimator(n=30, m=20, q=31, w=w, theta=theta)
print(E6.table(precision=1))