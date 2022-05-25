from mpkc import MQEstimator
from sage.functions.log import log

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


#Algorithm used (theoretical)
#Run the folowing lines of codes to generate the entries
#corresponding to the column Algorithm used (theoretical)


from mpkc import MQEstimator
E = MQEstimator(n=74, m=148, q=2, w=2.81, theta=2, h=10)
time = E.crossbred.time_complexity(k=22,D=4,d=1)
memory = E.crossbred.memory_complexity(k=22,D=4,d=1)
print("Row 1 (Type I)")
print(float(log(time,2)), float(log(memory,2)))

E = MQEstimator(n=99, m=66, q=2, w=2.81, theta=2, h=0)
time = E.exhaustive_search.time_complexity()
memory = E.exhaustive_search.memory_complexity()
print("Row 2 (Type IV)")
print(float(log(time, 2)), float(log(memory, 2)))

E = MQEstimator(n=103, m=69, q=2, w=2.81, theta=2, h=17)
time = E.crossbred.time_complexity(k=15, D=4, d=1)
memory = E.crossbred.memory_complexity(k=15, D=4, d=1)
print("Row 3 (Type IV)")
print(float(log(time, 2)), float(log(memory, 2)))

E = MQEstimator(n=37, m=74, q=2**8, w=2.81, theta=2, h=0)
time = E.f5.time_complexity()
memory = E.f5.memory_complexity()
print("Row 4 (Type II)")
print(float(log(time, 2)), float(log(memory, 2)))

E = MQEstimator(n=28, m=19, q=2**8, w=2.81, theta=2, h=1)
time = E.hybrid_f5.time_complexity()
memory = E.hybrid_f5 .memory_complexity()
print("Row 5 (Type V)")
print(float(log(time, 2)), float(log(memory, 2)))

E = MQEstimator(n=38, m=76, q=31, w=2.81, theta=2, h=0)
time = E.boolean_solve_fxl.time_complexity()
memory = E.boolean_solve_fxl.memory_complexity()
print("Row 6 (Type III)")
print(float(log(time, 2)), float(log(memory, 2)))

E = MQEstimator(n=30, m=20, q=31, w=2.81, theta=2, h=2)
time = E.hybrid_f5.time_complexity()
memory = E.hybrid_f5.memory_complexity()
print("Row 6 (Type III)")
print(float(log(time, 2)), float(log(memory, 2)))
