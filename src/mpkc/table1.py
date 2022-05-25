from mpkc import MQEstimator

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
