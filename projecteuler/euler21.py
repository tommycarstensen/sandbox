import math
set_amicable = set()
for n in range(1,10000):
    if n in set_amicable:
        continue
    d1 = 1+sum(i+n/i for i in range(2,int(math.sqrt(n))) if not n%i)
    if d1 <= n:
        continue
    d2 = 1+sum(i+d1/i for i in range(2,int(math.sqrt(d1))) if not d1%i)
    if n == d2:
        set_amicable |= set([n,d1])
print(sum(set_amicable))
