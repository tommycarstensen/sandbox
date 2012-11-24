import random

bool = False
while bool == False:
    no = int(10000000000000000*random.random())
    if int(str(no)[0]) not in [4,5,]:
        continue
    if len(str(no)) != 16:
        continue
    Sum = 0
    for i in range(0,16,2):
        digit = int(str(no)[i])
        double = str(2*digit)
        for j in range(len(double)):
            Sum += int(double[j])
    for i in range(1,16,2):
        digit = int(str(no)[i])
        Sum += digit
    print Sum
    if Sum%10 == 0:
        bool = True
        print no
