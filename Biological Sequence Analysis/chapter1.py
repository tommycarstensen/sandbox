## conditional probability, P(X|Y), prob of X given Y

## joint probability, P(X,Y)=P(X|Y)*P(Y), prob of X and Y

## marginal probability, P(X) = sumY(P(X,Y))

## posterior probability (Bayes' theorem relates conditional probabilities)
## P(X,Y)=P(Y,X)
## P(X|Y)*P(Y)=P(Y|X)*P(X)
## P(X)=P(X|Y)*P(Y)/P(Y|X)

import math
import sys
sys.path.append('../math')
import combinatorics

def main():

    print
    exercise_1_2()

    print
    exercise_1_4()

    print
    exercise_1_5()

    print
    combinatorics()

    print
    probability()

    print
    distribution_geometric()

    print
    distribution_poisson()

    print
    misc()

    return


def misc():

    return


def distribution_poisson():

    ## oligonucleotide occurence in DNA sequences

    ## palindromic sequence
    ## 6 nt
    ## circular sequence
    ## 84k nt
    p = .25**6
    Lambda = p*84000
    print p
    print Lambda
    stop
    f_1_11 = 0
    print '1.11)', f_1_11
    print '1.17)',    
    print '1.19)',    
    print '1.22)',    

    return


def distribution_geometric():

    ## length distribution of restriction fragments
    f_1_12 = 
    print '1.12)', f_1_12

    ## length distribution of ORFs
    print '1.14)',    

    return


def random_variables():

    print '1.13)',    

    print '1.18)',    

    return


def probability():

    ## Bernoulli trial
    freq_error = 0.1
    n = 3

    f_18a = (
        combination(3,3)*(0.9**3)*(0.1**0) ## CCC
        +
        combination(3,2)*(0.9**2)*(0.1**1) ## CCI, CIC, ICC
        )
    print '1.8a)', 'identified and correct (n=3)', f_18a

    print
    f_18b = (
        3*(
            combination(3,3)*((0.1/3)**3)*(1**0) ## I1I1I1
            +
            ## I1I1I2, I1I2I1, I2I1I1, I1I1I3, I1I3I1, I3I1I1
            combination(3,2)*combination(2,1)*((0.1/3)**2)*((0.1/3)**1)
            +
            combination(3,2)*combination(1,1)*((0.1/3)**2)*(0.9**1) ## I1I1C, I1CI1, CI1I1
            )
        )
    print '1.8b)', 'identified but incorrect (n=3)', f_18b
    print '1.8b) P3', combination(3,3)*((0.1/3)**3)*(1**0)
    print '1.8b) P2', combination(3,2)*combination(2,1)*((0.1/3)**2)*((0.1/3)**1),
    print '+', combination(3,2)*combination(1,1)*((0.1/3)**2)*(0.9**1),
    print '=', (
        combination(3,2)*combination(2,1)*((0.1/3)**2)*((0.1/3)**1)
        +
        combination(3,2)*combination(1,1)*((0.1/3)**2)*(0.9**1)
        )
    print '1.8b) 0.00356'

    print
    f_18c = (
        combination(3,1)*permutation(3,3)*(0.9**1)*((0.1/3)**2) ## CI1I2 (3 positions and 1 of 3 nucleotides left out)
        +
        permutation(3,3)*((0.1/3)**3) ## I1I2I3 (permutations)
        )
    print '1.8c)', 'not identified / all different (n=3)', f_18c
    f_18c = 1-f_18a-f_18b
    print '1.8c)', 'not identified / all different (n=3)', f_18c


    print
    f_18d = (
        combination(5,5)*(0.9**5)*(0.1**0) ## CCCCC
        +
        combination(5,4)*(0.9**4)*(0.1**1)
        +
        combination(5,3)*(0.9**3)*(0.1**2)
        )
    print '1.8d)', 'identified and correct (n=5)', f_18d
    f_18d = (
        combination(7,7)*(0.9**7)*(0.1**0) ## CCCCC
        +
        combination(7,6)*(0.9**6)*(0.1**1)
        +
        combination(7,5)*(0.9**5)*(0.1**2)
        +
        combination(7,4)*(0.9**4)*(0.1**3)
        )
    print '1.8d)', 'identified and correct (n=7)', f_18d

    print '1.16)',    

    return


def combinatorics():

    p_C = 0.35
    p_G = 0.35
    p_A = 0.15
    p_T = 0.15
    p_CG = p_C+p_G
    p_AT = p_A+p_T
    p_CG_8_of_15 = (p_CG**8)*(p_AT**7)
    combinations_8_from_15 = combination(15,8)
    p_CG_8_of_15 *= combinations_8_from_15
    print '1.6)', p_CG_8_of_15

    possibilities_mismatch0 = 1
    possibilities_mismatch1 = ((4-1)**1)*combination(8,7)
    possibilities_mismatch2 = ((4-1)**2)*combination(8,6)
    print '1.7)', possibilities_mismatch0+possibilities_mismatch1+possibilities_mismatch2

    print '1.9)'

    f_1_10 = (
        (1./19)*(5*1+9*2+4*3)
        +
        (1./38)*(1*1+1*2) ## Met, Trp
        )/3.
    print '1.10)', f_1_10

    return


def exercise_1_5():

    rolls = [1,3,4,2,4,6,2,1,2,2,]

    p_ML_2 = rolls.count(2)/float(len(rolls))
    print '1.5)', p_ML_2

    for i in range(1,6+1):
        rolls += [i]
    p_ML_2 = rolls.count(2)/float(len(rolls))
    print '1.5)', p_ML_2

    for pseudo_count in range(4):
        for i in range(1,6+1):
            rolls += [i]
    p_ML_2 = rolls.count(2)/float(len(rolls))
    print '1.5)', p_ML_2

    return


def exercise_1_4():

    print

    ## D disease
    ## H healthy

    p_D = 1./1000000
    p_H = 1-p_D
    p_conditional_positive_D = 1.
    p_conditional_positive_H = 0.0001
    p_joint_positive_D = p_conditional_positive_D*p_D
    p_joint_positive_H = p_conditional_positive_H*p_H
    p_marginal_positive = p_joint_positive_D+p_joint_positive_H
    p_conditional_D_positive = p_conditional_positive_D*p_D/p_marginal_positive

    print '1.4)', p_conditional_D_positive

    return


def exercise_1_2():

    ## exercise 1.1) conditional (if), joint (and), marginal probability (sumif)

    p_loaded = 0.01
    p_fair = 1-p_loaded
    p_conditional_six_loaded = 0.50
    p_conditional_six_fair = 1./6.
    p_joint_six_loaded = p_conditional_six_loaded*p_loaded
    p_joint_six_fair = p_conditional_six_fair*p_fair
    p_marginal_six = p_joint_six_loaded + p_joint_six_fair
    print '1.1)', p_marginal_six

    ## exercise 1.2) posterior probability

    ## loop/substitution method
    roll = 1
    while True:
##        p_joint_loaded_six = (.5**roll)*p_loaded + ((1./6.)**roll)*0.99
##        p_six_loaded = 0.5
##        p_loaded_six = (p_six_loaded**roll)*p_loaded/p_six
        p_conditional_sixes_loaded = 0.5**roll
        p_conditional_sixes_fair = (1./6.)**roll
        p_joint_sixes_loaded = p_conditional_sixes_loaded*p_loaded
        p_joint_sixes_fair = p_conditional_sixes_fair*p_fair
        p_marginal_sixes = p_joint_sixes_loaded + p_joint_sixes_fair
        p_conditional_loaded_sixes = p_conditional_sixes_loaded*p_loaded/p_marginal_sixes
        print '1.2)', roll, p_conditional_loaded_sixes
        if p_conditional_loaded_sixes > .5:
            break
        roll += 1
    print '1.2)', roll, p_conditional_loaded_sixes

    ## algebraic method
    p_conditional_loaded_sixes > 0.5
    p_conditional_sixes_loaded*p_loaded/p_marginal_sixes > 0.5
##    (0.5**roll)*p_loaded/(p_joint_sixes_loaded + p_joint_sixes_fair) > 0.5
##    (0.5**roll)*p_loaded/(p_conditional_sixes_loaded*p_loaded + p_conditional_sixes_fair*p_fair) > 0.5
##    (0.5**roll)*p_loaded/((0.5**roll)*p_loaded + ((1./6.)**roll)*p_fair) > 0.5
##    (0.5**roll)*p_loaded/((0.5**roll)*p_loaded + ((1./6.)**roll)*p_fair) > 0.5
##    1/(1 + ((1./6.)**roll)*p_fair / (0.5**roll)*p_loaded ) > 0.5
##    1/(1 + 11*(1/3.)**(roll-2)) > 0.5
##    2 > (1 + 11*(1/3.)**(roll-2))
##    1 > 11*(1/3.)**(roll-2)
##    1/11. > (1/3.)**(roll-2)
##    roll > math.log(1/11.)/math.log(1/3.)+2
    roll = math.log(1/11.)/math.log(1/3.)+2
    
    print '1.2)', roll

    return


def combination(n,r,):

    ## the number of ways to pick in random order r different objects from n different objects

    nCr = combinations = faculty(n)/(faculty(r)*faculty(n-r))

    return combinations


def permutation(n,r,):

    ## counting rule for permutations
    ## the number of ways to arrange in order n different objects within r positions

    nPr = permutations = faculty(n)/faculty(n-r) ## = n*(n-1)*...*(n-r+1)*(n-r)

    return permutations


def faculty(n):

    faculty = 1
    for i in range(1,n+1):
        faculty *= i

    return faculty


if __name__ == '__main__':
    main()
