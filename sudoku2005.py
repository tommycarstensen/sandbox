#!/bin/env python
#
# $Id$
#
# (C) University College Dublin 2005
#
class sudoku:

    def main(self):
        sudoku=[[0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0],
                [0,0,0,0,0,0,0,0,0]]
        sudoku=[[7,1,0,0,9,0,0,0,5],
                [3,6,5,0,2,0,9,0,1],
                [2,0,0,1,0,5,0,7,6],
                [1,4,3,0,0,2,5,9,8],
                [5,0,0,9,4,1,0,6,3],
                [9,2,6,5,0,0,4,1,7],
                [4,3,1,2,5,0,0,0,9],
                [8,0,2,0,0,0,1,5,4],
                [6,5,0,0,1,0,0,3,2]]

## len(n) = 1 er det ensomme tal

##        import Numeric
##        m = Numeric.zeros((9,9))
        m = [[[] for j in range(9)] for i in range(9)]
        for row in range(9):
            for column in range(9):
                m[row][column] = sudoku[row][column]
##        n = Numeric.zeros((9,9))
        n = [[[] for j in range(9)] for i in range(9)]
        set_all = set([0,1,2,3,4,5,6,7,8,9])
        n = [[set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all]]
        o = [[set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all],
             [set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all,set_all]]

        count = 0
        m_sum = sum([sum(m[i]) for i in range(9)])
##        while sum(sum(m)) < 9*(9*(1+9)/2):
        while m_sum < 9*(9*(1+9)/2):
            count += 1
            if count > 100:
                print m
                print
                for row in range(9):
                    print n[row]
                print
                for row in range(9):
                    print o[row]
                print count
                stop
            
##            print m
##            print n
##            n(elements of block, row, column).discard(value if determined right after determination)
##            n_element = o_element (if len(o_element) != 0 and len(o_element) < len(n_element))
            for row_global in range(3):
                for column_global in range(3):
                    for row_local in range(3):
                        for column_local in range(3):

                            if m[3*row_global+row_local][3*column_global+column_local] == 0:
                                set_block_m = self.f_set_block(m, n, row_global, column_global)
                                set_row_m = self.f_set_row(m, row_global, row_local)
                                set_column_m = self.f_set_column(m, column_global, column_local)
                            
                                set_field_m = set_block_m | set_row_m | set_column_m
                                set_field_n = set_all-set_field_m
                                set_field_o = self.f_set_field(n, row_global, row_local, column_global, column_local, set_field_n)
                                if len(set_field_m) == 9:
                                    m[3*row_global+row_local][3*column_global+column_local] = list(set_all-set_field_m)[0]
                                    n[3*row_global+row_local][3*column_global+column_local] = set([m[3*row_global+row_local][3*column_global+column_local]])
                                    o[3*row_global+row_local][3*column_global+column_local] = set([m[3*row_global+row_local][3*column_global+column_local]])
                                elif len(set_field_n) == 1:
                                    m[3*row_global+row_local][3*column_global+column_local] = list(set_field_n)[0]
                                    n[3*row_global+row_local][3*column_global+column_local] = set([m[3*row_global+row_local][3*column_global+column_local]])
                                    o[3*row_global+row_local][3*column_global+column_local] = set([m[3*row_global+row_local][3*column_global+column_local]])
                                elif len(set_field_o) == 1:
                                    m[3*row_global+row_local][3*column_global+column_local] = list(set_field_o)[0]
                                    n[3*row_global+row_local][3*column_global+column_local] = set([m[3*row_global+row_local][3*column_global+column_local]])
                                    o[3*row_global+row_local][3*column_global+column_local] = set([m[3*row_global+row_local][3*column_global+column_local]])
                                else:
                                    n[3*row_global+row_local][3*column_global+column_local] = set_field_n
                                    o[3*row_global+row_local][3*column_global+column_local] = set_field_o
                            elif m[3*row_global+row_local][3*column_global+column_local] != 0 and len(n[3*row_global+row_local][3*column_global+column_local]) != 1:
                                set_field_n = set([m[3*row_global+row_local][3*column_global+column_local]])
                                n[3*row_global+row_local][3*column_global+column_local] = set_field_n
                                set_field_o = self.f_set_field(n, row_global, row_local, column_global, column_local, set_field_n)
                                o[3*row_global+row_local][3*column_global+column_local] = set_field_o

        print ' _______ _______ _______'
        print '|       |       |       |'
        for i in range(9):
            print '|', m[i][0], m[i][1], m[i][2], '|', m[i][3], m[i][4], m[i][5], '|', m[i][6], m[i][7], m[i][8], '|'
            if i in [2,5]:
                print '|_______|_______|_______|'
                print '|       |       |       |'
        print '|_______|_______|_______|'
##        print m

    def f_set_field(self, n, row_global, row_local, column_global, column_local, set_field_n):
        set_field_o = set()
        for row in range(3):
            for column in range(3):
                if row == row_local and column == column_local:
                    continue
                set_field_o = set_field_o | n[3*row_global+row][3*column_global+column]
        return set_field_n - set_field_o

    def f_set_block(self, m, n, row_global, column_global):
        set_block_m = set([0])
        for row_local in range(3):
            for column_local in range(3):
                set_block_m.add(m[3*row_global+row_local][3*column_global+column_local])
        return set_block_m

    def f_set_row(self, m, row_global, row_local):
        set_row_m = set([0])
        for i in range(9):
            set_row_m.add(m[3*row_global+row_local][i])
        return set_row_m

    def f_set_column(self, m, column_global, column_local):
        set_column_m = set()
        for i in range(9):
            set_column_m.add(m[i][3*column_global+column_local])
        return set_column_m

if __name__=='__main__':
    instance_sudoku = sudoku()
    instance_sudoku.main()
