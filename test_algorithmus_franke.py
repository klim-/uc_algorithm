# -*- coding: utf-8 -*-

import unittest

import sympy as sp
from sympy import sin, cos, exp

import algorithmus_franke as af
from IPython import embed as IPS


class SymbToolsTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_srank(self):
        sp.var("a,b,c,d,e,f")
        M = sp.Matrix([
        [a-a,b,c,d],
        [0,b,c,c],
        [0,0,0,d]])

        self.assertTrue( M.rank(), af.srank(M) )

    def test_has_left_ortho_complement(self):
        pass

    def test_has_right_ortho_complement(self):
        pass

    def test_is_zero_matrix(self):
        M = sp.Matrix([
        [0.000000001,0,0,0,0],
        [0,0,0,0.00000001,0],
        [0,0,0,0,0]])

        self.assertTrue( af.is_zero_matrix(M) )

    def test_is_unit_matrix(self):
        M = sp.Matrix([
        [1.0000000, 0.000000001,0,0,0],
        [0,1,0,0.00000001,0],
        [0,0,1.00000001,0.00000001,0],
        [0,0,0,1.00000001,0],
        [0,0,0,0,1]])

        self.assertTrue( af.is_unit_matrix(M) )


def main():
    unittest.main()

if __name__ == '__main__':
    main()
