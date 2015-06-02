# -*- coding: utf-8 -*-
import unittest

import sympy as sp
from sympy import sin, cos, exp

import numpy as np

import unimodular as um

import symb_tools as st
from IPython import embed as IPS


class InteractiveConvenienceTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_is_unit_matrix_1(self):
        matrix_1 = sp.Matrix([
        [1.0, 5.63091748710959e-18],
        [  0,                  1.0]])

        res_1 = um.is_unit_matrix(matrix_1)
        self.assertEqual(res_1, True)
        
    def test_is_unit_matrix_2(self):
        matrix_2 = sp.Matrix([
        [   0.999999999999996, 7.105427357601e-15],
        [1.01389201174282e-13,  0.999999999999824]])

        res_2 = um.is_unit_matrix(matrix_2)
        self.assertEqual(res_2, True)

    def test_is_unit_matrix_3(self):
        matrix_3 = sp.Matrix([
        [                 1.0,   0],
        [1.04333429557543e-17, 1.0]])
        
        res_3 = um.is_unit_matrix(matrix_3)
        self.assertEqual(res_3, True)

    def test_is_unit_matrix_4(self):
        matrix_4 = sp.Matrix([
        [1.01389201174282e-13,   7.105427357601e-15],
        [1.04333429557543e-17, 5.63091748710959e-18]])
        
        res_4 = um.is_unit_matrix(matrix_4)
        self.assertEqual(res_4, False)


def main():
    unittest.main()

if __name__ == '__main__':
    main()
