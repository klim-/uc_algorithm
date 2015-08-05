# -*- coding: utf-8 -*-
import unittest

import sympy as sp
from sympy import sin, cos, exp

import numpy as np

import algebra as al

import symb_tools as st
from IPython import embed as IPS

"""
angelegt von Carsten, um Fehler zu finden und zukünftig zu vermeiden
"""


class InteractiveConvenienceTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_is_linearly_independent(self):
        
        
        A = st.symbMatrix(6,8)
        A[:-2,:]*=0
        res = al.is_linearly_independent(A[:, :2], A[:, 3])
        
        
        self.assertFalse(res)
        
def main():
    unittest.main()

if __name__ == '__main__':
    main()
