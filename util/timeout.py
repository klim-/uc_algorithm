from functools import wraps
import errno
import os
import signal
from time import sleep

import sympy as sp
from sympy import sin, cos

import symbtools as st



# This snippet is strongly inspired by
# http://stackoverflow.com/a/2282656/333403
class TimeoutError(Exception):
    pass

def timeout(seconds=10):
    
    def decorator(func):
        
        def _handle_timeout(signum, frame):
            """
            this function is called if time is over
            """
            
            raise TimeoutError("Scheduled time is over.")

        def wrapper(*args, **kwargs):
            """
            Wraps the original function such that the proper exception is
            raised if time is over.
            """
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wraps(func)(wrapper)

    return decorator 


def get_timed_simplify(t):
    
    # create decorator
    timeout_decorator = timeout(t)
    
    # apply decorator
    timeout_exception_simplify = timeout_decorator(sp.simplify)
    
    def wrapper(expr, **kwargs):
        
        try:
            result = timeout_exception_simplify(expr, **kwargs)
        except TimeoutError:
            result = expr
            print("simplification stopped for expr.count_ops:")
            print(st.count_ops(expr))
        return result
        
    return wrapper
    


if __name__ == "__main__":
    x, y = sp.symbols('x, y')

    expr = sp.expand( (x*sp.cos(x)**2 + x*sp.sin(x)**2)**(sp.cos(y)**2 + sp.sin(y)**2))

    expr = (x**2 + x)/(x*sin(y)**2 + x*cos(y)**2)

    simplify3 = get_timed_simplify(3)
    simplify1 = get_timed_simplify(1)

    print(expr)
    print(sp.simplify(expr))
    print(simplify3(expr))
    print(simplify1(expr))


    
if 0:    
    @timeout(3)
    def long_running_function(n):
        for i in range(n):
            sleep(.5)
            print(i)

            
    long_running_function(8)
