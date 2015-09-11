# set up system
s = sp.symbols("s")

x1, x2, x3, x4 = \
                    sp.symbols("x1, x2, x3, x4")
xdot1, xdot2, xdot3, xdot4 = \
                    sp.symbols("xdot1, xdot2, xdot3, xdot4")

dx1, dx2, dx3, dx4 = \
                    sp.symbols("dx1, dx2, dx3, dx4")
dxdot1, dxdot2, dxdot3, dxdot4 = \
                    sp.symbols("dxdot1, dxdot2, dxdot3, dxdot4")

c1, d1, m2 = sp.symbols("c1, d1, m2")

#vec_dxdot = sp.Matrix([dxdot1, dxdot2, dxdot3, dxdot4])
#vec_dx = sp.Matrix([dx1, dx2, dx3, dx4])
vec_x = sp.Matrix([x1, x2, x3, x4]) # muss in jedem schritt um die ableitungen ergänzt werden

P1i = sp.Matrix([ [1,0,0,0],[0,1,0,0],[0,0,0,1] ])
P0i = sp.Matrix([ [0,0,-1,0],[0,0,0,-1],[c1/m2,-c1/m2,d1/m2,-d1/m2] ])



