# -*- coding: utf-8 -*-

import sympy as sp
import symb_tools as st
import non_commutative_tools as nct

from IPython import embed

# Symbole fuer y
yy = sp.Matrix(sp.symbols('y1, y2', commutative=False))
s = sp.Symbol('s', commutative=False)
yyd1 = st.perform_time_derivative(yy, yy, order=1, commutative=False)
yyd2 = st.perform_time_derivative(yy, yy, order=2, commutative=False)
yyd3 = st.perform_time_derivative(yy, yy, order=3, commutative=False)
yyd4 = st.perform_time_derivative(yy, yy, order=4, commutative=False)

#Hilfssymbole zur Konstruktion von Gleichungen
a1, a2 = aa = sp.symbols('a1:3')

# Funktionen der Zeit
t = nct.t
yyf = st.symbs_to_func(yy, yy, nct.t)
yyfd1 = yyf.diff(t, 1)
yyfd2 = yyf.diff(t, 2)
yyfd3 = yyf.diff(t, 3)
yyfd4 = yyf.diff(t, 4)
yyaf = st.row_stack(yyf, yyfd1, yyfd2, yyfd3, yyfd4)

yya = st.row_stack(yy, yyd1, yyd2, yyd3, yyd4)

st.makeGlobal(yya, 1)

# einfacheres Beispiel:
# ich konstruiere mir zweistufig eine integrallos umkehrbare Trafo
# 1. Stufe
z1 = y1*ydot2
z2 = y2

ZZ = ZZ1 = sp.Matrix([z1, z2])

# 2. Stufe
xi1 = z1
xi2 = z2 * st.perform_time_derivative(z1, yya, commutative=False)

# damit wird weitergerechnet
ZZ = ZZ2 = sp.Matrix([xi1, xi2])

# Umkehrung (Konsistenz-Probe):

Z1 = xi1
Z2 = xi2/st.perform_time_derivative(xi1, yya, commutative=False)

Y1 = Z1/st.perform_time_derivative(Z2, yya, commutative=False)
Y2 = Z2

YY = sp.Matrix([Y1, Y2])
assert YY == yy


yya.reshape(5, 2)

def lie_baeck_jacobi_matrix(ZZ, yya, ny):
    
    # ny is the number of scalar components of y
    
    # determine the number of derivatives
    N = len(yya)*1.0/ny
    assert int(N) == N
    #print N
    
    P = sp.zeros(ny, ny)
    # yya is a col vector  -> convert to matrix
    # each row contains one derivative oder
    yyam = yya.reshape(int(N), ny)
    
    # now convert to list of rows    
    for i, vv in enumerate( yyam.tolist() ):
        #print vv, i
        P+=ZZ.jacobian(vv) *s**i
    return P

# diese matrix sollte unimodular sein. d.h. sie sollte eine polynommatirx
# in s als inverse haben --> mit der determinante kann man das scheinbar
# nicht zeigen
P = lie_baeck_jacobi_matrix(ZZ, yya, 2)

J1 = lie_baeck_jacobi_matrix(ZZ1, yya, 2)
J2 = lie_baeck_jacobi_matrix(ZZ2, yya, 2)

# Jetzt versuche ich eine Inverse für P zu konstruieren
# bzw. erstmal für J2

I0 = st.symbMatrix(2, 2, 'a', commutative=False)
I1 = st.symbMatrix(2, 2, 'b', commutative=False)
I2 = st.symbMatrix(2, 2, 'c', commutative=False)
st.makeGlobal(I0, 1)
st.makeGlobal(I1, 1)
st.makeGlobal(I2, 1)


M = I0 + I1*s + I2*s**2

# Gleichung die die Inverse definiert
# nc_mul: non commutative multiplication
R = nct.nc_mul(M, J2).subs(zip(yya, yyaf)) - sp.eye(M.shape[0])
eq = list(R)

tmp = M*J2 - nct.nc_mul(M, J2)
tmp.simplify()
assert tmp == tmp*0

#nct.right_shift_all(eq[0].expand(), s, t)
K = sp.Symbol('K')

# s überall nach rechts schieben
eq2 = [nct.right_shift_all(e, s, t) for e in eq]

back_subs = list(reversed(zip(yyaf, yya)))

eq1_symb = [e.subs(back_subs) for e in eq]
eq2_test = sp.Matrix([nct.right_shift_all(e, s, t, yya) for e in eq1_symb])

# Funktionen wieder durch Symbole ersetzen
eq2 = sp.Matrix(eq2).subs(back_subs)

w = eq2-eq2_test
w.simplify() # sollte den nullvektor ergeben


# Gleichungen für die einzelnen s-Potenzen aufstellen
eq3 = st.row_stack(eq2.subs(s,0), eq2.diff(s).subs(s,0), eq2.diff(s, 2).subs(s,0), eq2.diff(s, 3).subs(s,0) )

# jetzt spielt Nicht-Kommutativität keine Rolle mehr
eq3, st_c_nc = nct.make_all_symbols_commutative(eq3)

# Parameter alle in einem Vektor zusammenfassen
pars = sp.Matrix(list(I0) + list(I1) + list(I2))

pars_c, st_cnc_pars = nct.make_all_symbols_commutative(pars)


JE = eq3.jacobian(pars_c)

# Homogene gleichungen rausfinden
eq3_0 = eq3.subs(st.zip0(pars_c))

inhom_idcs = st.np.where(st.to_np(eq3_0) != 0)[0]
hom_idcs = st.np.where(st.to_np(eq3_0) == 0)[0]

# Homogene und Inhomogene Gleichungen trennen
eq4 = sp.Matrix(st.np.array(eq3)[hom_idcs])
eq5 = sp.Matrix(st.np.array(eq3)[inhom_idcs])

# Probe
assert eq4.subs(st.zip0(pars_c)) == eq4*0

JEh = eq4.jacobian(pars_c)

JEh = JEh.expand()

kk = st.nullspaceMatrix(JEh)

k1, k2 = st.col_split(kk)

K = a1*k1 + a2*k2
eq6 = eq5.subs(zip(pars_c, K))
eq6.simplify()

sol = sp.solve(eq6, aa)

tmp = eq6.subs(sol)
tmp.simplify()
assert tmp == eq6*0

KK = K.subs(sol)
KK.simplify()

# KK enthält die Lösung für den Parameter Vektor (kommutative Variablen)
# Jetzt fehlt noch dir Rücksubstitution auf nichtkommutative Variablen

KK_nc = KK.subs(st_c_nc)



# PROBE OB TATS. DIE EINHEITSMATRIX RAUSKOMMT
sol_c = zip(pars_c, KK)
sol_nc = zip(pars, KK_nc)
Inv = M.subs(sol_nc)
Invf = Inv.subs(zip(yya, yyaf))
J2f = J2.subs(zip(yya, yyaf))

rr = nct.nc_mul(J2f, Invf)
rr2 = nct.right_shift_all(rr)

test1 = st.subs_random_numbers(rr2) # sollte einheitsmatrix sein


ww = nct.nc_mul(Invf, J2f)
ww2 = nct.right_shift_all(ww)

test2 = st.subs_random_numbers(ww2) # sollte ebenfalls einheitsmatrix sein

r = nct.nc_mul(J2, Inv)
r2 = nct.right_shift_all(r, s, t, yy)
test3 = st.subs_random_numbers(r2) # sollte auch einheitsmatrix sein

embed()
