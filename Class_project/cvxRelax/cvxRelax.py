import numpy as np
import itertools
from sympy.polys.monomials import itermonomials, monomial_count
from sympy import symbols, diff, derive_by_array, Symbol, Mul, simplify
from sympy.polys.orderings import monomial_key
import abc
import IPython
import copy
import os
import math

"""
This file set up the fitting problem using SOS techniques.
"""

# Function to be fitted
# \Psi       = \mathbf{m}^T A \mathbf{m}
# \mathbf{m} = [u v w u^2 v^2 w^2 uv vw wu]
# A \in \mathbb{R}^{9 \times 9}
#
# Hessian
#          h = \mathbf{p}^T C \mathbf{p}
# \mathbf{p} = [x y z ux uy uz vx vy vz wx wy wz]
# A \in \mathbb{R}^{12 \times 12}

# NUM_VAR = 3
# POLY_DEGREE = 4 # Degree of \Psi
# MONO_DEGREE = int(POLY_DEGREE/2)


# uvec = [u, v, w]
class SOS:
    def __init__(self, uvec=None, POLY_DEGREE=None, MONO_LIST=None, SYMBOL_OF_A=None):
        # if POLY_DEGREE != 4 or NUM_VAR != 3 or len(uvec) != 3:
        #     raise ValueError('Assume POLY_DEGREE = 4 and NUM_VAR = 3')
        if POLY_DEGREE== None:
            self.poly_degree=4
        else:
            self.poly_degree=POLY_DEGREE

        if uvec==None:
            self.u = list(symbols('u0:%d'%3))
        else:
            self.u = uvec

        if SYMBOL_OF_A==None:
            self.symbol_of_A='a'
        else:
            self.symbol_of_A=SYMBOL_OF_A
        self.num_of_var = len(self.u)
        self.degree_of_mono = int(self.poly_degree/2)
        self._mono_list = MONO_LIST
        self.initAll()

    def initAll(self):
        self.initMono()
        self.initA()
        self.initPoly()
        self.initG()

    @property
    def mono(self):
        """Monomials vector for Psi"""
        return self._mono

    @property
    def A(self):
        return self._A

    @property
    def poly(self):
        return self._poly
    @property
    def alist(self):
        return self._alist
    @property
    def index_of_aij_by_col_major(self):
        return self._index_of_aij_by_col_major
    
    @property
    def G(self):
        return self._G

    def initMono(self):
        if self._mono_list==None:
            self._mono = sorted(itermonomials(self.u, self.degree_of_mono), \
                key=monomial_key('grlex', self.u))
        elif self._mono_list!=None:
            self._mono = self._mono_list

    def initA(self):
        self.dim_A = len(self._mono)
        self._num_of_aij_by_col_major = [self.dim_A-i for i in range(self.dim_A)]
        temlist = \
            [[(i+j, j) for i in range(self.dim_A-j)] for j in range(self.dim_A)]
        self._index_of_aij_by_col_major=[]
        for l in temlist:
            for item in l:
                self._index_of_aij_by_col_major.append(item)
        self._num_of_aij = sum(self._num_of_aij_by_col_major)
        self._alist = np.array(symbols(self.symbol_of_A+'0:%d'%self._num_of_aij))
        self._A = self.convert_alist_to_A(self._alist)

    def convert_alist_to_A(self, alist):
        temdic={}
        for n in range(self._num_of_aij):
            i, j = self._index_of_aij_by_col_major[n]
            temdic[(i,j)] = alist[n]
            if i>j:
                temdic[(j,i)] = alist[n]
        
        teml=[]
        for i in range(self.dim_A):
            temll = []
            for j in range(self.dim_A):
                temll.append(temdic[(i,j)])
            teml.append(temll)
        A = np.array(teml)
        return A

    def initPoly(self):
        self._poly = np.dot(self._mono, np.dot(self._A, self._mono))
    
    def initG(self):
        # f and its gradient g is linear in a.
        # g = G * a
        x = self.u
        a = self.alist
        f = self.poly
        g = [diff(f, x[i]) for i in range(len(x))]
        self._G = [[diff(g[i], a[j]) for j in range(len(a))] for i in range(len(x))]


def derive_SOS_mat(psi, h_flat, h_sos):
    print("Start calculating coef matching matrix")
    m = len(psi.alist)
    n = len(h_sos.alist)
    h_sos_flat = h_sos.poly
    temdict = {}
    teml = []
    for i in range(n):
        term = diff(h_sos_flat, h_sos.alist[i])
        temll = []
        for j in range(m):
            temterm = term
            freesymbols = list(temterm.free_symbols)
            temexpr = h_flat
            temvarr = None
            count = None
            while bool(temterm.free_symbols):
                # if power of temvar greater than one, coef * factorial.
                freesymbols = list(temterm.free_symbols)
                temvar = freesymbols[0]
                if temvar != temvarr:
                    count = 1
                    temvarr = temvar
                else:
                    count += 1
                temexpr = diff(temexpr, temvar)
                temterm = diff(temterm, temvar)
                temexpr /= count
            temexpr = simplify(temexpr)
            
            temexpr2 = diff(temexpr, psi.alist[j])
            
            temexpr2_fs = list(temexpr2.free_symbols)
            ins = {ind:0 for ind in temexpr2_fs}
            temexpr2 = temexpr2.subs(ins)

            

            k, l = h_sos.index_of_aij_by_col_major[i]
            # if k != l:
            #     temexpr2 = temexpr2/2
            temdict[i, j] = temexpr2
            temll.append(temexpr2)
        teml.append(temll)
        print('finished row '+str(i))
    print('Finished extract coeff matching matrix')
    temD = np.array(teml)

    temcount = count_mul_terms(h_sos)
    
    for i in range(len(h_sos.alist)):
        k, l = h_sos.index_of_aij_by_col_major[i]
        temD[i] /= temcount[(k, l)]

    return temD


def count_mul_terms(h_sos):
    termmat = np.outer(h_sos.mono, h_sos.mono)
    termlist = list(termmat.flatten())

    mm, nn = termmat.shape
    temcount = {}
    for i in range(mm):
        for j in range(nn):
            temcount[(i,j)] = termlist.count(termmat[i][j])
    return temcount



psi = SOS(SYMBOL_OF_A='a')

u = psi.u
f = psi.poly
g = derive_by_array(f, u)
H = [derive_by_array(g[i], u) for i in range(len(u))]

nvar = len(u)
x = list(symbols('x0:%d'%nvar))
xvec = np.array(x)
Hmat = np.array(H)
h_flat = simplify(np.dot(xvec, np.dot(Hmat, xvec)))


#TODO: feed SOS with a monomial vector
p = list(x) + list(xvec*u[0]) + list(xvec*u[1]) + list(xvec*u[2])

h_sos = SOS(uvec=list(x)+list(u), \
    POLY_DEGREE=4, MONO_LIST=p, SYMBOL_OF_A='c')

fname = 'D.npy'
if os.path.isfile(fname):
    D = np.load(fname)
else:
    D = derive_SOS_mat(psi, h_flat, h_sos)
    np.save(fname, D)

# verify D
clist = np.dot(D, psi.alist)
C = h_sos.convert_alist_to_A(clist)
h_induced = simplify(np.dot(p, np.dot(C, p)))
err1 =  simplify(h_induced - h_flat)

# verify G
G = psi.G
err2 = [simplify(g[i] - np.dot(G[i], psi.alist)) \
    for i in range(len(u))]