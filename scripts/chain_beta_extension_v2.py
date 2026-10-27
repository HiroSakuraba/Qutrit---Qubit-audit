#!/usr/bin/env python3
"""Chain Theory beta-function extension: reproducible, conditional qudit gauge/Dirac checks.

This is an EXPLICIT ADDITION to the paper's axioms, not a derivation from them.
Requires numpy and sympy. Run: python chain_beta_extension_v2.py
"""
from __future__ import annotations
import itertools
import math
import numpy as np
import sympy as sp


def matrices_two_qubit_clifford():
    sx=sp.Matrix([[0,1],[1,0]])
    sy=sp.Matrix([[0,-sp.I],[sp.I,0]])
    sz=sp.diag(1,-1)
    I=sp.eye(2)
    gam=[sp.kronecker_product(sx,s) for s in (sx,sy,sz)]
    gam.append(sp.kronecker_product(sy,I))
    eye=sp.eye(4)
    for a in range(4):
        for b in range(4):
            assert gam[a]*gam[b]+gam[b]*gam[a]==(2 if a==b else 0)*eye
    return gam


def exact_heat_coefficient(field,charge=1):
    """Return spin/orbital heat-kernel coefficients for Dslash=gamma_mu(partial+iqA).
    Local a4 density excludes the universal factor (4*pi)^(-2).
    P=-Dslash^2=-(covariant_laplacian+E), Omega=i q F I_4,
    E=(i q/4)[gamma_mu,gamma_nu] F_mu_nu.
    """
    gam=matrices_two_qubit_clifford()
    assert field.shape==(4,4) and field == -field.T
    E=sp.zeros(4)
    omega_sq=sp.zeros(4)
    for a in range(4):
        for b in range(4):
            E += sp.I*charge*sp.Rational(1,4)*(gam[a]*gam[b]-gam[b]*gam[a])*field[a,b]
            omega_sq += -charge*charge*field[a,b]*field[a,b]*sp.eye(4)
    spin=sp.simplify(sp.trace(E*E)/2)
    orbital=sp.simplify(sp.trace(omega_sq)/12)
    F2=sum(field[a,b]**2 for a in range(4) for b in range(4))
    return orbital,spin,sp.simplify((spin+orbital)/F2),sp.simplify(F2)


def test_gamma_and_heats():
    fields=[]
    F=sp.zeros(4);F[0,1]=1;F[1,0]=-1
    fields.append(F)
    F=sp.zeros(4)
    for a,b,k in [(0,1,3),(0,2,-2),(1,3,5),(2,3,7)]:
        F[a,b]=k;F[b,a]=-k
    fields.append(F)
    out=[]
    for F in fields:
        o,s,rat,F2=exact_heat_coefficient(F)
        assert o==-sp.Rational(1,3)*F2
        assert s==F2
        assert rat==sp.Rational(2,3)
        out.append((F2,o,s,rat))
    return out


def cyclic_link(K):
    """Clock Z and unitary cyclic-downshift X: X Z X† = exp(2pi i/K) Z."""
    omega=np.exp(2j*np.pi/K)
    Z=np.diag([omega**k for k in range(K)]).astype(complex)
    X=np.zeros((K,K),dtype=complex)
    for k in range(K):X[(k-1)%K,k]=1
    I=np.eye(K)
    return omega,Z,X,(np.linalg.norm(X@Z@X.conj().T-omega*Z),
                     np.linalg.norm(X@X.conj().T-I),
                     np.linalg.norm(Z@Z.conj().T-I))


def finite_ZK_gauge_matter(K):
    """Two spinless fermion modes as two qubits, one K-level gauge qudit.
    Gauge transformations G0,G1 commute exactly with hopping H.
    """
    omega,Z,X,errors=cyclic_link(K)
    assert max(errors)<5e-14
    ann=np.array([[0,1],[0,0]],dtype=complex)
    sigz=np.diag([1,-1]).astype(complex)
    I2=np.eye(2,dtype=complex)
    I4=np.eye(4,dtype=complex)
    c0=np.kron(ann,I2)
    c1=np.kron(sigz,ann)
    assert np.linalg.norm(c0@c1+c1@c0)<1e-14
    assert np.linalg.norm(c0@c1.conj().T+c1.conj().T@c0)<1e-14
    hopping=np.kron(c0.conj().T@c1,Z)
    H=hopping+hopping.conj().T
    n0=c0.conj().T@c0;n1=c1.conj().T@c1
    # exp(-i phi n) c exp(+i phi n) = exp(+i phi)c
    unit_phase=np.exp(-2j*np.pi/K)
    P0=I4+(unit_phase-1)*n0
    P1=I4+(unit_phase-1)*n1
    G0=np.kron(P0,X)
    G1=np.kron(P1,X.conj().T)
    residuals=dict(weyl=errors[0],
        local0=float(np.linalg.norm(G0@H@G0.conj().T-H)),
        local1=float(np.linalg.norm(G1@H@G1.conj().T-H)),
        commutation=float(np.linalg.norm(G0@G1-G1@G0)),
        hermiticity=float(np.linalg.norm(H-H.conj().T)),
        orderK=float(max(np.linalg.norm(np.linalg.matrix_power(G0,K)-np.eye(4*K)),
                         np.linalg.norm(np.linalg.matrix_power(G1,K)-np.eye(4*K)))))
    assert max(residuals.values())<1e-12
    return residuals


def discrete_wilson_gauge_covariance(K=7,L=2,r=1):
    """4D Wilson-Dirac 2^4-site matrix in Z_K link background; tests exact transformation."""
    gamma=[np.array(g.tolist(),dtype=complex) for g in matrices_two_qubit_clifford()]
    I4=np.eye(4,dtype=complex)
    sites=list(itertools.product(range(L),repeat=4))
    d={s:i for i,s in enumerate(sites)}
    rng=np.random.default_rng(2401024+K)
    omega=np.exp(2j*np.pi/K)
    def shift(s,mu,step):
        s=list(s);s[mu]=(s[mu]+step)%L;return tuple(s)
    edge={(s,mu):omega**int(rng.integers(0,K)) for s in sites for mu in range(4)}
    g={s:omega**int(rng.integers(0,K)) for s in sites}
    new={(s,mu):g[s]*edge[(s,mu)]*np.conj(g[shift(s,mu,1)]) for s in sites for mu in range(4)}
    def make_D(u):
        D=np.zeros((4*len(sites),4*len(sites)),dtype=complex)
        mass=sp.Rational(37,100)
        for s,i in d.items():
            sl=slice(4*i,4*i+4)
            D[sl,sl] += (float(mass)+4*r)*I4
            for mu in range(4):
                y=shift(s,mu,1); z=shift(s,mu,-1)
                D[sl,slice(4*d[y],4*d[y]+4)] += -0.5*(r*I4-gamma[mu])*u[(s,mu)]
                D[sl,slice(4*d[z],4*d[z]+4)] += -0.5*(r*I4+gamma[mu])*np.conj(u[(z,mu)])
        return D
    D=make_D(edge); Dp=make_D(new)
    G=np.kron(np.diag([g[s] for s in sites]),I4)
    rel=np.linalg.norm(Dp-G@D@G.conj().T)/np.linalg.norm(D)
    a=np.linalg.slogdet(D);b=np.linalg.slogdet(Dp)
    logdiff=float(abs(a[1]-b[1]));phasediff=float(abs(a[0]-b[0]))
    assert rel<1e-13 and logdiff<1e-11 and phasediff<1e-11
    return dict(K=K, wilson_r=r, covariance=float(rel),logdet=logdiff, phase=phasediff)


def landau_proper_time():
    """Independent spectrum calculation (assumes a continuum 4D Dirac operator).
    From Landau n>=0 and two spin states: K_B/V = 4/(4*pi*t)^2 x*coth(x).
    Extract 4*(x*coth x-1)/x^2 -> 4/3 for F_01=B=1.
    """
    x=sp.symbols('x')
    expression=4*x*sp.coth(x)
    series=sp.series(expression,x,0,8)
    assert sp.expand(series.removeO()).coeff(x,2)==sp.Rational(4,3)
    assert sp.expand(series.removeO()).coeff(x,4)==-sp.Rational(4,45)
    samples=[]
    for xval in (0.4,0.2,0.1,0.05,0.025):
        # n-sum of 2*(e^(-2nx)+e^(-2(n+1)x)): spin degeneracy and Landau ladder
        n=np.arange(0,1000)
        sumspin=float(2*np.sum(np.exp(-2*n*xval)+np.exp(-2*(n+1)*xval)))
        expected=2/math.tanh(xval)
        assert abs(sumspin-expected)<2e-12
        ratio=4*(xval/math.tanh(xval)-1)/(xval*xval)
        samples.append((xval,ratio))
    return series,samples


def electroweak_matter_trace():
    """Conditional ONE generation SM field roster; only Q_L (3,2) comes from paper's algebra.
    Compute U(1)_Y charges and SU(2),SU(3) Dynkin indices (T(fund)=1/2).
    """
    row=[
       # name,color,weak,Y, T_SU2 including color, T_SU3 including weak
       ('Q_L',3,2,sp.Rational(1,6),sp.Rational(3,2),sp.Integer(1)),
       ('u_R',3,1,sp.Rational(2,3),0,sp.Rational(1,2)),
       ('d_R',3,1,sp.Rational(-1,3),0,sp.Rational(1,2)),
       ('L_L',1,2,sp.Rational(-1,2),sp.Rational(1,2),0),
       ('e_R',1,1,sp.Integer(-1),0,0)
    ]
    Ys={name:sp.simplify(colors*weak*Y**2) for name,colors,weak,Y,t2,t3 in row}
    assert Ys['Q_L']==sp.Rational(1,6)
    assert sum(Ys.values())==sp.Rational(10,3)
    sum2=sum(v[4] for v in row);sum3=sum(v[5] for v in row)
    assert sum2==2 and sum3==2
    # Add 3 generations, one complex Higgs doublet with Y=1/2; SU(N) vector loops.
    bY=sp.Rational(2,3)*3*sum(Ys.values())+sp.Rational(1,3)*sp.Rational(1,2)
    b2=-sp.Rational(11,3)*2+sp.Rational(2,3)*3*sum2+sp.Rational(1,3)*sp.Rational(1,2)
    b3=-sp.Rational(11,3)*3+sp.Rational(2,3)*3*sum3
    assert (bY,b2,b3)==(sp.Rational(41,6),-sp.Rational(19,6),-sp.Integer(7))
    return Ys,sum2,sum3,(bY,b2,b3)


def main():
    h=test_gamma_and_heats()
    print('Two-qubit Clifford + heat-kernel tensors:')
    for F2,o,s,r in h:print('  F^2=',F2,'orbital=',o,'spin=',s,'ratio=',r)
    for K in (3,5,7,11): print('Z_K gauge-matter test',K,finite_ZK_gauge_matter(K))
    for K in (5,7):
        for r in (0,1):
            print('4D naive/Wilson Z_K background',discrete_wilson_gauge_covariance(K,r=r))
    # Both kinetic terms obey the same exact gauge symmetry, but their continuum
    # massless matter content is different: 16 tastes (naive) versus 1 (Wilson).
    corner_numbers={}
    for r in (0,1):
        n=0
        for p in itertools.product((0,sp.pi),repeat=4):
            scalar=sp.simplify(r*sum((1-sp.cos(x)) for x in p))
            if scalar==0:n+=1
        corner_numbers[r]=n
    assert corner_numbers=={0:16,1:1}
    print('Same qudit lattice, gauge symmetry; massless tastes naive/Wilson:',corner_numbers)
    ser,landau=landau_proper_time()
    print('Landau spectrum heat trace series:',ser)
    print('Extracted a4 per B^2 from heat trace:',landau)
    y,i2,i3,bs=electroweak_matter_trace()
    print('Single SM generation hypercharge-squared Weyl traces:',y)
    print('Single-generation nonabelian Dynkin index sums SU2, SU3:',i2,i3)
    print('Conditional complete SM one-loop b_Y, b_2, b_3:',bs)
    print('CONDITIONAL DIRAC RESULT: beta(e)=e^3/(12*pi^2), beta(alpha)=2*alpha^2/(3*pi)')
    print('ALL EXTENSION TESTS PASSED')

if __name__=='__main__':main()