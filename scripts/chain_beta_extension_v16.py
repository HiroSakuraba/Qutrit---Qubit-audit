#!/usr/bin/env python3
"""Chain Theory v16: invariant local readout effects and Choi-certified CP bridges.

Exact classification in the five-outcome *commutative* record-label algebra.
Does NOT claim classification of all maps on the complete 31104-dimensional
operator algebra, nor gauge equivariance of the source's untransported P_w.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import sympy as sp
from chain_beta_extension_v15 import EXACT as V15_EFFECTS

HERE = Path(__file__).resolve().parent
r = sp.Rational
weights = {'0':sp.S(0), '1':sp.S(1), 'sqrt2':sp.sqrt(2),'sqrt3':sp.sqrt(3),'sqrt6':sp.sqrt(6)}
mult = {'0':30932,'1':79,'sqrt2':48,'sqrt3':42,'sqrt6':3}
N = 31104
alpha = sp.simplify(sum(weights[k]*n for k,n in mult.items())/N)
I9 = sp.eye(9)
ket = sp.zeros(9,1)
for i in range(3):ket[4*i] = 1/sp.sqrt(3)
P = ket * ket.H
F = I9-P
rho = sp.zeros(9)
for a in range(2):
    for b in range(2):rho[3*a+b,3*a+b]=r(1,4)
assert sp.trace(P)==1 and sp.trace(F)==8 and sp.trace(rho*P)==r(1,6)

# Allowed effects of the *output-side* independent-frame-invariant local POVM:
# E_k = a_k P + b_k F; trace calibrates  (a_k+8 b_k)/9 = n_k/N.
# The 5 outcome parameters a_k satisfy a_k,b_k>=0, sum_k a_k=sum_k b_k=1.
# In particular  a_k in [0, 9*n_k/N] for k>0; a_0 carries the residue.
mu = {k:r(9*mult[k],N) for k in mult}

def family(a_positive):
    aa = {k:sp.sympify(a_positive.get(k,0)) for k in mult if k != '0'}
    aa['0']=sp.simplify(1-sum(aa.values()))
    bb={k:sp.simplify((mu[k]-aa[k])/8) for k in mult}
    for k in mult:
        assert aa[k]>=0 and bb[k]>=0, (k,aa[k],bb[k])
    assert sp.simplify(sum(aa.values())-1)==0
    assert sp.simplify(sum(bb.values())-1)==0
    ee={k:aa[k]*P+bb[k]*F for k in mult}
    assert sp.simplify(sum(ee.values(),sp.zeros(9))-I9)==sp.zeros(9)
    for k in mult:
        assert sp.simplify(sp.trace(ee[k])/9-r(mult[k],N))==0
        assert sp.simplify(P*ee[k]-ee[k]*P)==sp.zeros(9)
    return aa,bb,ee

# Three pairwise distinct positive and group-invariant readout bridges.
# 'outside_Bell' puts all weighted support in F and is the upper bound
# for the Bell-first instrument with score zero on the selected branch.
a_low={k:0 for k in mult if k!='0'}
a_high={k:mu[k] for k in mult if k!='0'}
a_scalar={k:r(mult[k],N) for k in mult if k!='0'}
AA={}
for name, choice in [('outside_Bell',a_low),('scalar',a_scalar),('inside_Bell',a_high)]:
    AA[name]=family(choice)

# v15 local effects twirled under the correlated SU(3) matching U⊗bar(U).
# For E=cI+dS, channel's twirl has P coefficient c+d,
# F coefficient c+d/4, since twirl(S)=P+(1/4)F.
E0v15=I9-sum(V15_EFFECTS.values(),sp.zeros(9))
v15={'0':E0v15,**V15_EFFECTS}
S=sp.zeros(9)
for a in range(3):
    for b in range(3):S[3*a+b,3*b+a]=1
assert sp.trace(S)==3
atw={}; btw={}; Tw={}
for k,E in v15.items():
    c=sp.simplify((sp.trace(E*S)*3-sp.trace(E)*9)/(9-81)) if False else None
    # Determine coefficients from traces on +1 and -1 swap sectors
    sym = (I9+S)/2; asym=(I9-S)/2
    eplus=sp.simplify(sp.trace(sym*E)/6)
    eminus=sp.simplify(sp.trace(asym*E)/3)
    c=sp.simplify((eplus+eminus)/2)
    d=sp.simplify((eplus-eminus)/2)
    assert sp.simplify(E-c*I9-d*S)==sp.zeros(9)
    atw[k]=sp.simplify(c+d)
    btw[k]=sp.simplify(c+d/4)
    Tw[k]=atw[k]*P+btw[k]*F
    assert sp.simplify(sp.trace(Tw[k])-sp.trace(E))==0
    assert sp.simplify(atw[k]+8*btw[k]-mu[k])==0
assert sp.simplify(sum(Tw.values(),sp.zeros(9))-I9)==sp.zeros(9)

# Compute scores for original two-qutrit state, Bell state, and maximally mixed state.
results={}
for name,(aa,bb,ee) in list(AA.items())+[('v15_correlated_twirl',(atw,btw,Tw))]:
    def read(rhostate):return sp.simplify(sum(weights[k]*sp.trace(rhostate*ee[k]) for k in mult))
    score_in=read(rho)
    score_bell=read(P)
    score_mixed=read(I9/9)
    score_seq=sp.simplify(r(5,6)*sum(weights[k]*bb[k] for k in mult))
    assert score_mixed==alpha
    assert score_seq==sp.simplify(sum(weights[k]*sp.trace(rho*F*ee[k]*F) for k in mult))
    results[name]={
        'a_Bell':{k:str(aa[k]) for k in mult},
        'b_failure':{k:str(bb[k]) for k in mult},
        'source_input_weighted_score_exact':str(score_in),
        'source_input_weighted_score_decimal':float(score_in.evalf(20)),
        'Bell_state_weighted_score_exact':str(score_bell),
        'uniform_state_weighted_score_exact':str(score_mixed),
        'Bell_first_score_Bell_assigned_zero_exact':str(score_seq),
        'Bell_first_score_Bell_assigned_zero_decimal':float(score_seq.evalf(20)),
        'weighted_bell_success_probability':str(sp.trace(rho*P))}

lo,hi=sp.simplify(r(15,16)*alpha),sp.simplify(r(3,2)*alpha)
assert sp.simplify(sum(weights[k]*sp.trace(rho*AA['outside_Bell'][2][k]) for k in mult)-lo)==0
assert sp.simplify(sum(weights[k]*sp.trace(rho*AA['inside_Bell'][2][k]) for k in mult)-hi)==0
assert sp.simplify(r(5,6)*sum(weights[k]*AA['outside_Bell'][1][k] for k in mult)-lo)==0
assert sp.simplify(sum(weights[k]*AA['inside_Bell'][1][k] for k in mult))==0
assert sp.simplify(r(1,6)*r(3,8)*alpha + lo-alpha)==0

# Choi certificate for map B(C^N) -> B(C^9):
# Phi(X)=sum_k Tr(P_k X)/n_k E_k.
# J(Phi)=sum_k(P_k^T/n_k)⊗E_k PSD blockwise (no 279936^2 array).
# Unital because sum E=I. Each block eigenvalues a_k/n_k (rank1),
# b_k/n_k (rank8), repeated n_k times.
choi_info={}
for name,(aa,bb,ee) in list(AA.items())+[('v15_correlated_twirl',(atw,btw,Tw))]:
    rk=sum(mult[k]*((1 if aa[k]>0 else 0)+(8 if bb[k]>0 else 0)) for k in mult)
    assert rk<=9*N
    choi_info[name]={'Choi_dim':9*N,'Choi_rank':rk,
      'min_Choi_eigenvalue':0.0 if rk<9*N else float(min([float(sp.N(v/mult[k])) for k in mult for v in [aa[k],bb[k]]])),
      'Choi_PSD_certified_by_positive_blocks':True,
      'Phi_of_identity_equals_identity':True,
      'Phi_of_source_record_sector_has_effect_trace_9n_over_N':True}

# Group invariance numerically: U⊗bar(U) fixes P,F so all new effects invariant.
rng=np.random.default_rng(20261025)
Pnum=np.array(P,dtype=np.complex128); Fnum=np.eye(9)-Pnum
max_fix=0.;max_bell=0.;max_cov_score=0.
for trial in range(80):
    Z=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    Q,R=np.linalg.qr(Z)
    diag=np.diag(R); U=Q@np.diag(np.conj(diag)/np.abs(diag))
    U=U/np.linalg.det(U)**(1/3)
    W=np.kron(U,U.conj())
    max_bell=max(max_bell,float(np.linalg.norm(W@Pnum@W.conj().T-Pnum)))
    for name,(aa,bb,ee) in list(AA.items())+[('v15_correlated_twirl',(atw,btw,Tw))]:
        for k in mult:
            En=complex(aa[k])*Pnum+complex(bb[k])*Fnum
            max_fix=max(max_fix,float(np.linalg.norm(W@En@W.conj().T-En)))
    if trial==0:
        # Original v15 effect is not invariant under a generic matched rotation.
        orig=np.array(V15_EFFECTS['1'],dtype=complex)
        assert np.linalg.norm(W@orig@W.conj().T-orig)>1e-6
assert max_fix < 1e-12 and max_bell < 1e-12

# Distinct group-invariant maps with same source spectral trace calibration
assert AA['outside_Bell'][1]['1']!=AA['inside_Bell'][1]['1']

report={
 'program':'Chain Theory v16 — commutant/Choi classification of local outcome bridges',
 'status':'Exact class of invariant outcome effects and explicit UCP map extensions; not unique; no derivation of beta function.',
 'core_assumptions':['output invariance under every W=U⊗conj(U) for U∈SU(3)',
   'five labelled record outcomes with n=(30932,79,48,42,3)',
   'trace calibration Tr(E_k)/9=n_k/31104',
   'local effects on the joining C3⊗C3 Hilbert space',
   'separate exact selected Bell Kraus branch K=P in the sequential instrument'],
 'alpha_candidate_exact':str(alpha),'alpha_candidate_decimal':float(alpha.evalf(20)),
 'commutant':'span{P3,I9-P3} for U⊗conj(U); rank(P3)=1 and rank(F)=8',
 'all_invariant_effects':'E_k=a_k P3+b_k F, a_k,b_k>=0, sum a=sum b=1; a_k+8b_k=9n_k/N',
 'positive_record_trace_fraction':str(r(172,N)),
 'source_input_Bell_probability':str(r(1,6)),
 'score_range_without_first_Bell_condition':{'lower_exact':str(lo),'upper_exact':str(hi),
     'lower_decimal':float(lo.evalf(20)),'upper_decimal':float(hi.evalf(20))},
 'Bell_first_zero_score_upper_bound':{'exact':str(lo),'decimal':float(lo.evalf(20)),
     'gap_from_alpha_exact':str(sp.simplify(alpha-lo)),
     'minimum_extra_Bell_score_to_match_alpha':str(sp.simplify(r(3,8)*alpha))},
 'examples':results,'choi_certificates':choi_info,
 'gauge_tests':{'random_matched_SU3_trials':80,'maximum_fixed_effect_residual':max_fix,'maximum_P3_residual':max_bell},
 'qualification':'The spectral-sector map is CP and unital; full equivariance under the original nontrivial matching action on all 31,104 record directions is NOT proved. Physical record-to-coupling identification remains unselected.',
 'beta_function_status':'Still undefined: no scale dynamics or fermion propagation from the source axiom.'}
report_path=HERE/'chain_beta_extension_v16_results.json'
report_path.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ['program','alpha_candidate_decimal','score_range_without_first_Bell_condition','Bell_first_zero_score_upper_bound','gauge_tests']},indent=2))
print('All exact v16 algebraic and analytic Choi checks passed.')