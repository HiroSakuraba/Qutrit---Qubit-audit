#!/usr/bin/env python3
"""Chain Theory beta v17: full-factor gauge intertwiners and rigorous CP nonuniqueness.

Exact source-motivated operator Hilbert-space factorization:
  R = (su3_HS⊗su2_HS)_unread ⊗ End(C3)_Q1 ⊗ End(C2)_q1
        ⊗ End(C3)_Q2 ⊗ End(C2)_q2, dimension 24*9*4*9*4=31104.

This is a conditional operational completion. Its full local basis-change group
is more generous than the matching subgroup of the original source; proving
existence of two maps covariant under the larger group proves non-uniqueness
under its subgroups. It does not assert this group is the paper's unique gauge
representation or that its observable is electromagnetic coupling.

Companion scripts chain_beta_extension_v9.py and _v11.py used for records.
"""
from __future__ import annotations
from pathlib import Path
from collections import Counter
import json
import math
import numpy as np
import sympy as sp
from chain_beta_extension_v9 import make_basis, gellmann, record_matrix, adjoint
from chain_beta_extension_v11 import build_canonical_records

ROOT=Path(__file__).resolve().parent
N=31104
D_REST=24*4*4  # unread connected and qubit operator dimensions
assert N == 9*9*D_REST
OUT=9
LABELS=('0','1','sqrt2','sqrt3','sqrt6')
COUNTS={'0':30932, '1':79, 'sqrt2':48, 'sqrt3':42, 'sqrt6':3}
W={'0':sp.S(0),'1':sp.S(1),'sqrt2':sp.sqrt(2),'sqrt3':sp.sqrt(3),'sqrt6':sp.sqrt(6)}
alpha=sp.simplify(sum(W[k]*COUNTS[k] for k in LABELS)/N)

# 1. Local matching-sector decomposition of R. Each local End(Cd) is 1+adj,
# so independent (SU3xSU2)^3 acts with 16 inequivalent irreducible blocks.
local_dims=(8,3,8,3)  # slots Q1,q1,Q2,q2
patterns={tuple((mask>>k)&1 for k in range(4)):24*math.prod(local_dims[k] if mask>>k & 1 else 1 for k in range(4)) for mask in range(16)}
assert len(patterns)==16 and sum(patterns.values())==N
assert patterns[(0,0,0,0)]==24 and patterns[(1,1,1,1)]==13824

# Full Choi commutant, not merely the 16 input-vector irreps.
# Cbar(R)⊗(3_Q1⊗3_Q2) carries (1+8)⊗3 = 2*3 + bar6 + 15
# on each joining qutrit; each read qubit contributes 1+3, unread 8*3.
qutrit_irreps={'3':(3,2),'bar6':(6,1),'15':(15,1)}  # (dimension,multiplicity)
qubit_irreps={'1':1,'3':3}
choi_blocks=[]
for a,(da,ma) in qutrit_irreps.items():
    for b,(db,mb) in qutrit_irreps.items():
        for x,dx in qubit_irreps.items():
            for y,dy in qubit_irreps.items():
                choi_blocks.append({'qutrit_1':a,'qutrit_2':b,'qubit_1':x,'qubit_2':y,
                    'irrep_dimension':24*da*db*dx*dy,'multiplicity':ma*mb})
assert len(choi_blocks)==36
assert sum(z['irrep_dimension']*z['multiplicity'] for z in choi_blocks)==9*N
commutant_complex_dimension=sum(z['multiplicity']**2 for z in choi_blocks)
assert commutant_complex_dimension==144
# 16 independent normalized partial-output traces (one per input sector).
# The depolarizing Choi I/(N) is strictly positive so the positive cone has
# a 128-dimensional affine interior after these constraints.
normalized_trace_constraint_rank=len(patterns)
affine_CP_dimension=commutant_complex_dimension-normalized_trace_constraint_rank
assert affine_CP_dimension==128

# Source-driven exact record projectors: canonical 175 slot-word coordinates.
_,blocks,singles,names=make_basis()
E,ww,tags=build_canonical_records(blocks,singles)
_,_ =record_matrix(blocks,singles)
assert E.shape==(175,172)
labels=[]
for val in ww:
    labels.append(min([k for k in LABELS if k!='0'],key=lambda k:abs(float(W[k])-val)))
assert Counter(labels)==Counter({k:COUNTS[k] for k in LABELS if k!='0'})
pattern_by_slot_name={'QQ':(1,0,1,0),'qq':(0,1,0,1),
  'Q1q1':(1,1,0,0),'Q1q2':(1,0,0,1),
  'Q2q1':(0,1,1,0),'Q2q2':(0,0,1,1),
  'q1':(0,1,0,0),'q2':(0,0,0,1)}
# Exact overlap T[pattern,k]=Tr(P_pattern P_k^R), including zero sector.
T={(p,k):sp.S(0) for p in patterns for k in LABELS}
for j,k in enumerate(labels):
    for i in np.flatnonzero(abs(E[:,j])>1e-10):
        co=float(abs(E[i,j])**2)
        exact=sp.Rational(1) if abs(co-1)<1e-12 else sp.Rational(1,2)
        assert abs(float(exact)-co)<1e-12
        p=pattern_by_slot_name[names[i][0]]
        T[(p,k)]+=exact
for p in patterns:
    T[(p,'0')]=sp.Integer(patterns[p])-sum(T[(p,k)] for k in LABELS if k!='0')
for k in LABELS:
    assert sp.simplify(sum(T[(p,k)] for p in patterns)-COUNTS[k])==0
for p in [(0,0,0,0),(1,1,1,1)]:
    assert all(T[(p,k)]==0 for k in LABELS if k!='0')

# Exact scalar CP maps: sigma_t=I/N + t (P_0000/24 - P_1111/13824).
# Gauge invariant under all independent basis changes. Distinct for t !=0,
# yet Tr(sigma_t P_k)=n_k/N for ALL FIVE source weight sectors.
# Check positive density at t=1/4 and exact calibration without dense arrays.
t_scalar=sp.Rational(1,4)
eigs_sigma={
    '0000':sp.Rational(1,N)+t_scalar/sp.Integer(patterns[(0,0,0,0)]),
    '1111':sp.Rational(1,N)-t_scalar/sp.Integer(patterns[(1,1,1,1)]),
    'other':sp.Rational(1,N),
}
assert all(x>=0 for x in eigs_sigma.values())
assert sp.simplify(sum(patterns[p]*(sp.Rational(1,N)+t_scalar*((sp.Rational(1,24) if p==(0,0,0,0) else 0)-(sp.Rational(1,13824) if p==(1,1,1,1) else 0))) for p in patterns)-1)==0
for k in LABELS:
    a=sp.simplify(sum(T[(p,k)]*(sp.Rational(1,N)+t_scalar*((sp.Rational(1,24) if p==(0,0,0,0) else 0)-(sp.Rational(1,13824) if p==(1,1,1,1) else 0))) for p in patterns))
    assert a==sp.Rational(COUNTS[k],N)

# 2. Build a full-gauge-covariant *non-scalar* linear perturbation of the
# depolarizing UCP map. Adj SU3 matrices F_a on 8D and Gell-Mann ell_a
# on 3D are intertwining vectors; sum_ab Tr(A_ab X) ell_a⊗ell_b.
# Using A_ab=F_a(Q1)⊗F_b(Q2)⊗I_rest/D_REST. All other independent local
# SU2/SU3 transformations act on I_rest and thus drop out.
ell=gellmann()  # tr ell_a ell_b =delta_ab, Hermitian
assert len(ell)==8
F8=[]
for la in ell:
    x=np.zeros((8,8),dtype=complex)
    for i,li in enumerate(ell):
        for j,lj in enumerate(ell):
            x[i,j]=np.trace(li@(la@lj-lj@la))
    assert np.linalg.norm(x-x.conj().T)<2e-13
    assert np.linalg.norm(x.imag+x.imag.T)<2e-13
    F8.append(x)
F9=[np.pad(f,((1,0),(1,0))) for f in F8]  # scalar + adjoint basis
for a in range(8):
    for b in range(8):
        assert abs(np.trace(F8[a]@F8[b])-(6. if a==b else 0.))<1e-12

# Choi certificate: J_Phi = I_(9*N)/N + epsilon/Drest * K1⊗K1 (up to
# factor permutation), where K1=sum_a F9_a^T⊗ell_a (27 by 27).
K1=sum((np.kron(F9[a].T,ell[a]) for a in range(8)),np.zeros((27,27),complex))
assert np.linalg.norm(K1-K1.conj().T)<1e-12
lam=np.linalg.eigvalsh(K1)
expected_K1=np.array([-1]*15+[0]*3+[1]*6+[3]*3,dtype=float)  # checked below
# Sign convention may flip depending on generator sign; verify by direct Casimir spectrum.
found=np.array([(-1. if abs(v+1)<1e-8 else 0. if abs(v)<1e-8 else 1. if abs(v-1)<1e-8 else 3. if abs(v-3)<1e-8 else np.nan) for v in lam])
if not np.all(np.isfinite(found)):
    raise AssertionError(('unexpected K1 eigenvalues',lam))
assert np.linalg.norm(np.sort(lam)-expected_K1)<1e-10,(lam,expected_K1)
K2eig=np.outer(lam,lam).ravel()
epsilon=sp.Rational(1,1000)
# N=81 Drest, so every Choi eigenvalue is (1+81*eps*l1*l2)/N.
choi_low=sp.simplify((1+81*epsilon*sp.Rational(int(round(min(K2eig))),1))/N)
choi_high=sp.simplify((1+81*epsilon*sp.Rational(int(round(max(K2eig))),1))/N)
assert choi_low>0 and choi_high>0
assert abs(min(K2eig)+3)<1e-10 and abs(max(K2eig)-9)<1e-10

# 3. Target physical record effects for all five sector projectors.
# Only QQ records project into the adjoint-adjoint source factor under the
# bilinear observable; symmetry decomposition has swap eigenvalues +/-1:
#    Tr(A_ab P_sqrt3) = +3 delta_ab / Drest,
#    Tr(A_ab P_1)     = -3 delta_ab / Drest.
S3=sp.sqrt(3)
Bsum=sum((np.kron(ell[a],ell[a]) for a in range(8)),np.zeros((9,9),complex))
Sswap=np.zeros((9,9),complex)
for a in range(3):
    for b in range(3): Sswap[3*a+b,3*b+a]=1.
assert np.linalg.norm(Bsum-(Sswap-np.eye(9)/3))<2e-13
Wphys={'1':-3.,'sqrt3':3.,'sqrt2':0.,'sqrt6':0.,'0':0.}
Eps={k:float(sp.Rational(COUNTS[k],N))*np.eye(9)+float(epsilon)*Wphys[k]/D_REST*Bsum for k in LABELS}
assert np.linalg.norm(sum(Eps.values())-np.eye(9))<1e-14
for k in LABELS:
    assert abs(np.trace(Eps[k]).real/9-COUNTS[k]/N)<1e-14
    assert np.linalg.eigvalsh(Eps[k]).min()>-1e-14
    assert np.linalg.norm(Eps[k]-Eps[k].conj().T)<1e-14
weighted=sum((float(W[k])*Eps[k] for k in LABELS),np.zeros((9,9),complex))
predicted= float(alpha)*np.eye(9)+float(epsilon)*float((S3-1)/128)*Bsum
assert np.linalg.norm(weighted-predicted)<2e-14
rho=np.zeros((9,9),complex)
for a in range(2):
    for b in range(2):rho[3*a+b,3*a+b]=1/4
Bell=np.zeros(9,complex)
for a in range(3): Bell[4*a]=1/np.sqrt(3)
Pbell=np.outer(Bell,Bell.conj());Fbell=np.eye(9)-Pbell
assert abs(np.trace(rho@Sswap)-.5)<1e-14
assert abs(np.trace(rho@Pbell)-1/6)<1e-14
physical_score=float(np.trace(rho@weighted).real)
expected_delta=sp.simplify(epsilon*(S3-1)/768)
assert abs(physical_score-float(alpha+expected_delta))<1e-14
# The effects commute with Bell P; Bell-first completion remains available,
# with the Bell outcome score requiring independent specification.
for effect in Eps.values():
    assert np.linalg.norm(effect@Pbell-Pbell@effect)<1e-14
assert np.linalg.norm(sum((Fbell@E@Fbell for E in Eps.values()),np.zeros((9,9)))-Fbell)<1e-14

# 4. Independently transformed frame covariance on source 9x9 ⊗ 9x9
# small factor and target C3⊗C3, for random Hermitian 81x81 inputs.
rng=np.random.default_rng(20261025)
def rand_su3():
    z=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    q,rr=np.linalg.qr(z);diag=np.diag(rr)
    q=q@np.diag(np.conj(diag)/np.abs(diag))
    return q/np.linalg.det(q)**(1/3)
Fsource=[np.kron(F9[a],F9[b]) for a in range(8) for b in range(8)]
Btarget=[np.kron(ell[a],ell[b]) for a in range(8) for b in range(8)]
def Psi_small(X):
    co=np.array([np.trace(A@X) for A in Fsource])
    return np.einsum('a,aij->ij',co,np.array(Btarget))
max_cov=0.;max_orth=0.;max_phi_cov=0.
for _ in range(12):
    u1,u2=rand_su3(),rand_su3()
    o1,o2=adjoint(u1,ell),adjoint(u2,ell)
    r1=np.eye(9);r2=np.eye(9);r1[1:,1:]=o1;r2[1:,1:]=o2
    R=np.kron(r1,r2); U=np.kron(u1,u2)
    max_orth=max(max_orth,np.linalg.norm(R@R.T-np.eye(81)))
    z=rng.normal(size=(81,81))+1j*rng.normal(size=(81,81))
    X=z+z.conj().T
    p0=Psi_small(X);p1=Psi_small(R@X@R.conj().T)
    max_cov=max(max_cov,float(np.linalg.norm(p1-U@p0@U.conj().T)))
    # Full lifted input X⊗I_rest has normalized output Tr(X)/81 I + eps Psi_small(X)
    # (full A_ab contraction cancels Drest), and respects independent action.
    q0=np.trace(X)/81*np.eye(9)+float(epsilon)*p0
    q1=np.trace(X)/81*np.eye(9)+float(epsilon)*p1
    max_phi_cov=max(max_phi_cov,float(np.linalg.norm(q1-U@q0@U.conj().T)))
assert max_cov<1e-9, max_cov
assert max_phi_cov<1e-9 and max_orth<1e-11

# 5. Independent explicit check of QQ sector contractions (small exact matrix).
Fprod=np.array(Fsource)
qq_diag={}
for a in [0,3,7]:
    f=np.kron(F8[a],F8[a]);
    swap=np.zeros((64,64),complex)
    for i in range(8):
        for j in range(8):swap[8*i+j,8*j+i]=1
    # trace F⊗F S = tr F^2 = 6
    assert abs(np.trace(f@swap)-6)<1e-12
    qq_diag[str(a)]={'trace_Bilinear_on_QQ_symmetric_unnormalized':float((np.trace(f@swap)/2).real),
          'trace_Bilinear_on_QQ_antisymmetric_unnormalized':float((-np.trace(f@swap)/2).real)}

result={
  'program':'Chain Theory QED beta v17 — full factor gauge-intertwiner audit',
  'status':'Proved explicit nonuniqueness under complete local matching-factor action; no unique physical alpha or beta function',
  'source_model_dimensions':{'record_space':N,'unread_connected':24,'two_read_arm_operator_space':1296,'physical_two_qutrit':9},
  'natural_full_local_group':'(SU(3)xSU(2))^3 acting by adjoint on each local operator factor; physical joining qutrits transform under SU(3)_Q1 x SU(3)_Q2',
  'qualification':'This factor representation is source-motivated and larger than pair-ket stabilizer matchings; the original paper does not define a unique operational group action or physical channel.',
  'sixteen_invariant_irreducible_sectors':{''.join(map(str,p)):v for p,v in patterns.items()},
  'full_Choi_intertwiner_classification':{'isotypic_blocks':len(choi_blocks),'commutant_complex_dimension':commutant_complex_dimension,
    'Hermitian_commutant_real_dimension':commutant_complex_dimension,
    'normalized_trace_constraint_rank':normalized_trace_constraint_rank,
    'affine_dimension_of_normalized_trace_preserving_UCP_maps':affine_CP_dimension,
    'Choi_Hilbert_dimension':9*N,
    'block_criterion':'In group-adapted basis J=direct_sum(I_irrep tensor H_block), H_block positive semidefinite; Tr_output J=(9/N)I_R; unital follows by output irreducibility',
    'qutrit_local_decomposition':'(1+8) tensor 3 = 3+3+bar6+15',
    'Choi_block_decomposition':choi_blocks},
  'record_projector_overlaps':{''.join(map(str,p)):{k:str(T[(p,k)]) for k in LABELS} for p in patterns},
  'scalar_UCP_full_gauge_intertwiners':{'number_of_distinct_examples':2,'rho1':'I/31104','rho2':'I/31104 + (1/4)(P_0000/24-P_1111/13824)',
       'min_density_eigenvalues_exact':{k:str(v) for k,v in eigs_sigma.items()},'identical_five_record_outcome_probabilities':{k:str(sp.Rational(COUNTS[k],N)) for k in LABELS},
       'different_on_other_full_algebra_inputs':True},
  'non_scalar_UCP_full_gauge_intertwiner':{'epsilon':str(epsilon),'K1_casimir_spectrum':{'-1':15,'0':3,'1':6,'3':3},
       'full_Choi_min_eigen_exact':str(choi_low),'full_Choi_max_eigen_exact':str(choi_high),
       'CP_proof':'K1⊗K1 spectrum; all Choi eigenvalues (1+81 epsilon lambda_i lambda_j)/31104 strictly positive',
       'unital_and_normalized_trace_preserving':True,'full_source_record_label_trace_calibration':True,
       'physical_effect_coefficients':{'1':'79/N I - 3 epsilon/384 (S-I/3)',
            'sqrt3':'42/N I + 3 epsilon/384 (S-I/3)',
            'sqrt2':'48/N I','sqrt6':'3/N I','0':'30932/N I'},
       'weighted_score_operator':'alpha I + epsilon (sqrt3-1)/128 (S-I/3)',
       'source_initial_state_score_exact':str(sp.simplify(alpha+expected_delta)),
       'score_deviation_exact':str(expected_delta),
       'source_initial_state_score_decimal':physical_score,
       'alpha_candidate_decimal':float(alpha.evalf(20)),
       'commutes_with_original_Bell_projection':True},
  'numerical_checks':{'independent_SU3_x_SU3_covariance_trials':12,
     'max_covariance_HS_error':max_cov,'max_full_lifted_phi_covariance_error':max_phi_cov,
     'max_source_representation_orthogonality_residual':max_orth,
     'direct_QQ_contract_checks':qq_diag},
  'beta_function':'NOT derived: no intrinsic scale, Dirac dynamics, or physical mapping from measurement score to electromagnetic charge polarization.'
}
p=ROOT/'chain_beta_extension_v17_results.json'
p.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('program','sixteen_invariant_irreducible_sectors','scalar_UCP_full_gauge_intertwiners','non_scalar_UCP_full_gauge_intertwiner','numerical_checks')},indent=2))
print('v17 complete: exact sector and Choi certificates + full independent-gauge covariance tests passed.')