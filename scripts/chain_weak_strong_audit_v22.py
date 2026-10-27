#!/usr/bin/env python3
"""Independent weak/strong coupling diagnostic audit of Chain Theory (2026-10-27).

Outputs JSON with provenance and strict, conditional interpretations.
All numbers involving electroweak parameters at MZ are inputs from PDG 2025,
not predictions of Chain Theory. One-loop RG illustration imports the SM.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import sympy as sp

OUT = Path('/mnt/data/chain_weak_strong_audit_v22_results.json')
R=sp.Rational
alpha_candidate_exact=(79+48*sp.sqrt(2)+42*sp.sqrt(3)+3*sp.sqrt(6))/31104
alpha_candidate=float(alpha_candidate_exact)
alpha_EM_MZ_inv=127.930 # PDG 2025 review MS electromagnetic 5-active flavor illustrative
sin2_MZ=0.23122 # PDG 2025 electroweak review MS
alpha_s_MZ=0.1180 # PDG 2025 QCD review
sigma_alpha_s=0.0009
mZ=91.1880 # GeV, benchmark used here
alpha_EM=1/alpha_EM_MZ_inv
alpha2=alpha_EM/sin2_MZ
alphaY=alpha_EM/(1-sin2_MZ)
alpha1=(5/3)*alphaY
alpha3=alpha_s_MZ
bs={'alpha1':R(41,10),'alpha2':-R(19,6),'alpha3':-R(7)}
alpha0={'alpha1':alpha1,'alpha2':alpha2,'alpha3':alpha3}
assert 1/alpha_EM-(1/alpha2+1/alphaY) < 1e-10
assert abs(alpha1/alphaY-5/3)<1e-13

# SU(5) embedding and trace normalization
T3=sp.diag(R(1,2),-R(1,2),0,0,0)
Y=sp.diag(R(1,2),R(1,2),-R(1,3),-R(1,3),-R(1,3))
assert sp.trace(T3*T3)==R(1,2)
assert sp.trace(Y*Y)==R(5,6)
GUT_hypercharge_scale=sp.sqrt(R(3,5))
assert sp.simplify(sp.trace((GUT_hypercharge_scale*Y)**2)-R(1,2))==0
pred_sin2_unified=R(3,8)

# Conditional: a single chiral (3,2)_1/6 quark multiplet as sole fermion
fields=[
  ('Q_L',3,2,R(1,6),1),
  ('u_R^c',3,1,-R(2,3),-1),
  ('d_R^c',3,1,R(1,3),-1),
  ('L_L',1,2,-R(1,2),0),
  ('e_R^c',1,1,R(1),0),
]
def anomaly(flds):
    cubic=sum(sp.Integer(n3*n2)*hy**3 for _,n3,n2,hy,_ in flds)
    grav=sum(sp.Integer(n3*n2)*hy for _,n3,n2,hy,_ in flds)
    # SU2^2 U1 (only doublets); Dynkin T2 = 1/2
    su2sq_y=sum(sp.Integer(n3)*hy*R(1,2) for _,n3,n2,hy,_ in flds if n2==2)
    # SU3^2 U1: for 3 and bar3 T3=1/2
    su3sq_y=sum(sp.Integer(n2)*hy*R(1,2) for _,n3,n2,hy,_ in flds if n3==3)
    su3cub=sum(sp.Integer(n2*sign) for _,n3,n2,hy,sign in flds if n3==3)
    doublets=sum(n3 for _,n3,n2,_,_ in flds if n2==2)
    return {'SU3_cubic':str(su3cub),'SU3_sq_U1':str(su3sq_y),'SU2_sq_U1':str(su2sq_y),'U1_cubic':str(cubic),'gravity_U1':str(grav),'SU2_doublet_count':int(doublets),'SU2_Witten_anomaly_if_odd':bool(doublets%2)}
a_link=anomaly(fields[:1]); a_full=anomaly(fields)
for k in ['SU3_cubic','SU3_sq_U1','SU2_sq_U1','U1_cubic','gravity_U1']:
    assert a_full[k]=='0'
assert a_full['SU2_doublet_count']==4 and not a_full['SU2_Witten_anomaly_if_odd']
assert a_link['SU2_doublet_count']==3 and a_link['SU2_Witten_anomaly_if_odd']

# Illustrative one-loop SM evolution, GUT U1 normalization.
def inv_at(name,log_mu_mz):
    return 1/alpha0[name] - float(bs[name])/(2*math.pi)*log_mu_mz
cross={}
for a,b in [('alpha1','alpha2'),('alpha1','alpha3'),('alpha2','alpha3')]:
    t=2*math.pi*((1/alpha0[a])-(1/alpha0[b]))/(float(bs[a])-float(bs[b]))
    mu=mZ*math.exp(t)
    cross[f'{a}_{b}']={'log_mu_over_MZ':t,'mu_GeV':mu,'common_alpha_inv':inv_at(a,t),'third_alpha_inv':{c:inv_at(c,t) for c in alpha0 if c not in (a,b)}}

# Exact rank/multiplicities of adjoints, and elementary dimension-only guesses.
ratios={
  'alpha2_over_alphaEMatMZ':alpha2/alpha_EM,
  'alpha3_over_alphaEMatMZ':alpha3/alpha_EM,
  'alpha3_over_alpha2':alpha3/alpha2,
  'alpha2_over_alpha_candidate':alpha2/alpha_candidate,
  'alpha3_over_alpha_candidate':alpha3/alpha_candidate,
  'alpha_candidate_times_3':alpha_candidate*3,
  'alpha_candidate_times_8':alpha_candidate*8,
  'dimension_3_8_guesses_relative_error_alpha2':(alpha_candidate*3/alpha2-1),
  'dimension_3_8_guesses_relative_error_alpha3':(alpha_candidate*8/alpha3-1),
}

out={
  'scope':'Weak/strong theory audit; Chain Theory original paper October 23 2026; external PDG 2025/2026 benchmarks.',
  'paper_source':'ChainTheoryAlpha.pdf, Letter eqs 3-7, Supplemental S3, S8, S9, S10.',
  'alpha_candidate_exact':str(alpha_candidate_exact),
  'alpha_candidate':alpha_candidate,
  'external_inputs_not_chain_theory_predictions':{
    'alpha_EM_MSbar_5flavor_MZ_inv':alpha_EM_MZ_inv,
    'sin2_thetaW_MSbar_MZ':sin2_MZ,
    'alpha_s_MSbar_MZ':alpha_s_MZ,
    'alpha_s_uncertainty':sigma_alpha_s,
    'MZ_GeV':mZ,
    'notes':'EM input is five-flavor MSbar as quoted PDG EW review; conversion to unbroken-SM electroweak couplings is illustrative at leading matching order. RG exercise one-loop with full SM matter and no thresholds above MZ.'},
  'couplings_at_MZ':{
    'alpha_em':alpha_EM,'alpha1_GUT_normalized':alpha1,'alpha2_weak':alpha2,'alpha3_strong':alpha3,
    'g1':math.sqrt(4*math.pi*alpha1),'g2':math.sqrt(4*math.pi*alpha2),'g3':math.sqrt(4*math.pi*alpha3),
    'alpha_inverse':{k:1/v for k,v in alpha0.items()}},
  'group_theory':{
    'dim_SU3_adjoint':8,'dim_SU2_adjoint':3,
    'Tr_Y_squared_in_fundamental_SU5':str(sp.trace(Y**2)),
    'Tr_T3_squared_in_fundamental_SU5':str(sp.trace(T3**2)),
    'canonical_SU5_hypercharge_scaling':str(GUT_hypercharge_scale),
    'unification_sin2_thetaW_at_unification_scale_conditional':str(pred_sin2_unified),
    'is_consequence_of_subgroup_alone':False,
  },
  'conditional_chiral_anomalies':{'isolated_link_as_single_left_chiral_quark_doublet':a_link,'full_one_generation_SM_with_left_handed_conjugates':a_full,'caveat':'The paper identifies a representation, not a dynamical chiral fermion spectrum. These are conditional consistency tests, not demonstrated anomalies of Chain Theory.'},
  'SM_one_loop_b_coefficients_gut_normalized':{k:str(v) for k,v in bs.items()},
  'SM_one_loop_pairwise_crossings':cross,
  'simple_multiplicity_guesses':ratios,
  'limits':[
    'SM beta function coefficients and matter representation content are externally supplied, not derived from Chain Theory.',
    'Gauge-algebra isomorphism fixes structure constants and representations but not kinetic coefficients 1/g_i^2.',
    'No inherent scale, propagator, chiral spectrum, Higgs mechanism or SU3 confinement in source construction.',
    'One-loop crossing scale numbers are illustrative, not precise GUT tests; thresholds and higher loops omitted.',
    'Do not compare alpha_cand at Q=0 directly to SU2 and SU3 coupling values at MZ as if same scale.'
  ]
}
OUT.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['couplings_at_MZ','conditional_chiral_anomalies','SM_one_loop_pairwise_crossings','simple_multiplicity_guesses']},indent=2))
print('PASS: exact group trace identities, benchmark couplings, anomaly cancellations, one-loop crossing diagnostics')