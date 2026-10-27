#!/usr/bin/env python3
"""Chain Theory beta program v21: source-exchange signs versus physical chronology.

Primary-source facts: two commuting record exchanges X,Y with signed sectors
(79,48,42,3), no time order for joined reads, single-copy acts, and Bell
selected contractions. None of the diagnostic quantum channels below is
claimed to follow as physical dynamics from the source axiom.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'chain_beta_extension_v21_results.json'
checks=[]
def check(name, condition):
    if not bool(condition):
        raise AssertionError(name)
    checks.append(name)

# Exact source-derived four signed exchange-character sectors.
sectors = [
    ('1',+1,+1,79,sp.Integer(1)),
    ('sqrt2',+1,-1,48,sp.sqrt(2)),
    ('sqrt3',-1,+1,42,sp.sqrt(3)),
    ('sqrt6',-1,-1,3,sp.sqrt(6)),
]
N=sum(s[3] for s in sectors)
check('172 source record eigenvectors',N==172)
x=np.concatenate([np.full(n,xv,dtype=np.int8) for _,xv,yv,n,w in sectors])
y=np.concatenate([np.full(n,yv,dtype=np.int8) for _,xv,yv,n,w in sectors])
weights=np.concatenate([np.full(n,float(w)) for _,xv,yv,n,w in sectors])
check('exchange involutions',np.all(x*x==1) and np.all(y*y==1))
check('exact commuting exchange action',np.all(x*y==y*x))
check('character traces',sum(n*xv for _,xv,yv,n,_ in sectors)==82 and sum(n*yv for _,xv,yv,n,_ in sectors)==70 and sum(n*xv*yv for _,xv,yv,n,_ in sectors)==-8)
# Each projector is (I +/- X)(I +/- Y)/4, exact spectral multiplicities.
for label,xx,yy,count,_ in sectors:
    mask=(1+xx*x)*(1+yy*y)//4
    check(f'character projector {label}',np.all((mask==0)|(mask==1)) and int(mask.sum())==count)
alpha=sp.simplify(sum(n*w for _,xx,yy,n,w in sectors)/31104)
check('exact 172-record weighted trace',sp.simplify(alpha-(79+48*sp.sqrt(2)+42*sp.sqrt(3)+3*sp.sqrt(6))/31104)==0)

# A diagnostic random-unitary semigroup on the *auxiliary* record Hilbert space.
# No physical time, no read repetition, and no stochastic exchange law is asserted by the paper.
rx=x[:,None]*x[None,:]
ry=y[:,None]*y[None,:]
mode_masks={
    'fixed':(rx==1)&(ry==1),
    'X_only':(rx==-1)&(ry==1),
    'Y_only':(rx==1)&(ry==-1),
    'both':(rx==-1)&(ry==-1),
}
mode_counts={k:int(a.sum()) for k,a in mode_masks.items()}
check('full 172 squared superoperator dimension',sum(mode_counts.values())==N*N)
check('exact conjugation-mode multiplicities',mode_counts=={'fixed':10318,'X_only':6924,'Y_only':7836,'both':4506})
check('X and Y mode dimension identities',mode_counts['X_only']+mode_counts['both']==2*127*45 and mode_counts['Y_only']+mode_counts['both']==2*121*51)

def semigroup_on_operator(mat,gamma_x,gamma_y,t):
    # exp[t gamma_x (Ad_X - I) + t gamma_y (Ad_Y - I)] acting entrywise
    eigen=gamma_x*(rx-1)+gamma_y*(ry-1)
    return np.exp(t*eigen)*mat

rng=np.random.default_rng(20261025)
z=rng.normal(size=(N,5))+1j*rng.normal(size=(N,5))
rho=z@z.conj().T
rho/=np.trace(rho)
check('record-state PSD',np.linalg.eigvalsh(rho).min()>-1e-12)
max_trace_err=0.0
max_semigroup_err=0.0
for gx,gy in [(0.,0.),(.3,.7),(.9,.1),(.5,.5)]:
    a=semigroup_on_operator(rho,gx,gy,0.6)
    check(f'CP semigroup density PSD gamma={gx},{gy}',np.linalg.eigvalsh(a).min()>-1e-11)
    trace_err=abs(np.trace(a)-np.trace(rho))
    max_trace_err=max(max_trace_err,float(trace_err))
    check(f'TP semigroup gamma={gx},{gy}',trace_err<2e-12)
    diff=semigroup_on_operator(semigroup_on_operator(rho,gx,gy,.4),gx,gy,.2)-a
    max_semigroup_err=max(max_semigroup_err,float(np.linalg.norm(diff)))
    check(f'composition gamma={gx},{gy}',np.linalg.norm(diff)<2e-12)
# An explicit Kraus random-unitary realization: independent Poisson parity flips.
def random_unitary_mix(mat,gx,gy,t):
    px=(1-np.exp(-2*gx*t))/2
    py=(1-np.exp(-2*gy*t))/2
    ans=np.zeros_like(mat)
    for ax in (0,1):
        for ay in (0,1):
            q=(px if ax else 1-px)*(py if ay else 1-py)
            u=(x if ax else np.ones(N,dtype=np.int8))*(y if ay else np.ones(N,dtype=np.int8))
            ans+=q*u[:,None]*mat*u[None,:]
    return ans
for gx,gy in [(.2,.8),(.8,.2),(.3,.3)]:
    check(f'exact random unitary parity process gx={gx}, gy={gy}',np.linalg.norm(random_unitary_mix(rho,gx,gy,.75)-semigroup_on_operator(rho,gx,gy,.75))<2e-12)
check('different rates cause different physical auxiliary dynamics',np.linalg.norm(semigroup_on_operator(rho,.2,.8,1)-semigroup_on_operator(rho,.8,.2,1))>1e-3)
# The given sector-weight observable M is diagonal, hence a fixed point of all rates.
check('all source record weights are unchanged by any rates',np.allclose(semigroup_on_operator(np.diag(weights),.7,.3,2),np.diag(weights)))

# Pauli read-order signs: AB=-BA, but each product as a single unknown-phase
# Kraus operator defines the same conventional CP channel.
I=sp.eye(2)
pauli=[sp.Matrix([[0,1],[1,0]]),sp.Matrix([[0,-sp.I],[sp.I,0]]),sp.diag(1,-1)]
for i in range(3):
    for j in range(i+1,3):
        a,b=pauli[i],pauli[j]
        check(f'Pauli anticommutes {i}{j}',a*b==-b*a)
        check(f'Pauli two-step product unitary {i}{j}',sp.simplify((a*b).H*(a*b))==I)
        r=sp.Matrix([[sp.Rational(2,3),sp.Rational(1,5)+sp.I/12],[sp.Rational(1,5)-sp.I/12,sp.Rational(1,3)]])
        check(f'Pauli-order channels identical {i}{j}',sp.simplify((a*b)*r*(a*b).H-(b*a)*r*(b*a).H)==sp.zeros(2))

# Explicit exact Bell contractions of three disjoint initial maximally
# entangled pairs. Two joins on (1,2) and (3,4) in zero-based 6-qudit slots,
# acting on distinct qudit-copy slots, commute. No selection of temporal order.
def bell(d):
    return np.eye(d,dtype=complex)/np.sqrt(d)

def pair_projection(tensor,i,j,d):
    moved=np.moveaxis(tensor,[i,j],[0,1])
    selected=np.einsum('aa...->...',moved)/np.sqrt(d)
    output=np.einsum('ab,...->ab...',bell(d),selected)
    return np.moveaxis(output,[0,1],[i,j])

bell_stats=[]
for d in (2,3,4):
    b=bell(d)
    incoming=np.einsum('ab,cd,ef->abcdef',b,b,b)
    first=pair_projection(incoming,1,2,d)
    second=pair_projection(first,3,4,d)
    first_reverse=pair_projection(incoming,3,4,d)
    second_reverse=pair_projection(first_reverse,1,2,d)
    target=np.einsum('af,bc,de->abcdef',b,b,b)
    eps_order=float(np.linalg.norm(second-second_reverse))
    eps_result=float(np.linalg.norm(second-d**(-2)*target))
    prob=float(np.linalg.norm(second.ravel())**2)
    check(f'chain d={d} order commutes',eps_order<1e-13)
    check(f'chain d={d} exact two-join amplitude',eps_result<1e-13)
    check(f'chain d={d} Born selected branch prob d^-4',abs(prob-d**(-4))<1e-13)
    check(f'chain d={d} normalized output order-independent',np.linalg.norm(second/np.linalg.norm(second)-target)<1e-13)
    bell_stats.append({'d':d,'two_joins_branch_amplitude':float(d**(-2)),'ordinary_born_selected_prob':float(d**(-4)),
          'order_residual':eps_order,'final_Bell_output_residual':eps_result})

# Invariance of the sequential *conditional* output does not equal a
# canonical probability distribution over the two unobserved descriptions.
# A putative scheduler can put any q ∈ [0,1] on the two indistinguishable orders.
for q in (0.,.2,.5,.9,1.):
    mixed_q=q*second+(1-q)*second_reverse
    check(f'arbitrary schedule q={q} yields same selected state',np.linalg.norm(mixed_q-second)<1e-13)

result={
    'status':'PASS', 'test_count':len(checks), 'checks':checks,
    'source_record_character_sectors':[{'weight':label,'X':xx,'Y':yy,'multiplicity':n} for label,xx,yy,n,w in sectors],
    'character_traces':{'TrX':82,'TrY':70,'TrXY':-8},
    'alpha_candidate_exact':str(alpha),'alpha_candidate_float':float(alpha),
    'auxiliary_exchange_superoperator_mode_multiplicities':mode_counts,
    'auxiliary_exchange_decay_eigenvalues':{'fixed':'0','X_only':'-2 gamma_X','Y_only':'-2 gamma_Y','both':'-2 (gamma_X+gamma_Y)'},
    'auxiliary_channel_trace_error_max':max_trace_err,
    'auxiliary_channel_semigroup_error_max':max_semigroup_err,
    'bell_chain_disjoint_joins':bell_stats,
    'conclusions':{
      'source_derived':'X,Y are commuting involutive exchanges labeling signed records; there is no time order for the joined reads; Born weights are not frequencies in the source ontology; one copy is read once',
      'exact_mathematics':'The char sectors allow continuum many auxiliary CP-TP rate assignments; Pauli read-order signs do not distinguish standalone channels; two disjoint Bell joins have identical selected final state in either execution order',
      'not_derived':'The original axiom does not choose a scheduler, stochastic law, an intrinsic clock, an electromagnetic propagator, or a QED beta function'
    }
}
OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'PASS','checks':len(checks),'alpha':float(alpha),'mode_counts':mode_counts,'bell_join_d':list(d['d'] for d in bell_stats),'max_semigroup_error':max_semigroup_err},indent=2))