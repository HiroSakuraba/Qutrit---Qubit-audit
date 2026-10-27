#!/usr/bin/env python3
"""Independent exact SymPy audit of v14 sector expectations for both source states.

Runs without floating-point matrices; uses the original per-slot construction and
172 sparse record directions. Depends on the v9/v11 indexing only.
"""
import json
from pathlib import Path
import sympy as S
from chain_beta_extension_v9 import make_basis
from chain_beta_extension_v11 import build_canonical_records

root=Path(__file__).resolve().parent
sq2=S.sqrt(2);sq3=S.sqrt(3);sq6=S.sqrt(6); I=S.I

def gens3():
    out=[]
    for a,b in ((0,1),(0,2),(1,2)):
        M=S.zeros(3);M[a,b]=M[b,a]=1/sq2;out.append(M)
        M=S.zeros(3);M[a,b]=-I/sq2;M[b,a]=I/sq2;out.append(M)
    out.extend([S.diag(1,-1,0)/sq2,S.diag(1,1,-2)/sq6])
    return out

def gens2():
    return [S.Matrix([[0,1],[1,0]])/sq2,
            S.Matrix([[0,-I],[I,0]])/sq2,
            S.diag(1,-1)/sq2]

B=(gens3(),gens2(),gens3(),gens2())
IDS=[S.eye(3)/sq3,S.eye(2)/sq2,S.eye(3)/sq3,S.eye(2)/sq2]
LOC={'QQ':(0,2),'qq':(1,3),'Q1q1':(0,1),'Q1q2':(0,3),
     'Q2q1':(2,1),'Q2q2':(2,3),'q1':(1,),'q2':(3,)}

def matrix_from_label(label):
    ops=list(IDS)
    slots=LOC[label[0]]
    for slot,index in zip(slots,label[1:]):ops[slot]=B[slot][index]
    return S.kronecker_product(*ops)

def states():
    b=S.zeros(6,1);b[0]=b[3]=1/sq2
    unjoined=S.kronecker_product(b,b)
    joined=S.zeros(36,1)
    for a in range(3):
        for b in range(2):joined[(((a*2+b)*3+a)*2+b),0]=1/sq6
    assert (unjoined.T*unjoined)[0]==1
    assert (joined.T*joined)[0]==1
    return unjoined,joined

def exact_coeff(x):
    if abs(x-1)<1e-12:return S.Integer(1)
    if abs(x+1)<1e-12:return S.Integer(-1)
    if abs(x-2**(-0.5))<1e-12:return 1/sq2
    if abs(x+2**(-0.5))<1e-12:return -1/sq2
    raise ValueError(x)

def main():
    _,blocks,singles,names=make_basis()
    records,w,tags=build_canonical_records(blocks,singles)
    psi_pre,psi_post=states()
    # Materialize each basis operator's action on 36-state vectors; not its full matrix.
    transformed=[[S.simplify(matrix_from_label(label)*psi) for label in names]
                for psi in (psi_pre,psi_post)]
    scores={k:[S.Rational(0),S.Rational(0)] for k in ('1','sqrt2','sqrt3','sqrt6')}
    inds={k:[] for k in scores}
    for j,t in enumerate(tags):
        key= ('sqrt6' if t=='comm' else 'sqrt2' if t=='mixed_sym' else
              'sqrt3' if t.endswith('_diag') or t.endswith('_sym') else '1')
        inds[key].append(j)
    for j in range(172):
        key=next(k for k in inds if j in inds[k])
        cs=[(i,exact_coeff(float(records[i,j]))) for i in range(175) if abs(records[i,j])>1e-12]
        for sidx in range(2):
            vs=sum((c*transformed[sidx][i] for i,c in cs),S.zeros(36,1))
            n2=S.simplify((vs.conjugate().T*vs)[0]);scores[key][sidx]+=n2
    # Cbar^2 expectation on unread |omega2> = 3/4; frame factor = 144.
    expected_pre={'1':S.Rational(103,9216),'sqrt2':S.Rational(1,144),
                  'sqrt3':S.Rational(175,27648),'sqrt6':S.Rational(1,2304)}
    expected_post={'1':S.Rational(55,6912),'sqrt2':S.Rational(1,108),
                   'sqrt3':S.Rational(25,3456),'sqrt6':S.Rational(1,3456)}
    result={}
    for si,(label,exp) in enumerate([('unjoined',expected_pre),('postselected',expected_post)]):
        probs={k:S.simplify(x[si]*S.Rational(3,4)/144) for k,x in scores.items()}
        assert all(S.simplify(probs[k]-exp[k])==0 for k in exp),(label,probs,exp)
        alpha=S.simplify(probs['1']+sq2*probs['sqrt2']+sq3*probs['sqrt3']+sq6*probs['sqrt6'])
        result[label]={'probabilities':{k:str(probs[k]) for k in probs},
                       'score_exact':str(alpha),'score_decimal':str(S.N(alpha,17))}
    print(json.dumps(result,indent=2))
    (root/'chain_beta_extension_v14_exact_results.json').write_text(json.dumps(result,indent=2))

if __name__=='__main__':main()