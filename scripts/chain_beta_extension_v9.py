#!/usr/bin/env python3
"""Chain Theory beta program v9: explicit 172-of-175 record gauge descent.

Constructs paper S7's 169 two-slot records and three qubit commutator
*difference* records in a 175D ambient read-operator space (169 products
plus six one-slot qubit letters). Performs exact block-algebra and numerical
operator-basis checks for matching-group gauge actions.

This is a record-geometry test, NOT a derivation of QED's beta function.
"""
import json
from pathlib import Path
import numpy as np
import sympy as sp

DIMS=(3,2,3,2)
SLOTS=(8,3,8,3)
S3=3**.5; S2=2**.5; S6=6**.5


def gellmann():
    out=[]
    for i,j in ((0,1),(0,2),(1,2)):
        a=np.zeros((3,3),complex); a[i,j]=a[j,i]=1;out.append(a/np.sqrt(2))
        a=np.zeros((3,3),complex); a[i,j]=-1j; a[j,i]=1j;out.append(a/np.sqrt(2))
    out.extend((np.diag([1,-1,0]).astype(complex)/np.sqrt(2),np.diag([1,1,-2]).astype(complex)/np.sqrt(6)))
    return out


def pauli():
    return [np.array([[0,1],[1,0]],complex)/np.sqrt(2),
            np.array([[0,-1j],[1j,0]],complex)/np.sqrt(2),
            np.diag([1,-1]).astype(complex)/np.sqrt(2)]


def kron_all(ops):
    out=ops[0]
    for op in ops[1:]:out=np.kron(out,op)
    return out


def make_basis():
    letters=(gellmann(),pauli(),gellmann(),pauli())
    idents=[np.eye(d,dtype=complex)/np.sqrt(d) for d in DIMS]
    records=[];names=[];blocks={}
    pair_order=[('QQ',0,2),('qq',1,3),('Q1q1',0,1),('Q1q2',0,3),('Q2q1',2,1),('Q2q2',2,3)]
    for name,a,b in pair_order:
        inds=[]
        for i in range(SLOTS[a]):
            for j in range(SLOTS[b]):
                x=list(idents); x[a]=letters[a][i];x[b]=letters[b][j]
                inds.append(len(records));names.append((name,i,j));records.append(kron_all(x))
        blocks[name]=inds
    assert len(records)==169
    singles={}
    for slot,name in [(1,'q1'),(3,'q2')]:
        inds=[]
        for j in range(3):
            x=list(idents); x[slot]=letters[slot][j]
            inds.append(len(records));names.append((name,j));records.append(kron_all(x))
        singles[name]=inds
    assert len(records)==175
    # These operators are a set of 175 orthonormal directions in the
    # 1296-dimensional operator space of two (3,2) read arms.
    vectors=np.array(records).reshape(175,-1)
    gram=vectors.conj()@vectors.T
    assert np.max(np.abs(gram-np.eye(175)))<1e-13
    return records,blocks,singles,names


def swap_weight(M,idx,n):
    # X=-S on these same-end product words; symmetric swap => sqrt(3).
    start=idx[0]
    assert idx==list(range(start,start+n*n))
    for i in range(n):
        for j in range(n):
            a=start+n*i+j;b=start+n*j+i
            M[a,a]+=(S3+1)/2
            M[a,b]+=(S3-1)/2


def record_matrix(blocks,singles):
    M=np.zeros((175,175),float)
    swap_weight(M,blocks['QQ'],8)
    swap_weight(M,blocks['qq'],3)
    # Y=-swap: symmetric combinations of the q1/q2 slots get sqrt(2).
    for a,b in [('Q1q1','Q1q2'),('Q2q1','Q2q2')]:
        for i,j in zip(blocks[a],blocks[b]):
            M[i,i]=M[j,j]=(S2+1)/2
            M[i,j]=M[j,i]=(S2-1)/2
    q1,q2=singles['q1'],singles['q2']
    for i,j in zip(q1,q2):
        M[i,i]=M[j,j]=S6/2
        M[i,j]=M[j,i]=-S6/2
    assert np.max(np.abs(M-M.T))==0
    assert abs(np.trace(M)-(79+48*S2+42*S3+3*S6))<2e-12
    ev=np.linalg.eigvalsh(M)
    for eigen,multi in [(0,3),(1,79),(S2,48),(S3,42),(S6,3)]:
        assert np.sum(np.abs(ev-eigen)<1e-9)==multi,(eigen,multi)
    P=np.eye(175)
    for i,j in zip(q1,q2):
        P[i,i]=P[j,j]=0.5
        P[i,j]=P[j,i]=-0.5
    assert np.linalg.norm(P@P-P)<1e-13 and abs(np.trace(P)-172)<1e-13
    assert np.linalg.norm(M-P@M@P)<1e-12
    return M,P


def adjoint(U,letters):
    return np.array([[np.trace(a.conj().T@U@b@U.conj().T).real for b in letters] for a in letters])


def transport(blocks,singles,uq1=None,uQ1=None,uq2=None,uQ2=None):
    uq1=np.eye(2) if uq1 is None else uq1
    uq2=np.eye(2) if uq2 is None else uq2
    uQ1=np.eye(3) if uQ1 is None else uQ1
    uQ2=np.eye(3) if uQ2 is None else uQ2
    a1=adjoint(uQ1,gellmann()); a2=adjoint(uQ2,gellmann())
    b1=adjoint(uq1,pauli()); b2=adjoint(uq2,pauli())
    W=np.zeros((175,175))
    for name,x in [('QQ',np.kron(a1,a2)),('qq',np.kron(b1,b2)),
                   ('Q1q1',np.kron(a1,b1)),('Q1q2',np.kron(a1,b2)),
                   ('Q2q1',np.kron(a2,b1)),('Q2q2',np.kron(a2,b2)),
                   ('q1',b1),('q2',b2)]:
        inds=blocks[name] if name in blocks else singles[name]
        W[np.ix_(inds,inds)]=x
    assert np.linalg.norm(W@W.T-np.eye(175))<1e-12
    return W


def verify_physical_action(records,W,uq1,uQ1,uq2,uQ2):
    U=kron_all([uQ1,uq1,uQ2,uq2]); assert U.shape==(36,36)
    # Compare exact conjugation on every one of 175 operator directions.
    vectors=np.array(records).reshape(175,-1)
    moved=np.array([U@o@U.conj().T for o in records]).reshape(175,-1)
    matrix=vectors.conj()@moved.T
    err=np.linalg.norm(matrix-W)
    assert err<3e-12,err
    return err


def exact_symbolic_commutator():
    I=sp.eye(3)
    H=sp.diag(-1,-1,1)
    Pm=sp.Rational(1,2)*sp.BlockMatrix([[I,-I],[-I,I]]).as_explicit()
    W=sp.diag(-1,-1,1,1,1,1)
    moved=W*Pm*W.T
    hs2=sp.trace((Pm-moved)**2)
    overlap=sp.trace(Pm*moved)
    assert hs2==4 and overlap==1
    assert sp.trace(Pm)==3
    assert sp.simplify(sp.trace(Pm*(sp.sqrt(6)*moved)*Pm)-sp.sqrt(6))==0
    plus=sp.Rational(1,2)*sp.BlockMatrix([[I,I],[I,I]]).as_explicit()
    assert Pm*plus==sp.zeros(6)
    assert Pm+plus==sp.eye(6)
    return dict(projector_shift_HS_squared=str(hs2),fixed_commutator_trace_after_shift=str(sp.sqrt(6)*overlap),
        missing_sector_rank=str(plus.trace()),exact=True)


def sample_rotations(blocks,singles,M,P):
    trial=[]
    for theta in [0.0,np.pi/4,np.pi/2,np.pi]:
        uq=np.diag([np.exp(1j*theta/2),np.exp(-1j*theta/2)])
        uQ=np.eye(3,dtype=complex);uQ[:2,:2]=uq.conj()
        W=transport(blocks,singles,uq1=uq,uQ1=uQ)
        delta=float((np.trace(M)-np.trace(P@W@M@W.T@P))/31104)
        expected=float(S6*(1-np.cos(theta))/31104)
        assert abs(delta-expected)<1e-15, (delta,expected)
        trial.append({'theta_over_pi':theta/np.pi,'delta_alpha_fixed_projection':delta,
                      'formula':'sqrt(6)*(1-cos(theta))/31104'})
    return trial


def random_matched_SU2(blocks,singles,P,seed=317,steps=20):
    rng=np.random.default_rng(seed)
    max_pred_err=0; max_ket_err=0;max_orth=0
    bell=np.zeros(6,complex);bell[0]=bell[3]=1/np.sqrt(2)
    for _ in range(steps):
        Z=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        q,r=np.linalg.qr(Z)
        q=q@np.diag(np.conjugate(np.diag(r)/np.abs(np.diag(r))))
        q=q/np.sqrt(np.linalg.det(q))
        U3=np.eye(3,dtype=complex);U3[:2,:2]=np.conj(q)
        max_ket_err=max(max_ket_err,float(np.linalg.norm(np.kron(U3,q)@bell-bell)))
        W=transport(blocks,singles,uq1=q,uQ1=U3)
        R=adjoint(q,pauli())
        expected=float(3-np.trace(R))
        observed=float(np.linalg.norm(P-W@P@W.T)**2)
        max_pred_err=max(max_pred_err,abs(observed-expected))
        max_orth=max(max_orth,float(np.linalg.norm(R@R.T-np.eye(3))))
    assert max_ket_err<1e-13 and max_pred_err<1e-12
    return {'trials':steps,'max_pair_ket_invariance_error':max_ket_err,
            'max_leakage_formula_error':max_pred_err,'max_SO3_adjoint_orthogonality_error':max_orth}


def main():
    records,blocks,singles,names=make_basis()
    M,P=record_matrix(blocks,singles)
    id2=np.eye(2, dtype=complex);id3=np.eye(3, dtype=complex)
    # This paired matching leaves |omega_2>_{Q1 q1} EXACTLY fixed.
    uq1=np.diag([1j,-1j]);uQ1=np.diag([-1j,1j,1])
    assert abs(np.linalg.det(uq1)-1)<1e-14
    assert abs(np.linalg.det(uQ1)-1)<1e-14
    bell=np.zeros(6,complex);bell[0]=bell[3]=1/np.sqrt(2)
    match=np.kron(uQ1,uq1)
    assert np.linalg.norm(match@bell-bell)<1e-14
    W=transport(blocks,singles,uq1=uq1,uQ1=uQ1)
    phys_err=verify_physical_action(records,W,uq1,uQ1,id2,id3)
    WP=W@P@W.T
    gauge_leak_sq=float(np.linalg.norm(WP-P)**2)
    assert abs(gauge_leak_sq-4)<1e-11,gauge_leak_sq
    # Two commutator differences map exactly into the missing symmetric sums.
    q1,q2=singles['q1'],singles['q2']
    original_diffs=[(np.eye(175)[i]-np.eye(175)[j])/np.sqrt(2) for i,j in zip(q1,q2)]
    projected_l2=[float(np.linalg.norm(P@W@r)**2) for r in original_diffs]
    assert np.linalg.norm(np.array(projected_l2)-[0,0,1])<1e-12,projected_l2
    assert abs(np.trace(W@M@W.T)-np.trace(M))<2e-12
    projected_readout=float(np.trace(P@W@M@W.T@P))
    untransformed=float(np.trace(M))
    assert abs(projected_readout-(untransformed-2*S6))<2e-12
    # The 175D covariant closure remains the same under this exact matching.
    # Full independent SU(2) rotations on q1 and q2 force closure of 6 single-letter directions.
    rot_x=1j*np.array([[0,1],[1,0]],complex)
    rot_z=uq1
    sub=np.column_stack(original_diffs)
    space=[sub]
    for rot in [rot_x,rot_z]:
        gQ=(rot.conj())
        G=np.eye(3,dtype=complex);G[:2,:2]=gQ
        T=transport(blocks,singles,uq1=rot,uQ1=G)
        space.append(T@sub)
    closure_rank=int(np.linalg.matrix_rank(np.column_stack(space),tol=1e-10))
    assert closure_rank==6,closure_rank
    # Correlated independent read-frame transformations protect the full 172 record space,
    # as long as the same end-frame transformation acts on each read.
    Wcorr=transport(blocks,singles,uq1=uq1,uQ1=uQ1,uq2=uq1,uQ2=uQ1)
    assert np.linalg.norm(Wcorr@P@Wcorr.T-P)<1e-12
    assert np.linalg.norm(Wcorr@M@Wcorr.T-M)<1e-12
    # A finite four-rotation local twirl kills the off-diagonal q1/q2
    # part of the three commutator projectors. Its unique invariant value
    # in this sector is sqrt(6)/2 on each of the SIX singles.
    paulis=[np.eye(2,dtype=complex),
            1j*np.sqrt(2)*pauli()[0],1j*np.sqrt(2)*pauli()[1],1j*np.sqrt(2)*pauli()[2]]
    Es=np.zeros((6,6),float)
    inds=q1+q2
    for t in paulis:
        b=adjoint(t,pauli())
        U=np.block([[b,np.zeros((3,3))],[np.zeros((3,3)),np.eye(3)]])
        Es+=U@M[np.ix_(inds,inds)]@U.T/4
    assert np.linalg.norm(Es-(S6/2)*np.eye(6))<1e-12
    alpha=(79+48*S2+42*S3+3*S6)/31104
    delta=2*S6/31104
    symbolic=exact_symbolic_commutator()
    rotations=sample_rotations(blocks,singles,M,P)
    random_matches=random_matched_SU2(blocks,singles,P)
    outputs={
        'symbolic_commutator_audit':symbolic,
        'continuous_matching_rotation_sample':rotations,
        'random_matched_SU2_tests':random_matches,
        'source':'ChainTheoryAlpha.pdf, supplemental S7 equations S3--S5',
        'status':'exact finite-dimensional record test, not a QED beta derivation',
        'ambient_dimension':175,'paper_record_dimension':172,
        'paper_record_eigenvalue_multiplicities':{'1':79,'sqrt2':48,'sqrt3':42,'sqrt6':3,'zero_extra':3},
        'original_readout_exact':'(79 + 48*sqrt(2) + 42*sqrt(3) + 3*sqrt(6))/31104',
        'original_readout_decimal':alpha,
        'record_subspace_projector_rank':round(np.trace(P)),
        'explicit_operator_HS_orthonormality':'verified for 175 operators on C3*C2*C3*C2',
        'matched_local_SU3_SU2_exactly_preserves_link_ket':True,
        'matching_action_on_175_directions_residual':phys_err,
        'independent_matching_record_projector_HS_squared':gauge_leak_sq,
        'commutator_differences_survival_projections':projected_l2,
        'minimum_single_letter_gauge_closure_dimension':closure_rank,
        'full_ambient_invariant_trace':untransformed,
        'fixed_172_projected_trace_after_matching':projected_readout,
        'fixed_172_trace_difference_exact':'2*sqrt(6)',
        'fixed_172_normalized_readout_difference_exact':'2*sqrt(6)/31104',
        'fixed_172_normalized_readout_difference_decimal':delta,
        'correlated_identical_matching_M_and_record_projector_invariant':True,
        'commutator_sector_four_rotation_twirl':'sqrt(6)/2 * identity_6',
        'key_qualification':'All scalar trace shifts are for a fixed-coordinate re-projection prescription, NOT gauge-invariant physical predictions. Covariantly transported M has unchanged trace.'
    }
    out=Path('/mnt/data/chain_beta_extension_v9_results.json')
    out.write_text(json.dumps(outputs,indent=2,sort_keys=True)+'\n')
    print(json.dumps(outputs,indent=2,sort_keys=True))

if __name__=='__main__':main()