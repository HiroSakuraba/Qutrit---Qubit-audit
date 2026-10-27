#!/usr/bin/env python3
"""Chain Theory beta research v11.

An explicit, *source-motivated* 172-direction embedding in the complete
24 x 36 x 36 = 31104D operator Hilbert space, and a joint gauge action audit.

This is NOT an independent derivation of a QED beta function or a proof of
uniqueness of the physical representation in Wissner-Gross's axioms.

Requirements: numpy, sympy and chain_beta_extension_v9.py in this directory.
"""
import json
from pathlib import Path
import numpy as np
import sympy as sp
from chain_beta_extension_v9 import (
    make_basis, record_matrix, pauli, gellmann, adjoint, transport,
)

ROOT=Path(__file__).resolve().parent
sqrt=np.sqrt

def exact_connected():
    w=sp.zeros(6,1)
    w[0]=w[3]=sp.sqrt(sp.Rational(1,2))
    rho=w*w.T
    rho3=sp.Matrix(3,3,lambda a,b:sum(rho[2*a+i,2*b+i] for i in range(2)))
    rho2=sp.Matrix(2,2,lambda i,j:sum(rho[2*a+i,2*a+j] for a in range(3)))
    C=sp.simplify(rho-sp.kronecker_product(rho3,sp.eye(2)/2)
                  -sp.kronecker_product(sp.eye(3)/3,rho2)
                  +sp.eye(6)/6)
    assert sp.trace(C*C)==sp.Rational(3,4)
    assert sp.trace(C)==0
    # All zero partial traces prove C sits in su(3)\otimes su(2).
    assert all(sum(C[2*a+i,2*b+i] for i in range(2))==0
               for a in range(3) for b in range(3))
    assert all(sum(C[2*a+i,2*a+j] for a in range(3))==0
               for i in range(2) for j in range(2))
    return rho, C

def build_canonical_records(blocks,singles):
    """Write all 172 normalized read-arm records in the v9 175D basis.

    Matrices in the 175D basis are valid operator Hilbert space directions,
    not 175 distinct physical Hilbert-space basis states.
    """
    vectors=[]; weights=[]; tags=[]
    def one(items,fac,tag):
        x=np.zeros(175)
        for index,coeff in items:x[index]=coeff
        vectors.append(x);weights.append(float(fac));tags.append(tag)
    for name,n in [('QQ',8),('qq',3)]:
        inds=blocks[name]
        for i in range(n):
            one([(inds[i*n+i],1)],sqrt(3),name+'_diag')
        for i in range(n):
            for j in range(i+1,n):
                one([(inds[i*n+j],1/sqrt(2)),(inds[j*n+i],1/sqrt(2))],sqrt(3),name+'_sym')
                one([(inds[i*n+j],1/sqrt(2)),(inds[j*n+i],-1/sqrt(2))],1,name+'_anti')
    for left,right in [('Q1q1','Q1q2'),('Q2q1','Q2q2')]:
        for li,ri in zip(blocks[left],blocks[right]):
            one([(li,1/sqrt(2)),(ri,1/sqrt(2))],sqrt(2),'mixed_sym')
            one([(li,1/sqrt(2)),(ri,-1/sqrt(2))],1,'mixed_anti')
    for li,ri in zip(singles['q1'],singles['q2']):
        one([(li,1/sqrt(2)),(ri,-1/sqrt(2))],sqrt(6),'comm')
    E=np.array(vectors).T
    values=np.array(weights)
    assert E.shape==(175,172)
    assert np.max(np.abs(E.T@E-np.eye(172)))<1e-14
    from collections import Counter
    assert Counter(tags)==Counter({'QQ_diag':8,'QQ_sym':28,'QQ_anti':28,
      'qq_diag':3,'qq_sym':3,'qq_anti':3,
      'mixed_sym':48,'mixed_anti':48,'comm':3})
    return E,values,tags

def qutrit_match(u):
    G=np.eye(3,dtype=complex)
    G[:2,:2]=np.conjugate(u)
    return G

def bell3():
    b=np.zeros(9,complex)
    for i in range(3):b[3*i+i]=1/sqrt(3)
    return b

def scalar_delta(s1,s2):
    # 2-norm squared of two rank-1 projectors, avoiding an extra 9x9 matrix.
    return float(2*(1-abs(np.vdot(s1,s2))**2))

def main():
    rho_sym, C_sym=exact_connected()
    C=np.array(C_sym, dtype=complex)
    cn=C/sqrt(0.75)
    bases=[np.kron(a,b) for a in gellmann() for b in pauli()]
    cvec=np.array([np.trace(b.conj().T@cn) for b in bases],complex)
    reconstructed=sum(z*b for z,b in zip(cvec,bases))
    assert np.linalg.norm(reconstructed-cn)<5e-15
    assert abs(np.vdot(cvec,cvec)-1)<5e-15

    arms,blocks,singles,_=make_basis()
    M,P=record_matrix(blocks,singles)
    E,weights,tags=build_canonical_records(blocks,singles)
    assert np.linalg.norm(M@E-E*weights)<2e-13
    assert np.linalg.norm(E@E.T-P)<5e-13
    arm_op=np.stack([A.reshape(-1) for A in arms],axis=1)  # 1296x175
    assert np.linalg.norm(arm_op.conj().T@arm_op-np.eye(175))<2e-13
    read_dirs=arm_op@E # 1296x172
    # An explicit 31104x172 array: each column is a concrete operator-HS
    # direction in R=(traceless3 tensor traceless2) tensor End(C6) tensor End(C6).
    full_dirs=np.einsum('u,vj->uvj',cvec,read_dirs).reshape(31104,172)
    full_norms=np.sum(abs(full_dirs)**2,axis=0)
    assert np.max(abs(full_norms-1))<5e-14
    test_inds=np.array([0,8,35,70,100,150,168,171])
    assert np.linalg.norm(full_dirs[:,test_inds].conj().T@full_dirs[:,test_inds]
         -np.eye(len(test_inds)))<1e-13
    assert 24*36*36==31104
    expected_trace=79+48*sqrt(2)+42*sqrt(3)+3*sqrt(6)
    assert abs(np.sum(weights)-expected_trace)<2e-13
    assert sum(tags.count(t) for t in set(tags))==172

    # Explicit full-tensor direction for single two-qubit commutator record.
    # The ambient matrix itself need not be formed: M=full_dirs diag(w) full_dirs^H.
    # Trace and nonzero eigenvalues follow from orthonormality of columns.

    # Gauge covariance for unread connected factor: G3 \otimes u3.
    phi=np.array([0.37,-0.13,-0.24])
    G3=np.diag(np.exp(1j*phi))
    u3=np.diag(np.exp(1j*np.array([0.26,-0.26])))
    U3=np.kron(G3,u3)
    Cmove=U3@cn@U3.conj().T
    cvec_move=np.array([np.trace(b.conj().T@Cmove) for b in bases])
    Wu=np.kron(adjoint(G3,gellmann()),adjoint(u3,pauli()))
    unread_cov=float(np.linalg.norm(cvec_move-Wu@cvec))
    assert unread_cov<2e-14
    u_pair=np.diag([1j,-1j]); G_pair=qutrit_match(u_pair)
    assert np.linalg.norm(np.kron(G_pair,u_pair)@np.array([1,0,0,1,0,0])/sqrt(2)
                         -np.array([1,0,0,1,0,0])/sqrt(2))<2e-14
    assert np.linalg.norm(np.kron(G_pair,u_pair)@cn@np.kron(G_pair,u_pair).conj().T-cn)<2e-14

    W=transport(blocks,singles,uq1=u_pair,uQ1=G_pair)
    full_subspace_shift=float(np.linalg.norm(P-W@P@W.T)**2)
    assert abs(full_subspace_shift-4)<2e-13
    readout_delta=float((np.trace(M)-np.trace(P@W@M@W.T@P))/31104)
    assert abs(readout_delta-2*sqrt(6)/31104)<2e-15
    weighted_trace_cov=float(abs(np.trace(M)-np.trace(W@M@W.T)))
    assert weighted_trace_cov<2e-13

    # The two distinct correlations protect *different* structures:
    # Same qubit rotations: selected record space invariant, Bell join not.
    # Conjugate qubit rotations: Bell join invariant, record differences not.
    join=bell3()
    u=np.diag(np.exp(1j*np.array([np.pi/2,-np.pi/2])))  # theta = pi
    A=qutrit_match(u)
    Wdiag=transport(blocks,singles,uq1=u,uQ1=A,uq2=u,uQ2=A)
    diag_record_shift=float(np.linalg.norm(P-Wdiag@P@Wdiag.T)**2)
    diag_join=np.kron(A,A)@join
    diag_join_fid=float(abs(np.vdot(join,diag_join))**2)
    assert diag_record_shift<2e-13
    assert abs(diag_join_fid-sp.Rational(1,9))<2e-14
    # Independent exact SymPy witness for the two stabilizer scenarios.
    zz=(1+sp.I)/sp.sqrt(2)
    uu=sp.diag(zz,sp.conjugate(zz))
    uubar=sp.conjugate(uu)
    GA=sp.diag(sp.conjugate(zz),zz,1)
    GB=sp.diag(zz,sp.conjugate(zz),1)
    B3=sp.zeros(9,1)
    for ii in range(3): B3[3*ii+ii]=1/sp.sqrt(3)
    B2=sp.Matrix([1,0,0,1])/sp.sqrt(2)
    assert all(sp.simplify(t)==0 for t in sp.kronecker_product(GA,GB)*B3-B3)
    assert all(sp.simplify(t)==0 for t in sp.kronecker_product(uu,uubar)*B2-B2)
    R1=sp.Matrix([[0,-1,0],[1,0,0],[0,0,1]])
    R2=R1.T
    W6=sp.diag(R1,R2)
    I3=sp.eye(3)
    P6=sp.Rational(1,2)*sp.BlockMatrix([[I3,-I3],[-I3,I3]]).as_explicit()
    assert sp.trace((P6-W6*P6*W6.T)**2)==4
    Gdiag=sp.diag(-sp.I,sp.I,1)
    assert sp.simplify((B3.T*sp.kronecker_product(Gdiag,Gdiag)*B3)[0]+sp.Rational(1,3))==0

    u=np.diag(np.exp(1j*np.array([np.pi/4,-np.pi/4])))  # theta = pi/2
    Ac=qutrit_match(u)
    uc=np.conjugate(u)
    Ac2=qutrit_match(uc)
    Wanti=transport(blocks,singles,uq1=u,uQ1=Ac,uq2=uc,uQ2=Ac2)
    anti_join_fid=float(abs(np.vdot(join,np.kron(Ac,Ac2)@join))**2)
    anti_record_shift=float(np.linalg.norm(P-Wanti@P@Wanti.T)**2)
    assert abs(anti_join_fid-1)<2e-14
    assert abs(anti_record_shift-4)<2e-13
    # Verify the ENTIRE postselected source-inspired three-pair state:
    # |omega3>_{Q1Q2} tensor |omega2>_{q1q2} tensor |omega2>_{Q3f}.
    b2=np.array([1,0,0,1],complex)/sqrt(2)
    unread_bell=np.array([1,0,0,1,0,0],complex)/sqrt(2)
    Omega=np.kron(np.kron(join,b2),unread_bell)
    Uomega=np.kron(np.kron(np.kron(Ac,Ac2),np.kron(u,uc)),np.eye(6))
    joint_stabilizer_error=float(np.linalg.norm(Uomega@Omega-Omega))
    assert Omega.shape==(216,) and joint_stabilizer_error<2e-14
    # Original selected record subspace still changes despite full Omega fixed.
    assert anti_record_shift>3.999999999

    # For a real SU(2) matrix u, u*=u, both relations coincide (small
    # common subgroup). This does not select a preferred relationship.
    theta=0.419
    u_real=np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
    Ar=qutrit_match(u_real)
    Wreal=transport(blocks,singles,uq1=u_real,uQ1=Ar,uq2=u_real,uQ2=Ar)
    assert np.linalg.norm(P-Wreal@P@Wreal.T)<2e-13
    assert scalar_delta(join,np.kron(Ar,Ar)@join)<2e-13

    # Direct operator covariance on full R for a representative subset; no
    # dense 31104x31104 matrix is formed.
    Uread=np.kron(np.kron(np.kron(G_pair,u_pair),np.eye(3)),np.eye(2))
    physical_err=max(np.linalg.norm((Uread@op@Uread.conj().T).reshape(-1)-arm_op@(W@np.eye(175)[:,j]))
                     for j,op in enumerate(arms))
    assert physical_err<2e-13
    # All three virtual-copy arms: explicit covariance of the FULL 31104D
    # record directions, using direct conjugation on a sample of 12 records.
    sampled=[0,3,8,27,44,63,70,85,108,135,168,171]
    full_covariance_err=0.
    for j in sampled:
        readop=read_dirs[:,j].reshape((36,36))
        physical_read=(Uread@readop@Uread.conj().T).reshape(-1)
        physical_full=np.kron(cvec_move,physical_read)
        predicted_full=np.kron(Wu@cvec, arm_op@W@E[:,j])
        full_covariance_err=max(full_covariance_err,float(np.linalg.norm(physical_full-predicted_full)))
    assert full_covariance_err<3e-13

    data={
      'description':'v11 source-motivated full-31104 operator-HS embedding and gauge-action audit',
      'source_verification':{
        'source_unread_link_dimension':24,'source_read_arm_dimension':1296,
        'full_readout_dimension':31104,'weighted_record_rank':172,
        'record_multiplicities':{'1':79,'sqrt(2)':48,'sqrt(3)':42,'sqrt(6)':3},
        'weighted_trace':float(expected_trace),
        'alpha_candidate':float(expected_trace/31104),
        'full_record_array_shape':list(full_dirs.shape),
        'max_full_column_norm_error':float(np.max(abs(full_norms-1))),
        'subset_full_gram_error':float(np.linalg.norm(
          full_dirs[:,test_inds].conj().T@full_dirs[:,test_inds]-np.eye(len(test_inds))))
      },
      'connected_unread':{'tr_C_squared_exact':'3/4','C_entries_exact':str(C_sym),
        'coordinate_reconstruction_error':float(np.linalg.norm(reconstructed-cn)),
        'passive_matched_frame_covariance_error':unread_cov},
      'matchings':{'independent_single_read_matched_link_stabilizer_projector_HS_squared':full_subspace_shift,
        'fixed_subspace_reprojection_readout_loss':readout_delta,
        'weighted_trace_covariant_error':weighted_trace_cov,
        'same_read_frames':{'record_projector_HS_squared':diag_record_shift,
                            'qutrit_join_fidelity':diag_join_fid},
        'conjugate_read_frames':{'qutrit_join_fidelity':anti_join_fid,
                                 'record_projector_HS_squared':anti_record_shift,
                                 'entire_joint_state_stabilizer_error':joint_stabilizer_error},
        'shared_real_SU2_subgroup_preserves_both':True,
        'exact_sympy_join_stabilizer_and_commutator_displacement':True,
        'sample_physical_operator_covariance_error':float(physical_err),
        'full_31104D_three_arm_covariance_error_12_records':full_covariance_err},
      'scope_limits':[
        'The 31104D embedding is a concrete source-motivated realization, not a proof the original axiom selects it uniquely.',
        'Invariance of trace is weaker than gauge descent of individual record operators.',
        'Fixed join and fixed commutator-difference projector have different symmetry stabilizers.',
        'Matching-basis dependence is not alone a source contradiction: S7 places matching invariance on reads rather than fixed individual records.',
        'No mass, momentum, spacetime evolution, QED beta function, or physical Wilson loop was derived.'
      ]
    }
    out=ROOT/'chain_beta_extension_v11_results.json'
    out.write_text(json.dumps(data,indent=2,default=str)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k!='connected_unread'},indent=2))
    print('PASS — v11 exact and numerical checks, results:',out)

if __name__=='__main__':main()