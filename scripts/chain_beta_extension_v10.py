#!/usr/bin/env python3
"""Chain Theory -> QED: v10 three-copy frame-transport and holonomy audit.

Primary-source inputs: ChainTheoryAlpha.pdf, Supplemental Sections S1/S7/S8.
Extends the explicit 175D ambient record embedding of v9, clearly distinguished
from a derivation of gauge dynamics or from an operational definition of holonomy.
"""
import importlib.util
import json
from pathlib import Path
import numpy as np
import sympy as sp

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('v9',HERE/'chain_beta_extension_v9.py')
v9=importlib.util.module_from_spec(spec);spec.loader.exec_module(v9)

I2=np.eye(2,dtype=complex);I3=np.eye(3,dtype=complex)
sigx=np.array([[0,1],[1,0]],complex)
sigy=np.array([[0,-1j],[1j,0]],complex)
sigz=np.diag([1,-1]).astype(complex)


def u2_rotation(theta, axis):
    axis=np.asarray(axis,dtype=float);axis=axis/np.linalg.norm(axis)
    mat=sum((a*s for a,s in zip(axis,[sigx,sigy,sigz])),np.zeros((2,2),complex))
    return np.cos(theta/2)*I2-1j*np.sin(theta/2)*mat


def su2_haar(rng):
    z=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
    q,r=np.linalg.qr(z)
    q=q@np.diag(np.conjugate(np.diag(r)/np.abs(np.diag(r))))
    return q / np.sqrt(np.linalg.det(q))


def matched_SU3(u):
    # For SU2 u, pair-preserving qutrit operation is conjugate u on occupied 2D support.
    return np.block([[np.conj(u),np.zeros((2,1))],[np.zeros((1,2)),np.ones((1,1))]])


def bell32():
    z=np.zeros(6,complex);z[0]=z[3]=1/np.sqrt(2)
    return z


def three_pair_state():
    return np.kron(np.kron(bell32(),bell32()),bell32())


def adjoint3(u):
    return v9.adjoint(u,v9.pauli())


def hs_projector(t):
    """Difference plane {(v,-t.T v)} in two adjoint-3 read frames."""
    return np.block([[I3.real,-t],[-t.T,I3.real]])/2


def two_copy_matching_test():
    rot=u2_rotation(np.pi/2,[0,0,1]); m3=matched_SU3(rot)
    ps=bell32()
    assert np.linalg.norm(np.kron(m3,rot)@ps-ps)<1e-14
    assert abs(np.linalg.det(m3)-1)<1e-14
    # A physically identical local pair admits a nontrivial read-frame stabilizer.
    # Hence a coordinate-fixed antisymmetric inter-read subspace is not invariant.
    r=adjoint3(rot);p=hs_projector(I3.real);pt=hs_projector(r)
    d=float(np.linalg.norm(p-pt)**2)
    assert abs(d-2)<1e-12, d
    # The projection is transported rather than left fixed.
    g=np.block([[r,np.zeros((3,3))],[np.zeros((3,3)),np.eye(3)]])
    assert np.linalg.norm(g@p@g.T-pt)<1e-13
    return {'local_pair_ket_stabilized':True,'fixed_record_projector_HS_distance_squared':d,
            'transported_record_projector_rank':int(round(np.trace(pt)))}


def exact_sympy_test():
    eye=sp.eye(3); rz=sp.Matrix([[0,-1,0],[1,0,0],[0,0,1]])
    def pd(t): return sp.BlockMatrix([[eye,-t],[-t.T,eye]]).as_explicit()/2
    p=pd(eye);pt=pd(rz)
    assert sp.simplify(p*p-p)==sp.zeros(6)
    assert sp.trace(pt)==3
    assert sp.trace((p-pt)**2)==2
    assert sp.trace(sp.sqrt(6)*pt)==3*sp.sqrt(6)
    assert sp.trace(p*sp.sqrt(6)*pt*p)==2*sp.sqrt(6) # angle pi/2 -> 2+cos theta
    # Exact SU(2) Wilson loops, flat versus a pi/2 twist.
    z=sp.diag((1+sp.I)/sp.sqrt(2),(1-sp.I)/sp.sqrt(2))
    assert sp.simplify(sp.trace(z)-sp.sqrt(2))==0
    assert sp.trace(rz)==1
    return {'projector_rank':3,'projector_distance_squared':2,
            'weighted_trace':str(3*sp.sqrt(6)),
            'fixed_space_weighted_trace':str(2*sp.sqrt(6)),
            'SU2_twist_wilson_trace':str(sp.sqrt(2)),
            'SO3_twist_wilson_trace':1}


def three_copy_gauge_test():
    rng=np.random.default_rng(10110)
    ket=three_pair_state()
    assert ket.shape==(216,)
    max_state_error=0;max_flat_loop_error=0;max_transport_error=0
    max_loop_trace_error=0;max_nonflat_unchanged=0
    theta=np.pi/2
    twist=u2_rotation(theta,[0,0,1])
    ts0=[I2.copy(),I2.copy(),twist.copy()]
    h0=ts0[0]@ts0[1]@ts0[2]
    nonflat_wilson=float(np.real(np.trace(h0)))
    assert abs(nonflat_wilson-np.sqrt(2))<1e-14
    for _ in range(60):
        us=[su2_haar(rng) for _ in range(3)]
        a=[matched_SU3(u) for u in us]
        pair_ops=[np.kron(t,u) for t,u in zip(a,us)]
        moved=pair_ops[0]@bell32()
        for op in pair_ops[1:]: moved=np.kron(moved,op@bell32())
        max_state_error=max(max_state_error,float(np.linalg.norm(moved-ket)))
        # All transports of the form Ti<-j=gi gj^dagger have trivial holonomy.
        flat=[us[0]@us[1].conj().T,us[1]@us[2].conj().T,us[2]@us[0].conj().T]
        hf=flat[0]@flat[1]@flat[2]
        max_flat_loop_error=max(max_flat_loop_error,float(np.linalg.norm(hf-I2)))
        # An independent choice of connection transforms with the same vertex gauge action.
        gauged=[us[0]@ts0[0]@us[1].conj().T,
                us[1]@ts0[1]@us[2].conj().T,
                us[2]@ts0[2]@us[0].conj().T]
        hn=gauged[0]@gauged[1]@gauged[2]
        max_transport_error=max(max_transport_error,float(np.linalg.norm(hn-us[0]@h0@us[0].conj().T)))
        max_loop_trace_error=max(max_loop_trace_error,float(abs(np.trace(hn)-np.trace(h0))))
        max_nonflat_unchanged=max(max_nonflat_unchanged,float(abs(np.trace(hn)-nonflat_wilson)))
    assert max_state_error<3e-14
    assert max_flat_loop_error<2e-14
    assert max_transport_error<2e-14
    assert max_loop_trace_error<2e-14
    return {'trials':60,'max_three_pair_ket_invariance_error':max_state_error,
            'max_pure_gauge_holonomy_error':max_flat_loop_error,
            'max_twisted_holonomy_covariance_error':max_transport_error,
            'max_twisted_wilson_trace_error':max_loop_trace_error,
            'flat_SU2_wilson_trace':2.,'twisted_SU2_wilson_trace':nonflat_wilson,
            'flat_SO3_wilson_trace':3.,'twisted_SO3_wilson_trace':1.,
            'caveat':'Triangle transports are added relational labels, not links specified in the source three-arm object.'}


def pair_diff_transport_test():
    rng=np.random.default_rng(2241);max_cov=0;max_flat=0;max_trace=0
    for _ in range(70):
        u1,u2=(su2_haar(rng),su2_haar(rng))
        r1,r2=(adjoint3(u1),adjoint3(u2))
        t=u2_rotation(0.8,[1,1,2]);p=adjoint3(t)
        q=r1@p@r2.T
        d=np.block([[r1,np.zeros((3,3))],[np.zeros((3,3)),r2]])
        pa=hs_projector(p);pb=hs_projector(q)
        max_cov=max(max_cov,float(np.linalg.norm(d@pa@d.T-pb)))
        max_flat=max(max_flat,float(np.linalg.norm(pb@pb-pb)))
        max_trace=max(max_trace,abs(float(np.trace(pb)-3)))
    assert max_cov<2e-13 and max_flat<2e-13
    return {'trials':70,'max_projector_covariance_error':max_cov,
            'max_projector_idempotence_error':max_flat,
            'max_rank_trace_error':max_trace,
            'spectrum':'three eigenvalues 1, three 0',
            'weighted_trace':'3 sqrt(6) independent of relative transport'}


def source_record_unread_test():
    records,blocks,singles,names=v9.make_basis()
    m,p=v9.record_matrix(blocks,singles)
    assert len(records)==175
    rot1=u2_rotation(.41,[1,2,3]); rot2=u2_rotation(.83,[1,3,2]);rot3=u2_rotation(.97,[2,3,1]);
    uqs=[matched_SU3(u) for u in (rot1,rot2,rot3)]
    w=v9.transport(blocks,singles,uq1=rot1,uQ1=uqs[0],uq2=rot2,uQ2=uqs[1])
    mt=w@m@w.T;pt=w@p@w.T
    assert np.linalg.norm(mt-pt@mt@pt)<2e-12
    alpha=float(np.trace(m)/31104)
    tr_err=float(abs(np.trace(mt)-np.trace(m)))
    assert tr_err<1e-11
    # Connected/traceless unread-link factor has dimension 8x3=24.
    # A single fixed unit vector in that factor produces an illustrative
    # 172-direction embedding in the readout, NOT a full reconstruction
    # of how the source assigns the unread component to each record.
    w3=np.kron(v9.adjoint(uqs[2],v9.gellmann()),v9.adjoint(rot3,v9.pauli()))
    assert w3.shape==(24,24)
    v=np.zeros(24);v[0]=1
    vprime=w3@v
    overlap=float(abs(v@vprime)**2)
    gram_shift_squared=float(2*319*(1-overlap))
    # 319 is sum of squared 172 weights: 79+48*2+42*3+3*6.
    assert abs(np.trace(m@m)-319)<1e-10
    assert abs(np.linalg.norm(vprime)-1)<1e-13
    assert gram_shift_squared>=0
    return {'original_record_count':int(round(np.trace(p))),
            'weighted_trace':float(np.trace(m)), 'readout_alpha':alpha,
            'two_read_arm_covariant_trace_error':tr_err,
            'unread_connected_dimension':24,
            'unread_test_factor_overlap_squared':overlap,
            'illustrative_full_record_operator_HS_shift_squared_due_only_to_unread_matching':gram_shift_squared,
            'scope':'Unchanged 172-record spectrum under transported two-arm action; unread-link result is a model embedding, not proof of gauge descent for the source exact readout.'}


def rank_deficiency_test():
    # Overlap of two rank-2 qutrit supports is a contraction, sometimes singular.
    V1=np.array([[1,0],[0,1],[0,0]],complex)
    V2=np.array([[0,0],[1,0],[0,1]],complex)
    overlap=V1.conj().T@V2
    assert np.linalg.matrix_rank(overlap)==1
    assert np.allclose(np.linalg.svd(overlap,compute_uv=False),[1,0])
    # The polar partial isometry does not give a unique unitary extension:
    # U_phi=[[0,exp(i phi)],[1,0]], all act identically on supported column.
    extensions=[]
    for phi in [0,np.pi/3,np.pi/2,np.pi]:
        U=np.array([[0,np.exp(1j*phi)],[1,0]],complex)
        assert np.linalg.norm(U.conj().T@U-I2)<1e-14
        assert np.linalg.norm(U[:,0]-np.array([0,1]))<1e-14
        extensions.append(U)
    assert np.linalg.norm(extensions[0]-extensions[1])>.1
    return {'intersection_support_dimension':1,
            'overlap_rank':1,'singular_values':[1.,0.],
            'distinct_unitary_extensions_tested':len(extensions),
            'result':'Support-overlap polar transport is not a uniquely completed U(2) connection at rank-deficient overlaps.'}


def main():
    out={'program':'Chain Theory QED beta, v10: relational three-copy transport / holonomy',
         'date':'2026-10-24','source':'ChainTheoryAlpha.pdf supplemental S1 S7 S8',
         'tests':{
            'exact_symbolic':exact_sympy_test(),
            'link_stabilizer':two_copy_matching_test(),
            'three_copy_and_holonomy':three_copy_gauge_test(),
            'covariant_difference_projector':pair_diff_transport_test(),
            'full_172_and_unread_link':source_record_unread_test(),
            'partial_isometry_completion':rank_deficiency_test()},
         'verdict':'The pair axioms allow a gauge-covariant family of record spaces once a connection is provided; they do not define or uniquely select a three-copy connection or its loop holonomy. No beta function follows.'}
    path=HERE/'chain_beta_extension_v10_results.json'
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()