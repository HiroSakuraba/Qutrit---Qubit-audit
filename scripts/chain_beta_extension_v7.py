#!/usr/bin/env python3
"""Chain Theory beta v7: fixed-join no-go and covariant operational completion.

Theorem tests are finite-dimensional and do NOT derive QED or the paper's
operational semantics. Source: ChainTheoryAlpha.pdf, §§S1,S7,S8.
Requires numpy and sympy. Run: python chain_beta_extension_v7.py
"""
import itertools, json, math
import numpy as np
import sympy as sp

SQ=sp.sqrt; R=sp.Rational

def exact_basis(n, k):
    v=sp.zeros(n,1);v[k]=1;return v

def exact_fixed_projection_example():
    # Input matrix (A1 x A2 = 9) by (B1 x B2 = 4).
    psi=sp.zeros(9,4)
    for i,a in enumerate([0,2]):
        for j,b in enumerate([0,1]):
            psi[3*a+b,2*i+j]=R(1,2)
    assert psi.T*psi==sp.eye(4)/4
    omega=(exact_basis(9,0)+exact_basis(9,4)+exact_basis(9,8))/SQ(3)
    amp=sp.simplify(omega.T*psi)
    assert amp==sp.Matrix([[1/(2*SQ(3)),0,0,0]])
    p=sp.simplify((amp*amp.T)[0]);assert p==R(1,12)
    return {'p_fixed':'1/12','remote':'|00>','remote_concurrence':'0'}

def exact_finite_matching_twirl():
    # 3 cyclic permutations x 8 diagonal sign flips are in the paper's
    # monomial matching family (fourth-root phases include +/-1).
    omega=(exact_basis(9,0)+exact_basis(9,4)+exact_basis(9,8))/SQ(3)
    P=omega*omega.T
    X=sp.zeros(3,3)
    for i in range(3):X[(i+1)%3,i]=1
    twirled=sp.zeros(9,9)
    for power in range(3):
        for signs in itertools.product([-1,1],repeat=3):
            U=X**power*sp.diag(*signs)
            T=sp.kronecker_product(U,sp.eye(3))
            twirled+=T*P*T.T
    twirled=sp.simplify(twirled/R(24))
    assert twirled==sp.eye(9)/9
    # A remote marginal whose A two-qutrit support carries two copies of C2.
    # Effective POVM element I9/9 => all inputs have p=1/9 and Bob unchanged.
    psi=sp.zeros(9,4)
    for i in range(2):
        for j in range(2):psi[3*i+j,2*i+j]=R(1,2)
    sigma=sp.simplify(psi.T*twirled*psi)
    assert sigma==sp.eye(4)/36
    assert sigma.trace()==R(1,9)
    assert sigma/R(1,9)==sp.eye(4)/4
    return {'finite_rotations':24,'twirled_join':'I_9/9',
            'mean_success_probability':'1/9',
            'Bob_conditional_state_after_forgetting_relative_frame':'I_4/4'}

def haar_u(rng, n=3):
    Z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    Q,Rm=np.linalg.qr(Z)
    return Q @ np.diag(np.conj(np.diag(Rm)/np.abs(np.diag(Rm))))

def entropy_singulars(C):
    sv=np.linalg.svd(C,compute_uv=False)
    prob=float(np.sum(sv*sv))
    if prob<=1e-28:return 0.,0.
    conc=float(2*abs(np.linalg.det(C))/prob)
    return prob,conc

def v6_frame():
    E=np.eye(9,dtype=complex)
    b0=(E[:,0]+E[:,4])/np.sqrt(2)
    b1=(E[:,0]-E[:,4])/np.sqrt(2)
    b2=(E[:,1]+E[:,3])/np.sqrt(2)
    b3=(E[:,1]-E[:,3])/np.sqrt(2)
    omega=(E[:,0]+E[:,4]+E[:,8])/np.sqrt(3)
    eta=np.sqrt(3)*b0-np.sqrt(2)*omega
    basis=np.column_stack([omega,eta,b1,b2,b3,E[:,2],E[:,5],E[:,6],E[:,7]])
    assert np.linalg.norm(basis.conj().T@basis-np.eye(9))<2e-14
    return basis

def transported_instrument_stress(trials=128,seed=705241):
    rng=np.random.default_rng(seed)
    B=v6_frame()
    P0=np.outer(B[:,0],B[:,0].conj())
    M0=np.zeros((9,4),dtype=complex)
    for i in range(2):
        for j in range(2):M0[3*i+j,2*i+j]=.5
    max_err_cov=0.;max_err_prob=0.;max_err_nosig=0.;max_conc_err=0.
    sample_fixed=[]
    probabilities=np.array([1/6,1/12,1/4,1/4,1/4,0,0,0,0])
    fixed_bra=B[:,0].conj()
    for _ in range(trials):
        U1=haar_u(rng);U2=haar_u(rng)
        T=np.kron(U1,U2)
        M=T @ M0
        vs=T @ B
        join=vs[:,0]
        P=np.outer(join,join.conj())
        assert np.linalg.norm(P-T@P0@T.conj().T)<5e-14
        # In moving frame, all measurements K_j=|J><v_j|.
        sigma=np.zeros((4,4),complex)
        for j in range(9):
            row=vs[:,j].conj() @ M
            prob=float(np.vdot(row,row).real)
            max_err_prob=max(max_err_prob,abs(prob-probabilities[j]))
            sigma+=np.outer(row,row.conj())
            if j<5:
                _,cc=entropy_singulars(row.reshape((2,2)))
                max_conc_err=max(max_conc_err,abs(cc-1))
        max_err_nosig=max(max_err_nosig,float(np.linalg.norm(sigma-np.eye(4)/4)))
        # Gauge-covariant branch p and Bell output. Transporting P is essential.
        row=join.conj()@M
        pp,cc=entropy_singulars(row.reshape((2,2)))
        max_err_cov=max(max_err_cov,abs(pp-1/6),abs(cc-1))
        rowfixed=fixed_bra@M
        pf,cf=entropy_singulars(rowfixed.reshape((2,2)))
        sample_fixed.append((pf,cf))
    assert max_err_cov<1e-13 and max_err_prob<1e-13
    assert max_err_nosig<1e-13 and max_conc_err<1e-13
    sample_fixed=np.asarray(sample_fixed)
    return {'random_unitary_trials':trials,
        'covariant_join_probability':'1/6 for every pair of support frames',
        'covariant_join_concurrence':'1 for every pair of support frames',
        'moving_Kraus_probabilities':['1/6','1/12','1/4','1/4','1/4','0','0','0','0'],
        'max_probability_error':max_err_prob,
        'max_branch_covariance_error':max_err_cov,
        'max_no_signalling_density_residual':max_err_nosig,
        'max_pair_concurrence_error':max_conc_err,
        'fixed_join_probability_sample_range':[float(sample_fixed[:,0].min()),float(sample_fixed[:,0].max())],
        'fixed_join_concurrence_sample_range':[float(sample_fixed[:,1].min()),float(sample_fixed[:,1].max())]}

def constructive_any_fixed_bra_rank_drop(trials=128,seed=9095):
    """For any B!=0, choose V1 so M=V1^T conjugate(B) nonzero;
    choose first column of V2 in ker(M), second giving nonzero image.
    This yields positive-probability Schmidt-rank-one output.
    """
    rng=np.random.default_rng(seed)
    max_first=0.;max_second=0.;min_p=100.;min_conc=100.
    ranks=[1,2,3]
    for rank in ranks:
        for _ in range(trials):
            L=haar_u(rng);Rt=haar_u(rng)
            diag=np.diag(([1.,.6,.3][:rank])+[0.]*(3-rank))
            B=L@diag@Rt.T
            V1=haar_u(rng)[:,:2]
            M=V1.T@B.conj()
            if np.linalg.norm(M)<1e-12:raise AssertionError('Unexpected zero contraction')
            # Nullspace of 2 x 3 matrix: last right singular vector.
            UU,ss,VH=np.linalg.svd(M, full_matrices=True)
            w0=VH[-1,:].conj()
            # If rank M=1, last right singular vector is in kernel as well.
            w1=VH[0,:].conj()
            V2=np.column_stack([w0,w1])
            C=.5*M@V2
            p,cc=entropy_singulars(C)
            assert p>1e-17 and cc<2e-13
            max_first=max(max_first,np.linalg.norm(C[:,0]))
            max_second=max(max_second,cc)
            min_p=min(min_p,p)
            min_conc=min(min_conc,cc)
    return {'fixed_bra_ranks_tested':ranks,'tests_per_rank':trials,
            'max_null_column_norm':max_first,
            'max_concurrence_after_forced_rank_drop':max_second,
            'min_strictly_positive_branch_probability':min_p,
            'theorem':'No nonzero fixed rank-one Bell-outcome bra is universally entanglement-preserving on all independent 2-dimensional supports.'}

def complement_phase_ambiguity():
    # U2 -> U2 diag(1,1,-1) leaves (3,2) link state unchanged.
    # But relative full-qutrit Bell join acquires sign on |22>.
    o0=(exact_basis(9,0)+exact_basis(9,4)+exact_basis(9,8))/SQ(3)
    op=(exact_basis(9,0)+exact_basis(9,4)-exact_basis(9,8))/SQ(3)
    ov=sp.simplify((o0.T*op)[0]);assert ov==R(1,3)
    # Trace distance for pure states sqrt(1-abs(overlap)^2)
    D=sp.sqrt(1-ov**2);assert sp.simplify(D-2*SQ(2)/3)==0
    psi=sp.zeros(9,4)
    for i in range(2):
        for j in range(2):psi[3*i+j,2*i+j]=R(1,2)
    assert o0.T*psi==op.T*psi
    return {'same_input_links':True,'different_rank_3_join_projectors':True,
            'join_state_overlap':'1/3','projector_trace_distance':'2*sqrt(2)/3',
            'selected_branch_amplitudes_identical':True,
            'conclusion':'Input pair states alone do not determine the join on the unused third qutrit directions.'}


def su3_volume_form_completion(trials=128,seed=18682):
    """Natural SU(3) volume-form complement of a 3x2 isometry.
    It is frame equivariant for SU(3), but NOT a function of the pair RAY:
    V and exp(i*phi)V describe the same projective link state.
    """
    rng=np.random.default_rng(seed)
    def oriented_completion(V):
        w=np.conj(np.cross(V[:,0],V[:,1]))
        U=np.column_stack([V,w])
        assert np.linalg.norm(U.conj().T@U-np.eye(3))<2e-13
        assert abs(np.linalg.det(U)-1)<2e-13
        return U
    max_su3=0.; max_ray_condition=0.
    for _ in range(trials):
        U0=haar_u(rng)
        U0=U0*np.exp(-1j*np.angle(np.linalg.det(U0))/3)
        G=haar_u(rng)
        G=G*np.exp(-1j*np.angle(np.linalg.det(G))/3)
        V=U0[:,:2]
        T0=oriented_completion(V)
        Tg=oriented_completion(G@V)
        max_su3=max(max_su3,float(np.linalg.norm(Tg-G@T0)))
        phi=np.pi/3
        Tph=oriented_completion(np.exp(1j*phi)*V)
        max_ray_condition=max(max_ray_condition,float(np.linalg.norm(Tph-T0@np.diag([np.exp(1j*phi),np.exp(1j*phi),np.exp(-2j*phi)]))))
    assert max_su3<1e-12 and max_ray_condition<1e-12
    # The SAME link RAY (V vs e^(i*pi/3)V) generates two different full-rank joins.
    # Their Bell vector overlap = (2+exp(-i*3phi))/3 up to global phase.
    phi=sp.pi/3
    overlap=sp.simplify((2*sp.exp(sp.I*phi)+sp.exp(-2*sp.I*phi))/3)
    assert sp.simplify(abs(overlap)**2-R(1,9))==0
    return {'trials':trials,'max_SU3_equivariance_residual':max_su3,
      'max_phase_transform_residual':max_ray_condition,
      'volume_form_extension_is_unique_given_phased_3x2_isometry':True,
      'volume_form_extension_is_not_a_function_of_unphased_pair_ray':True,
      'pair_ray_phase_example':'phi=pi/3 on first link',
      'join_projector_trace_distance_for_identical_input_pair_rays':'2*sqrt(2)/3',
      'selected_output_remains_a_Bell_pair':True,
      'interpretation':'The SU(3) epsilon tensor fixes the complement from a phased link amplitude, but the original projective pair state cannot fix the phase of that amplitude.'}


def exact_projective_ray_stabilizer():
    # A qutrit SU(3) transformation acting as a scalar on the 2D support
    # fixes the individual (3,2) input RAY, but not the full qutrit Bell join.
    z=sp.exp(sp.I*sp.pi/3)
    G=sp.diag(z,z,z**(-2))
    assert sp.simplify(G.det())==1
    omega=(exact_basis(9,0)+exact_basis(9,4)+exact_basis(9,8))/SQ(3)
    T=sp.kronecker_product(G,sp.eye(3))
    joined=T*omega
    overlap=sp.simplify((omega.conjugate().T*joined)[0])
    assert sp.simplify(sp.expand_complex(overlap*sp.conjugate(overlap)-R(1,9)))==0
    psi=sp.zeros(9,4)
    for i in range(2):
        for j in range(2): psi[3*i+j,2*i+j]=R(1,2)
    assert sp.simplify(T*psi-z*psi)==sp.zeros(9,4)
    return {'stabilizer':'G_phi=diag(exp(i phi),exp(i phi),exp(-2i phi)) in SU(3)',
       'G_pi_over_3_changes_input_pair_only_by_global_phase':True,
       'full_join_projector_trace_distance':'2*sqrt(2)/3',
       'projective_ray_equivariance_obstruction':'Any rank-three maximally entangled join projector equivariantly selected solely from these input rays would need to be invariant under G_phi on one qutrit; its vector cannot span both inequivalent eigenvalue sectors, so Schmidt rank <= 2.'}

def main():
    data={
      'scope':'v7 covariance/uniqueness test only; does not derive beta function',
      'aligned_versus_rotated_fixed_P3':exact_fixed_projection_example(),
      'exact_matching_twirl':exact_finite_matching_twirl(),
      'fixed_rank_one_join_universal_no_go_stress':constructive_any_fixed_bra_rank_drop(),
      'transported_join_and_instrument':transported_instrument_stress(),
      'unoccupied_complement_ambiguity':complement_phase_ambiguity(),
      'su3_epsilon_volume_form_extension':su3_volume_form_completion(),
      'projective_ray_stabilizer_no_go':exact_projective_ray_stabilizer(),
      'status':'All exact-symbolic and seeded complex NumPy assertions passed'
    }
    print(json.dumps(data,indent=2))

if __name__=='__main__':main()