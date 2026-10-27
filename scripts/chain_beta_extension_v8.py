#!/usr/bin/env python3
"""Chain Theory QED beta program v8: gauge descent of qutrit join and 64 qutrit records.
Source: ChainTheoryAlpha.pdf (2026-10-23). Independent exact and seeded numerical checks.
Nothing here independently derives QED running.
"""
import json
import math
import hashlib
from pathlib import Path
import numpy as np
import sympy as s

sqrt=s.sqrt
I3=s.eye(3); I2=s.eye(2); I9=s.eye(9)
e=lambda n,i:s.eye(n)[:,i]
z=s.exp(s.I*s.pi/3)
G=s.diag(z,z,z**(-2))
assert s.simplify(s.det(G))==1

# Fixed basis of two qutrit ends and two remote qubits.
B2=(s.kronecker_product(e(3,0),e(3,0))+s.kronecker_product(e(3,1),e(3,1)))/sqrt(2)
C3=s.kronecker_product(e(3,2),e(3,2))
B3=(sqrt(2)*B2+C3)/sqrt(3)
P3=s.simplify(B3*B3.conjugate().T)
Pocc=s.diag(1,1,0)
Socc=s.kronecker_product(Pocc,Pocc)
Psi=s.zeros(9,4)
for a in range(2):
    for b in range(2):Psi[3*a+b,2*a+b]=s.Rational(1,2)


def exact_join_and_stabilizer():
    # The projective input ray transforms by a global phase z.
    input_pair=(s.kronecker_product(e(3,0),e(2,0))+s.kronecker_product(e(3,1),e(2,1)))/sqrt(2)
    assert all(s.simplify(x)==0 for x in ((s.kronecker_product(G,I2)*input_pair)-z*input_pair))
    T=s.kronecker_product(G,I3)
    Bj=s.simplify(T*B3)
    Bjs=(sqrt(2)*B2-C3)/sqrt(3)
    Pj=s.simplify(Bjs*Bjs.T)
    overlap=s.simplify((B3.conjugate().T*Bj)[0])
    # All phases in Bj: z on |00>,|11>, and -z on |22>.
    assert s.simplify(Bj-z*(sqrt(2)*B2-C3)/sqrt(3))==s.zeros(9,1)
    assert s.simplify(s.expand_complex(s.conjugate(overlap)*overlap-s.Rational(1,9)))==0
    assert s.simplify(Socc*(Pj-P3)*Socc)==s.zeros(9)
    # The full rank-three effect changes, but the occupied restriction is unchanged.
    assert Pj!=P3
    Cfull=s.simplify((B3.conjugate().T*Psi))
    Crank2=s.simplify((B2.conjugate().T*Psi))
    pf=s.simplify((Cfull*Cfull.conjugate().T)[0])
    pr=s.simplify((Crank2*Crank2.conjugate().T)[0])
    assert pf==s.Rational(1,6) and pr==s.Rational(1,4)
    bell2=(e(4,0)+e(4,3))/sqrt(2)
    assert s.simplify(Cfull/sqrt(pf)-bell2.T)==s.zeros(1,4)
    assert s.simplify(Crank2/sqrt(pr)-bell2.T)==s.zeros(1,4)
    # Pure-state projective fidelity trace distance.
    dist=s.simplify(s.sqrt(1-s.expand_complex(s.conjugate(overlap)*overlap)))
    assert s.simplify(dist-2*sqrt(2)/3)==0
    # Q_center (hypercharge, not electromagnetic charge) is scalar on (3,2) tensor link.
    Qanti=s.I*s.diag(3,3,-2,-2,-2)
    Q_on_link=s.I*(3-2)*s.eye(6)
    assert s.simplify(Q_on_link-s.I*s.eye(6))==s.zeros(6)
    return dict(projective_ray_stabilizer='G=diag(exp(i*pi/3),exp(i*pi/3),exp(-2*i*pi/3)) in SU(3)',
       identical_projective_input=True,full_join_projector_changed=True,
       full_join_trace_distance='2*sqrt(2)/3',
       occupied_block_difference='zero exactly',
       full_qutrit_join_success='1/6',support_only_join_success='1/4',
       both_selected_remote_states='|omega_2>',
       central_generator_Q_on_link='i*I_6 (not Q_em)')


def exact_ray_stabilizer_twirl():
    # Average independently admissible ray phases; terms within S={0,1} keep the same phase.
    # Analytic 3-phase twirl makes occupied/complement interference disappear.
    E=s.Rational(2,3)*(B2*B2.conjugate().T)+s.Rational(1,3)*(C3*C3.conjugate().T)
    assert s.simplify(E.trace()-1)==0
    assert s.simplify(E*E-E)!=s.zeros(9)
    assert s.simplify(Socc*(E-P3)*Socc)==s.zeros(9)
    assert s.simplify(E*E).trace()==s.Rational(5,9)
    # For original aligned input, effective selected-branch unnormalized Bob state.
    sigma=s.simplify(Psi.conjugate().T*E*Psi)
    b=(e(4,0)+e(4,3))/sqrt(2)
    assert s.simplify(sigma-s.Rational(1,6)*b*b.conjugate().T)==s.zeros(4)
    # An independent rotated-support input is NOT fixed by this mixed effect.
    rotated=s.zeros(9,4)
    Afirst=[0,2];Asecond=[0,1]
    for i in range(2):
        for j in range(2):rotated[3*Afirst[i]+Asecond[j],2*i+j]=s.Rational(1,2)
    rotated_sigma=s.simplify(rotated.T*E*rotated)
    assert rotated_sigma==s.Rational(1,12)*e(4,0)*e(4,0).T
    # Verify finite 3-phase average numerically with 0, 2*pi/9, 4*pi/9.
    P3n=np.array(P3,complex); En=np.array(E,complex)
    avg=np.zeros((9,9),complex)
    for k in range(3):
        phi=2*np.pi*k/9
        gn=np.diag([np.exp(1j*phi),np.exp(1j*phi),np.exp(-2j*phi)])
        W=np.kron(gn,np.eye(3))
        avg+=W@P3n@W.conj().T/3
    residual=float(np.linalg.norm(avg-En))
    assert residual < 3e-15
    return dict(twirl_group='three phase representatives phi=0,2*pi/9,4*pi/9',
       twirled_effect='E=(2/3)|omega2_occ><omega2_occ|+(1/3)|22><22|',
       rank=2,purity='5/9',
       same_restricted_effect_on_occupied_block=True,
       original_input_probability='1/6',
       original_input_conditional_remote_state='|omega_2>',
       rotated_support_example_probability='1/12',
       rotated_support_example_remote_state='|00> (product)',
       finite_phase_twirling_residual=residual,
       qualification='E is a mixed effect, not a rank-one qutrit Bell projector')


def gellmann():
    arr=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        v=np.zeros((3,3),complex);v[i,j]=v[j,i]=1;arr.append(v)
        v=np.zeros((3,3),complex);v[i,j]=-1j;v[j,i]=1j;arr.append(v)
    arr.append(np.diag([1,-1,0]).astype(complex))
    arr.append(np.diag([1,1,-2]).astype(complex)/np.sqrt(3))
    assert len(arr)==8
    assert np.max(np.abs(np.array([[np.trace(a@b) for b in arr] for a in arr])-2*np.eye(8)))<1e-14
    return arr


def exact_sector_swap():
    # G has adjoint eigenvalues + four times and - four times on su(3).
    # First pair (01) unchanged, while 02 and 12 flip both.
    signs=[1,1,-1,-1,-1,-1,1,1]
    D=s.diag(*signs)
    assert D.trace()==0
    S=s.zeros(64)
    for i in range(8):
        for j in range(8):S[8*j+i,8*i+j]=1
    W=s.kronecker_product(D,s.eye(8))
    Pm=(s.eye(64)-S)/2
    Pp=(s.eye(64)+S)/2
    # Source S7: X = -swap on same-end products; even/symmetric is X=-1, weight sqrt(3).
    M=s.simplify(sqrt(3)*Pp+Pm)
    Mprime=W*M*W.T
    assert (S*S)==s.eye(64)
    assert s.simplify(M.trace()-(28+36*sqrt(3)))==0
    assert s.simplify(Mprime.trace()-M.trace())==0
    sq=s.simplify(s.trace((M-Mprime)*(M-Mprime)))
    assert s.simplify(sq-32*(sqrt(3)-1)**2)==0
    # Correlated action preserves the swap; for this G, Ad(conj G) = Ad(G) = D.
    Vcorr=s.kronecker_product(D,D)
    assert Vcorr*S*Vcorr.T==S
    # Distinguish state/probe covariance from active independent-frame invariance.
    u=(e(8,0)+e(8,2))/sqrt(2)
    probe=s.kronecker_product(u,u)
    val=s.simplify((probe.T*M*probe)[0])
    valp=s.simplify((probe.T*Mprime*probe)[0])
    assert val==sqrt(3)
    assert s.simplify(valp-(1+sqrt(3))/2)==0
    assert s.simplify(((W*probe).T*Mprime*(W*probe))[0]-val)==0
    # Check actual action via Gell-Mann adjoint of original SU3 subgroup.
    Gnum=np.array(G.evalf(),complex)
    basis=gellmann()
    A=np.array([[np.trace(a@Gnum@b@Gnum.conj().T).real/2 for b in basis] for a in basis])
    err=float(np.linalg.norm(A-np.diag(signs)))
    assert err<1e-14
    return dict(qutrit_two_slot_sector_size=64,
      symmetric_swap_eigenvectors=36,antisymmetric_swap_eigenvectors=28,
      sector_trace='28+36*sqrt(3)',
      independent_ray_gauge_adjoint='diag(+1,+1,-1,-1,-1,-1,+1,+1)',
      adjoint_numeric_residual=err,
      gauge_changed_M_HS_squared='32*(sqrt(3)-1)^2',
      unchanged_trace=True,
      matched_correlated_gauge_preserves_swap=True,
      fixed_probe_before='sqrt(3)',fixed_probe_after='(1+sqrt(3))/2',
      cotransported_probe_expectation='sqrt(3)',
      qualification='64-direction sub-sector only; not a full transport construction of all 172 record operators')


def haar_su3(rng):
    x=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    q,r=np.linalg.qr(x)
    q=q@np.diag(np.conjugate(np.diag(r)/np.abs(np.diag(r))))
    return q*np.exp(-1j*np.angle(np.linalg.det(q))/3)


def random_covariant_rays(trials=128,seed=81631):
    rng=np.random.default_rng(seed)
    ket2=np.array(B2,complex).reshape(9)
    ket3=np.array(B3,complex).reshape(9)
    psin=np.array(Psi,complex)
    max_delta=0;max_phase=0
    for _ in range(trials):
        U1,U2=haar_su3(rng),haar_su3(rng)
        ph1,ph2=rng.uniform(-np.pi,np.pi,2)
        V1=U1[:,:2];V2=U2[:,:2]
        V1p=np.exp(1j*ph1)*V1
        V2p=np.exp(1j*ph2)*V2
        # A support-joined (rank-2 Schmidt) state, intrinsically phase independent as a ray.
        J2=sum(np.kron(V1[:,i],V2[:,i]) for i in range(2))/np.sqrt(2)
        J2p=sum(np.kron(V1p[:,i],V2p[:,i]) for i in range(2))/np.sqrt(2)
        max_phase=max(max_phase, float(np.linalg.norm(J2p-np.exp(1j*(ph1+ph2))*J2)))
        # The input is the product of (3,2) links, with A1A2 versus B1B2 ordering.
        psim=np.zeros((9,4),complex)
        for i in range(2):
            for j in range(2):psim[:,2*i+j]+=.5*np.kron(V1[:,i],V2[:,j])
        row=J2.conj()@psim
        norm2=float(np.vdot(row,row).real)
        bell=(np.eye(4)[0]+np.eye(4)[3])/np.sqrt(2)
        max_delta=max(max_delta, abs(norm2-.25),np.linalg.norm(row/.5-bell))
        # Corresponding transported full qutrit join also works with p=1/6.
        J3=np.kron(U1,U2)@ket3
        row3=J3.conj()@psim
        max_delta=max(max_delta, abs(np.vdot(row3,row3).real-1/6),np.linalg.norm(row3/np.sqrt(1/6)-bell))
    assert max_delta<2e-14 and max_phase<2e-14
    return dict(random_SU3_support_pairs=trials,
       support_only_join_success='1/4',transformed_full_join_success='1/6',
       maximum_probability_and_output_error=max_delta,
       maximum_pair_ray_phase_residual=max_phase,
       caveat='Relative identification of the two qubit Schmidt bases is fixed in this test.')


def main():
    data={'program':'ChainTheory beta v8','source':'ChainTheoryAlpha.pdf',
          'status':'conditional mathematical tests; no beta derivation',
          'exact_join':exact_join_and_stabilizer(),
          'ray_stabilizer_twirl':exact_ray_stabilizer_twirl(),
          '64_record_sector':exact_sector_swap(),
          'random_covariant_inputs':random_covariant_rays()}
    print(json.dumps(data,indent=2,ensure_ascii=False))

if __name__=='__main__':main()