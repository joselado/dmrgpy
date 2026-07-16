# Golden values computed with the pre-optimization MultiOperator
# implementation (git commit bbf88b7, the parent of the "Remove
# accidental O(n^2) blowup in MultiOperator algebra" commit), i.e.
# before multioperator.py's copy()/clean()/get_dagger()/jordan_wigner()
# were rewritten for speed. The test suite in this directory checks
# that the current (optimized) implementation still reproduces these
# numbers exactly (or, for iterative DMRG results, to within DMRG
# convergence tolerance).
#
# Regenerate by running gen_reference.py (kept alongside this file's
# history in the PR description) with "--orig" against a checkout of
# the pre-optimization commit, and without against the current tree;
# the two outputs should match before updating this file.

HEISENBERG_DMRG_ENERGY = -4.25803520681935

HEISENBERG_ED_ENERGY = -4.258035207282873

SET_EXCHANGE_NTERMS = 30  # after an explicit .clean()

SET_EXCHANGE_DMRG_ENERGY = -4.98715426777585

SET_EXCHANGE_ED_ENERGY = -4.987154267775847

FERMION_CORRELATOR_U0 = [0.4999999999999999,
 -0.4355596199317578,
 -6.938893903907228e-18,
 0.19384226684174413,
 -6.938893903907228e-18,
 -0.15070830458391019]

FERMION_CORRELATOR_U2 = [0.4999999999999979,
 -0.39957120881017805,
 7.945033519973776e-16,
 0.15202937529797725,
 -4.510281037539698e-16,
 -0.09909989440229558]

FERMION_CORRELATOR_U10 = [0.50000000000004,
 -0.17058307465618394,
 -3.52062129449493e-15,
 0.01205263580785773,
 2.409097227262791e-16,
 -0.001489343178748384]

# Jordan-Wigner transform of Cdag_i * C_j (spinless fermions, to_terms()
# output), keyed by "i_j". Each value is a tuple of (real, imag, ops)
# terms, sorted for order-independent comparison.
JW_TERMS_CDAGC = {'0_0': ((1.0, 0.0, (('Adag', 1), ('A', 1))),),
 '0_1': ((1.0, 0.0, (('Adag', 1), ('A', 2))),),
 '0_3': ((1.0, 0.0, (('Adag', 1), ('F', 2), ('F', 3), ('A', 4))),),
 '1_0': ((1.0, 0.0, (('A', 1), ('Adag', 2))),),
 '2_5': ((1.0, 0.0, (('Adag', 3), ('F', 4), ('F', 5), ('A', 6))),),
 '3_0': ((1.0, 0.0, (('A', 1), ('F', 2), ('F', 3), ('Adag', 4))),)}

# Jordan-Wigner transform of Cdag_0 * C_2 * Cdag_3 * C_1
JW_TERMS_4FERMION = ((1.0, 0.0, (('Adag', 1), ('F', 2), ('A', 3), ('A', 2), ('F', 3), ('Adag', 4))),)

# to_terms() of (Cdag_1 * C_4).get_dagger()
DAGGER_TERMS_CDAG1_C4 = ((1.0, 0.0, (('A', 2), ('F', 3), ('F', 4), ('Adag', 5))),)

FERMION_HERMITIAN_H_IS_HERMITIAN = True

SPIN_COMMUTATOR_IS_ZERO = True

SPIN_PICTURE_ENERGY = -1.616025403784438

SPINFUL_EXCHANGE_ENERGY = -2.2360679774997867
