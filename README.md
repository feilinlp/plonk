# PLONK Zero-Knowledge Proof System Implementation

A C++ implementation of the PLONK (Permutations over Lagrange-bases for Oecumenical Noninteractive arguments of Knowledge) zero-knowledge proof system, based on the original paper by Ariel Gabizon, Zachary J. Williamson, and Oana Ciobotaru.

## Overview

PLONK is a universal and updateable zk-SNARK that enables efficient zero-knowledge proofs for arbitrary computation. This implementation provides a complete proof system including circuit setup, proof generation, and verification.

## Paper Reference

This implementation is based on the original PLONK paper:

- **Title**: "PLONK: Permutations over Lagrange-bases for Oecumenical Noninteractive arguments of Knowledge"
- **Authors**: Ariel Gabizon, Zachary J. Williamson, Oana Ciobotaru
- **URL**: https://eprint.iacr.org/2019/953.pdf
- **Published**: 2019

## Features

### Implemented Components

- **Circuit Representation**: Support for arithmetic circuits with multiplication, linear, and constant gates
- **Permutation Arguments**: Copy constraints using permutation polynomials over cosets
- **KZG Polynomial Commitments**: Kate-Zaverucha-Goldberg commitments for hiding polynomials
- **Fiat-Shamir Heuristic**: Non-interactive proof generation using cryptographic hashing
- **Five-Round Protocol**: Complete implementation of all PLONK rounds
- **Verification**: Full verifier with pairing-based checks

### Core Algorithms

1. **Circuit Setup**: Transform arithmetic circuits into PLONK constraint system
2. **Trusted Setup**: Generate structured reference string (SRS) for polynomial commitments
3. **Proof Generation**: Five-round protocol generating succinct proofs
4. **Verification**: Efficient verification using bilinear pairings

## Dependencies

- **MCL Library**: For elliptic curve operations and finite field arithmetic
- **OpenSSL**: For cryptographic hash functions (SHA-256)

## Algorithm Details

### Five-Round Protocol

1. **Round 1**: Commit to witness polynomials a(X), b(X), c(X)
2. **Round 2**: Generate permutation polynomial z(X) using β, γ challenges
3. **Round 3**: Construct quotient polynomial t(X) and split into low/mid/high degree parts
4. **Round 4**: Open polynomials at evaluation point ζ
5. **Round 5**: Batch opening proofs using linearization polynomial

### Key Cryptographic Components

- **Permutation Arguments**: Enforce copy constraints using grand product arguments
- **Polynomial Commitments**: Hide polynomial evaluations while allowing verification
- **Random Linear Combinations**: Batch multiple polynomial checks efficiently
- **Fiat-Shamir**: Derive verifier challenges from proof transcript

## Current Limitations & Potential Improvements

### 1. NTT Evaluations on Lagrange Basis

**Current State**: NTT operations performed in coefficient form

```cpp
// Current: Convert between coefficient and evaluation form
ntt_inverse(e, omega); // Coefficient to evaluation
```

**Improvement Needed**: Direct NTT operations on Lagrange basis polynomials

```cpp
// Target: Native Lagrange basis operations
lagrange_ntt(lagrange_poly, omega); // Direct Lagrange NTT
```

## References

1. Gabizon, A., Williamson, Z. J., & Ciobotaru, O. (2019). PLONK: Permutations over Lagrange-bases for Oecumenical Noninteractive arguments of Knowledge. Cryptology ePrint Archive, Report 2019/953. https://eprint.iacr.org/2019/953.pdf

2. Kate, A., Zaverucha, G. M., & Goldberg, I. (2010). Constant-size commitments to polynomials and their applications. International Conference on the Theory and Application of Cryptology and Information Security. https://www.iacr.org/archive/asiacrypt2010/6477178/6477178.pdf
