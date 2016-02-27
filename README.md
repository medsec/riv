# RIV

Robust Initialization Vector, a flexible and performant framework for robust
authenticated encryption. 

This repository contains two versions:
- RIV with XOR-CTR using AES-128 and 128-bit CLHASH, and 128-bit tags.
- RIV with CTRT using Deoxys-BC-128-128 and 256-bit CLHASH, and 256-bit tags.

The second is an experimental version. More information regarding RIV can be
found soon in the corresponding paper.

## Usage
Each version contains:
- A reference version that should run and produce consistent results on
  32-/64-bit as well as on little- or big-endian architectures.
- An optimized version which employs AES new instructions (including aesenc for
  Counter-mode AES and pclmulqdq for carry-less multiplications in CLHASH) for
  benchmarking.
- A utils directory with tests and benchmarks.
- A Makefile for compilation.

	make 
	./all
	./ref-test
	./ni-test
	./ni-bench
	./clean

## References
### CLHASH
RIV uses a Toeplitz-extended version of the fast universal hashing family
CLHASH:
Daniel Lemire and Owen Kaser, Faster 64-bit universal hashing using
carry-less multiplications, Journal of Cryptographic Engineering (to appear)
http://arxiv.org/abs/1503.03465 Repo:
https://github.com/lemire/StronglyUniversalStringHashing

### CTRT
RIV-256 uses a version of the Counter-in-Tweak (CTRT) mode with the (AES-round-
based) Deoxys-BC-128-128 tweakable block cipher from the TWEAKEY framework. Our
version of CTRT differs in the fact that we employed XOR instead of modular
addition to increment the tweak:

	T_i = T xor <ctr>

We reimplemented the Deoxys-BC-128-128 block cipher according to the optimized
first-round CAESAR implementation of Deoxys that was contained in the SUPERCOP
framework:
http://hyperelliptic.org/ebats/supercop-20141124.tar.bz2

Thomas Peyrin and Yannick Seurin: Counter-in-Tweak: Authenticated Encryption
Modes for Tweakable Block Ciphers. https://eprint.iacr.org/2015/1049

### Deoxys-BC
Jérémy Jean and Ivica Nikolić and Thomas Peyrin:  Tweaks and Keys for Block
Ciphers: the TWEAKEY Framework 
http://eprint.iacr.org/2014/831
http://www1.spms.ntu.edu.sg/~syllab/m/index.php/Deoxys

## Disclaimer/License
Use on your own risk. Our code may be susceptible to side-channel attacks. If
you reuse our code, a link to our repo would be nice. 

The original CLHASH by Lemire and Kaser is under Apache License 2.0.

