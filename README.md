# Optimization-of-NB-LDPC-row-codes-using-Fossorier-Declerq
The GitHub repository contains a tool for optimizing the row non-binary labeling of Low-Density Parity-Check (LDPC) codes using the Fossorier-Declerq approach. The tool improves the weight spectrum of binary images of row codes using the companion matrix of GF(q) field, leading to an improvement in the code distance of non-binary codes according to various bounds.

The tool also includes a MATLAB implementation of the Union Bound and Weight Distribution Estimation, based on Morelos-Zaragoza's MATLAB source code available at http://the-art-of-ecc.com/.

The general idea behind the optimization of row non-binary labeling is taken from the paper "Design of regular (2,d/sub c/)-LDPC codes over GF(q) using their binary images" by C. Poulliat, M. Fossorier, and D. Declercq, published in IEEE Transactions on Communications in 2008.

This tool can also be applied to optimize Generalized LDPC (GLDPC) codes, which can be considered as a relaxed version of NB-LDPC codes. GLDPC codes allow for more freedom in representing the symbol binary image code, not only through the companion matrix but also allowing for other symbol-to-symbol mappings.

Overall, this repository provides a useful tool for optimizing the row non-binary labeling of LDPC codes and improving the weight spectrum of binary images of row codes using the companion matrix of GF(q) field.
