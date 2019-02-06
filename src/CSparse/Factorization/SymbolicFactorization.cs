﻿// -----------------------------------------------------------------------
// <copyright file="SymbolicFactorization.cs">
// Copyright (c) 2006-2016, Timothy A. Davis
// Copyright (c) 2012-2016, Christian Woltering
// </copyright>
// -----------------------------------------------------------------------

namespace CSparse.Factorization
{
    /// <summary>
    /// Symbolic Cholesky, LU, or QR factorization storage.
    /// </summary>
    public class SymbolicFactorization
    {
        public int[] pinv;     // inverse row perm. for QR, fill red. perm for Chol
        public int[] q;        // fill-reducing column permutation for LU and QR
        public int[] parent;   // elimination tree for Cholesky and QR
        public int[] cp;       // column pointers for Cholesky, row counts for QR
        public int[] leftmost; // leftmost[i] = min(find(A(i,:))), for QR
        public int m2;         // # of rows for QR, after adding fictitious rows
        public int lnz;    // # entries in L for LU or Cholesky; in V for QR
        public int unz;    // # entries in U for LU; in R for QR
    }
}
