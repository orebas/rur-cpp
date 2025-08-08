// This file contains the fix for the F4 NORMAL macro bug
// The bug is in the delayed reduction optimization for primes <= 2^31-1

// Original buggy code in f4_reduce_import and f4_reduce_vector:
// NORMAL(y,p2);  // BUG: Using p² instead of p!

// The fix is to use NORMAL(y,p) instead of NORMAL(y,p2)

// To apply this fix to axf4_lib.c, search for all occurrences of:
// NORMAL(y,p2)
// and replace with:
// NORMAL(y,p)

// There are two locations that need fixing:

// 1. In f4_reduce_import (around line 1710):
/*
    y = vec[z] + x - p2;
    NORMAL(y,p2);  // CHANGE TO: NORMAL(y,p);
*/

// 2. In f4_reduce_vector (around line 1750):  
/*
    y = vec[z] - x;
    NORMAL(y,p2);  // CHANGE TO: NORMAL(y,p);
*/

// Note: The variable p2 = p*p is used for delayed reduction optimization,
// but the NORMAL macro should still use p, not p2, to properly handle
// negative values in the range [-p, 0).