// Copyright 2025 The University of Hong Kong, Department of Computer Science
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

package tools

func EditingDistance(s1, s2 string) (int, int, int, int) {
	m, n := len(s1), len(s2)

	dp := make([][]int, m+1)
	for i := range dp {
		dp[i] = make([]int, n+1)
	}

	for i := 0; i <= m; i++ {
		dp[i][0] = i
	}
	for j := 0; j <= n; j++ {
		dp[0][j] = j
	}

	for i := 1; i <= m; i++ {
		for j := 1; j <= n; j++ {
			if s1[i-1] == s2[j-1] {
				dp[i][j] = dp[i-1][j-1]
			} else {
				dp[i][j] = min(
					dp[i-1][j-1]+1,
					dp[i][j-1]+1,
					dp[i-1][j]+1,
				)
			}
		}
	}

	i, j := m, n
	substitutions, insertions, deletions := 0, 0, 0

	for i > 0 || j > 0 {
		if i > 0 && j > 0 && s1[i-1] == s2[j-1] {
			i--
			j--
		} else if i > 0 && j > 0 && dp[i][j] == dp[i-1][j-1]+1 {
			substitutions++
			i--
			j--
		} else if j > 0 && dp[i][j] == dp[i][j-1]+1 {
			insertions++
			j--
		} else {
			deletions++
			i--
		}
	}

	return dp[m][n], substitutions, insertions, deletions
}

func Hamming_Distrance(dna1, dna2 []rune) int {
	if len(dna1) != len(dna2) {
		return -1
	}
	count := 0
	for i := 0; i < len(dna1); i++ {
		if dna1[i] != dna2[i] {
			count++
		}
	}
	return count
}
