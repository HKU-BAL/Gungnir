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
