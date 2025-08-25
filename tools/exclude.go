package tools

import "fmt"

const HomopolymerLen = 3

var Exclude_motifs = []string{
	// *******************
	// "GTGCAC",
	// "TCTAGA",
	// "CTCGAG",
	// "GAAGAC",
	// "GGTCTC",
	// "GTCTTC",
	// "GAGACC",
	// "GATATC",
	// "GAATTC",
	// "AAGCTT",
	// "CCATGG",
	// "GCGGGC",
	// "CTGCAG",
	// "CAGCTG",
	// "GTCGAC",
	// "CCCGGG",
	// "TTTAAA",
	// *******************
}

var Exclude_long_motifs = []string{
	// *******************
	// "GCCNNNNNGGC",
	// "CCTNNNNNAGG",
	// *******************
}

func Existance(r rune, s []rune) bool {
	for i := 0; i < len(s); i++ {
		if r == s[i] {
			return true
		}
	}
	return false
}

func CheckMotifs(s []rune, params Params) []rune {
	res := make([]rune, 0)
	if len(s) != params.PreviousNuc {
		fmt.Println("Invalid PreviousNuc Length!")
		return nil
	}
	lasts := s[params.PreviousNuc-1]
	homosymbol := true
	for i := 0; i < HomopolymerLen-1; i++ {
		if s[params.PreviousNuc-2-i] != lasts {
			homosymbol = false
			break
		}
	}
	if homosymbol {
		res = append(res, lasts)
	}
	for i := 0; i < len(Exclude_motifs); i++ {
		fixed_motif_pre := len(Exclude_motifs[i]) - 1
		if Exclude_motifs[i][:fixed_motif_pre] == string(s[params.PreviousNuc-fixed_motif_pre:]) {
			temp := []rune(Exclude_motifs[i])[fixed_motif_pre]
			if !Existance(temp, res) {
				res = append(res, temp)
			}
		}
	}

	for i := 0; i < len(Exclude_long_motifs); i++ {
		long_motif_pre := len(Exclude_long_motifs[i]) - 1
		pattern := []rune(Exclude_long_motifs[i])
		match := true
		for j := 0; j < long_motif_pre; j++ {
			if pattern[i] == 'N' {
				continue
			}
			if pattern[i] != s[params.PreviousNuc-long_motif_pre+i] {
				match = false
				break
			}
		}
		if match {
			temp := pattern[long_motif_pre]
			if !Existance(temp, res) {
				res = append(res, temp)
			}
		}
	}

	return res
}
