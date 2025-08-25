package tools

import (
	"math"
	"math/rand"
)

var (
	subRand = rand.New(rand.NewSource(1))
	insRand = rand.New(rand.NewSource(2))
	delRand = rand.New(rand.NewSource(3))
)

func ExistanceInt(r int, s []int) bool {
	for i := 0; i < len(s); i++ {
		if r == s[i] {
			return true
		}
	}
	return false
}

func deepCopyString(s string) string {
	newBytes := make([]byte, len(s))
	copy(newBytes, s)
	return string(newBytes)
}

func Substitutions(strand string) string {
	runes := []rune(strand)
	pos := subRand.Intn(len(strand))
	base_Int := Nuc2Int(runes[pos])
	base := (subRand.Intn(3) + base_Int + 1) & 3
	runes[pos] = Int2Nuc(base)
	return string(runes)
}

func Insertions(strand string) string {
	pos := insRand.Intn(len(strand) + 1)
	base := insRand.Intn(4)
	runes := []rune(strand)
	if pos == len(strand) {
		runes = append(runes, Int2Nuc(base))
	} else {
		runes = append(runes[:pos], append([]rune{Int2Nuc(base)}, runes[pos:]...)...)
	}

	return string(runes)
}

func Deletions(strand string) string {
	pos := delRand.Intn(len(strand))
	runes := []rune(strand)
	runes = append(runes[:pos], runes[pos+1:]...)
	return string(runes)
}

func SplitThree(n int) ([]int, []int, []int) {
	splitrand := rand.New(rand.NewSource(9))

	nums := make([]int, n)
	for i := range nums {
		nums[i] = i
	}
	splitrand.Shuffle(n, func(i, j int) {
		nums[i], nums[j] = nums[j], nums[i]
	})

	base := n / 3
	remainder := n % 3
	sizes := [3]int{base, base, base}
	r := splitrand.Intn(3)

	if remainder > 0 {
		if remainder == 1 {
			sizes[r] += 1
		} else {
			sizes[0] += 1
			sizes[1] += 1
			sizes[2] += 1
			sizes[r] -= 1
		}
	}

	offset := 0
	set1 := nums[offset : offset+sizes[0]]
	offset += sizes[0]
	set2 := nums[offset : offset+sizes[1]]
	offset += sizes[1]
	set3 := nums[offset : offset+sizes[2]]

	return set1, set2, set3
}

func AddErrorforseqs(s []string, errornum int) []string {
	res := make([]string, len(s))
	base_err := errornum / 3
	reminder := errornum % 3
	if reminder == 0 {
		for i := 0; i < len(s); i++ {
			res[i] = AddErrorWithFixNum(s[i], base_err, base_err, base_err)
		}
	} else if reminder == 1 {
		set1, set2, set3 := SplitThree(len(s))
		for i := 0; i < len(set1); i++ {
			res[set1[i]] = AddErrorWithFixNum(s[set1[i]], base_err+1, base_err, base_err)
		}
		for i := 0; i < len(set2); i++ {
			res[set2[i]] = AddErrorWithFixNum(s[set2[i]], base_err, base_err+1, base_err)
		}
		for i := 0; i < len(set3); i++ {
			res[set3[i]] = AddErrorWithFixNum(s[set3[i]], base_err, base_err, base_err+1)
		}
	} else if reminder == 2 {
		set1, set2, set3 := SplitThree(len(s))
		for i := 0; i < len(set1); i++ {
			res[set1[i]] = AddErrorWithFixNum(s[set1[i]], base_err, base_err+1, base_err+1)
		}
		for i := 0; i < len(set2); i++ {
			res[set2[i]] = AddErrorWithFixNum(s[set2[i]], base_err+1, base_err, base_err+1)
		}
		for i := 0; i < len(set3); i++ {
			res[set3[i]] = AddErrorWithFixNum(s[set3[i]], base_err+1, base_err+1, base_err)
		}
	}
	return res
}

func AddErrorWithFixNum(s string, sub, del, ins int) string {
	res := deepCopyString(s)
	for i := 0; i < ins; i++ {
		res = Insertions(res)
	}
	for i := 0; i < sub; i++ {
		res = Substitutions(res)
	}
	for i := 0; i < del; i++ {
		res = Deletions(res)
	}
	return res
}

func DistributeErrs(n, sub, ins, del int) ([]int, []int, []int) {
	minerr := (sub + ins + del) / n

	subpos := make([]int, 0)
	inspos := make([]int, 0)
	delpos := make([]int, 0)

	Newrand := rand.New(rand.NewSource(7))

	perm := Newrand.Perm(n)
	count := 0

	for i := 0; i < n; i++ {
		if count < ins {
			inspos = append(inspos, perm[i])
			count++
			continue
		} else if count < (ins + del) {
			delpos = append(delpos, perm[i])
			count++
			continue
		}
	}

	for i := 0; i < n; i++ {
		if count < (ins + del) {
			delpos = append(delpos, perm[i])
			count++
			continue
		}
	}

	perm = Newrand.Perm(n)
	count = 0

	for i := 0; i < n; i++ {
		if count < sub {
			errnum := 0
			if ExistanceInt(perm[i], inspos) {
				errnum++
			}
			if ExistanceInt(perm[i], delpos) {
				errnum++
			}
			if errnum < minerr {
				subpos = append(subpos, perm[i])
				count++
			}
		}
	}

	for i := 0; i < n; i++ {
		if count < sub {
			if !ExistanceInt(perm[i], subpos) {
				subpos = append(subpos, perm[i])
				count++
			}
		}
	}

	return subpos, inspos, delpos
}

func AddError(s []string, subrate, delrate, insrate, errrate float64) []string {
	res := make([]string, len(s))
	for i := 0; i < len(s); i++ {
		res[i] = deepCopyString(s[i])
	}
	nuc_num := len(s) * len(s[0])
	Ins := int(math.Ceil(float64(nuc_num) * insrate))
	Del := int(math.Ceil(float64(nuc_num) * delrate))
	Sub := int(math.Ceil(float64(nuc_num)*errrate)) - Ins - Del
	I_average := Ins / len(s)
	I_reminder := Ins - len(s)*I_average
	D_average := Del / len(s)
	D_reminder := Del - len(s)*D_average
	S_average := Sub / len(s)
	S_reminder := Sub - len(s)*S_average

	subpos, inspos, delpos := DistributeErrs(len(s), S_reminder, I_reminder, D_reminder)

	for i := 0; i < len(s); i++ {
		ins := I_average
		if ExistanceInt(i, inspos) {
			ins += 1
		}
		del := D_average
		if ExistanceInt(i, delpos) {
			del += 1
		}
		sub := S_average
		if ExistanceInt(i, subpos) {
			sub += 1
		}
		res[i] = AddErrorWithFixNum(s[i], sub, del, ins)
	}

	return res
}
