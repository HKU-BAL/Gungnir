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

import (
	"fmt"
	"log"
	"sort"
)

type Block struct {
	BlockID int
	Hash    [][]int
	Payload [][]int
}

func BitsintoBlocks(bits []int, params Params) []Block {
	BlockNum := (len(bits) + params.PayloadLen - 1) / params.PayloadLen
	newbits := make([]int, params.PayloadLen*BlockNum)
	for i := 0; i < len(newbits); i++ {
		if i < len(bits) {
			newbits[i] = bits[i]
		} else {
			newbits[i] = 0
		}
	}

	Data := make([]Block, BlockNum)
	for i := 0; i < BlockNum; i++ {
		bits_temp := newbits[i*params.PayloadLen : (i+1)*params.PayloadLen]
		Data[i].BlockID = i
		Data[i].Payload = make([][]int, len(params.PayloadLenSet))
		sum := 0
		for j := 0; j < len(params.PayloadLenSet); j++ {
			Data[i].Payload[j] = bits_temp[sum : sum+params.PayloadLenSet[j]]
			sum += params.PayloadLenSet[j]
		}

		Data[i].Hash = make([][]int, len(params.HashLenSet))
		Data[i].Hash[0] = Hash0Payload(Data[i].BlockID, bits_temp, params)

		for j := 1; j < len(params.HashLenSet); j++ {
			Data[i].Hash[j] = HashIndex(Data[i].BlockID, Data[i].Payload[j], Data[i].Hash[0], j, params)
		}

	}
	return Data
}

func BlocksintoBits(Data []Block, params Params) (bits []int) {
	bits = make([]int, len(Data)*params.PayloadLen)
	for i := 0; i < len(Data); i++ {
		temp_payload := ReshapePayload(Data[i].Payload, params)
		for j := 0; j < params.PayloadLen; j++ {
			bits[i*params.PayloadLen+j] = temp_payload[j]
		}
	}
	return
}

func InitPattern(p int) string {
	temp := []rune(Primer)
	res := temp[len(temp)-p:]
	return string(res)
}

func Extend_rule(kmerPre []rune, previous int, blockID int, index int, data map[string]Kmer, forbid rune, params Params) []rune {

	Patterns := make([]string, 0)
	for i := 0; i < 4; i++ {
		if Int2Nuc(i) != forbid {
			temp := append(kmerPre, Int2Nuc(i))
			Patterns = append(Patterns, string(temp))
		}
	}

	sort.Slice(Patterns, func(i, j int) bool {
		return data[Patterns[i]].E_rate < data[Patterns[j]].E_rate
	})

	res := make([]rune, 2)
	res[0] = []rune(Patterns[0])[KmerSize-1]
	res[1] = []rune(Patterns[1])[KmerSize-1]

	return res
}

func Main_rule(kmerPre []rune, previous int, blockID int, index int, data map[string]Kmer) []rune {
	res := make([]rune, 2)
	scale := HashScale(blockID, index, previous)

	sA := string(append(kmerPre, 'A'))
	sC := string(append(kmerPre, 'C'))
	sG := string(append(kmerPre, 'G'))
	sT := string(append(kmerPre, 'T'))

	if scale == 0 {
		if data[sA].E_rate < data[sC].E_rate {
			res[0] = 'A'
		} else {
			res[0] = 'C'
		}
		if data[sG].E_rate < data[sT].E_rate {
			res[1] = 'G'
		} else {
			res[1] = 'T'
		}
	} else if scale == 1 {
		if data[sA].E_rate < data[sG].E_rate {
			res[0] = 'A'
		} else {
			res[0] = 'G'
		}
		if data[sC].E_rate < data[sT].E_rate {
			res[1] = 'C'
		} else {
			res[1] = 'T'
		}

	} else if scale == 2 {
		if data[sA].E_rate < data[sT].E_rate {
			res[0] = 'A'
		} else {
			res[0] = 'T'
		}
		if data[sC].E_rate < data[sG].E_rate {
			res[1] = 'C'
		} else {
			res[1] = 'G'
		}
	} else if scale == 3 {
		if data[sA].E_rate < data[sT].E_rate {
			res[1] = 'A'
		} else {
			res[1] = 'T'
		}
		if data[sC].E_rate < data[sG].E_rate {
			res[0] = 'C'
		} else {
			res[0] = 'G'
		}
	} else if scale == 4 {
		if data[sA].E_rate < data[sG].E_rate {
			res[1] = 'A'
		} else {
			res[1] = 'G'
		}
		if data[sC].E_rate < data[sT].E_rate {
			res[0] = 'C'
		} else {
			res[0] = 'T'
		}
	} else if scale == 5 {
		if data[sA].E_rate < data[sC].E_rate {
			res[1] = 'A'
		} else {
			res[1] = 'C'
		}
		if data[sG].E_rate < data[sT].E_rate {
			res[0] = 'G'
		} else {
			res[0] = 'T'
		}
	}

	return res
}

func Hash_rule(previous int, blockID int, index int) []rune {
	res := make([]rune, 2)
	val := uint64((blockID << 48) + (index << 40) + previous)
	val = Ran_hash(val) & 3
	temp := int(val)
	res[0] = Int2Nuc(temp)
	res[1] = Int2Nuc((temp + 1) & 3)
	return res
}

func Balance_Rule(gc, at int, kmerPre []rune, previous int, blockID int, index int, data map[string]Kmer, params Params) []rune {
	res := make([]rune, 2)
	symbol := JudgeGC(gc, at, params)
	if symbol == 1 {
		res[0] = 'C'
		res[1] = 'G'
	} else {
		res[0] = 'T'
		res[1] = 'A'
	}

	kmer0 := string(append(kmerPre, res[0]))
	kmer1 := string(append(kmerPre, res[1]))
	if data[kmer1].E_rate < data[kmer0].E_rate {
		temp := res[0]
		res[0] = res[1]
		res[1] = temp
	}

	return UpdateBalance(res, kmerPre, previous, blockID, index, data, params)
}

func Extend_rule_Hash(previous int, blockID int, index int, forbid rune) []rune {
	res := make([]rune, 0)
	val := uint64((blockID << 48) + (index << 40) + previous)
	temp := int(Ran_hash(val) & 3)

	for i := 0; i < 4; i++ {
		c := Int2Nuc((i + temp) & 3)
		if c != forbid {
			res = append(res, c)
		}
	}

	res = res[:2]

	return res
}

func Three_rule(previous []int, blockID int, index int, forbid rune) []rune {
	res := make([]rune, 0)
	val := make([]uint64, len(previous))
	// fmt.Println(index, previous)
	val[0] = uint64((blockID << 48) + (index << 40) + previous[0])
	for i := 1; i < len(previous); i++ {
		val[i] = uint64(previous[i])
	}
	temp := int(Ran_hash_multi(val) & 3)

	for i := 0; i < 4; i++ {
		c := Int2Nuc((i + temp) & 3)
		if c != forbid {
			res = append(res, c)
		}
	}

	res = res[:3]

	return res
}

func GenNextC(gc, at int, pattern []rune, previous int, blockID int, index int, data map[string]Kmer, params Params) []rune {
	if params.Option == Gungnir_ONT_Params {
		return GenNextC_ONT(gc, at, pattern, previous, blockID, index, data, params)
	}
	return GenNextC_Default(gc, at, pattern, previous, blockID, index, data, params)
}

func GenNextC_Default(gc, at int, pattern []rune, previous int, blockID int, index int, data map[string]Kmer, params Params) []rune {
	Exclude := CheckMotifs(pattern, params)

	if len(Exclude) == 1 {
		return Extend_rule_Hash(previous, blockID, index, Exclude[0])
	} else if len(Exclude) == 2 {
		res := make([]rune, 0)
		for i := 0; i < 4; i++ {
			c := Int2Nuc(i)
			if !Existance(c, Exclude) {
				res = append(res, c)
			}
		}
		return res
	} else if len(Exclude) == 3 {
		for i := 0; i < 4; i++ {
			c := Int2Nuc(i)
			if !Existance(c, Exclude) {
				return []rune{c, c}
			}
		}
	}

	symbol := JudgeGC(gc, at, params)
	if symbol == 1 {
		return []rune{'G', 'C'}
	} else if symbol == 0 {
		return []rune{'A', 'T'}
	}

	return Hash_rule(previous, blockID, index)
}

// previous0 no more than 40 bits
func GenNextC_Three(gc, at int, pattern []rune, previous []int, blockID int, index int, params Params) []rune {
	Exclude := CheckMotifs(pattern, params)

	if len(Exclude) == 1 {
		return Three_rule(previous, blockID, index, Exclude[0])
	} else if len(Exclude) == 2 {
		res := make([]rune, 0)
		for i := 0; i < 4; i++ {
			c := Int2Nuc(i)
			if !Existance(c, Exclude) {
				res = append(res, c)
			}
		}
		res = append(res, res[0])
		return res
	} else if len(Exclude) == 3 {
		for i := 0; i < 4; i++ {
			c := Int2Nuc(i)
			if !Existance(c, Exclude) {
				return []rune{c, c, c}
			}
		}
	}

	symbol := JudgeGC(gc, at, params)
	if symbol == 1 {
		val := make([]uint64, len(previous))
		val[0] = uint64((blockID << 48) + (index << 40) + previous[0])
		for i := 0; i < len(previous); i++ {
			val[i] = uint64(previous[i])
		}
		hashbit := Ran_hash_multi(val) & 1
		if hashbit == 0 {
			return Three_rule(previous, blockID, index, 'A')
		} else {
			return Three_rule(previous, blockID, index, 'T')
		}
	} else if symbol == 0 {
		val := make([]uint64, len(previous))
		val[0] = uint64((blockID << 48) + (index << 40) + previous[0])
		for i := 0; i < len(previous); i++ {
			val[i] = uint64(previous[i])
		}
		hashbit := Ran_hash_multi(val) & 1
		if hashbit == 0 {
			return Three_rule(previous, blockID, index, 'C')
		} else {
			return Three_rule(previous, blockID, index, 'G')
		}
	}

	return Three_rule(previous, blockID, index, 'N')
}

func GenNextC_ONT(gc, at int, pattern []rune, previous int, blockID int, index int, data map[string]Kmer, params Params) []rune {
	kmerPre := pattern[params.PreviousNuc-KmerSize+1:]
	Exclude := CheckMotifs(pattern, params)

	if len(Exclude) == 1 {
		if Exclude[0] == 'G' || Exclude[0] == 'C' {
			if JudgeGC(gc, at, params) == 0 {
				temp := Balance_Rule(gc, at, kmerPre, previous, blockID, index, data, params)
				if temp != nil {
					return temp
				}
			}
		} else {
			if JudgeGC(gc, at, params) == 1 {
				temp := Balance_Rule(gc, at, kmerPre, previous, blockID, index, data, params)
				if temp != nil {
					return temp
				}
			}
		}
		return Extend_rule(kmerPre, previous, blockID, index, data, Exclude[0], params)
	} else if len(Exclude) == 2 {
		res := make([]rune, 0)
		for i := 0; i < 4; i++ {
			c := Int2Nuc(i)
			if !Existance(c, Exclude) {
				res = append(res, c)
			}
		}
		kmer0 := string(append(kmerPre, res[0]))
		kmer1 := string(append(kmerPre, res[1]))
		if data[kmer1].E_rate < data[kmer0].E_rate {
			temp := res[0]
			res[0] = res[1]
			res[1] = temp
		}
		return res
	} else if len(Exclude) == 3 {
		for i := 0; i < 4; i++ {
			c := Int2Nuc(i)
			if !Existance(c, Exclude) {
				return []rune{c, c}
			}
		}
	}

	if JudgeGC(gc, at, params) != -1 {
		temp := Balance_Rule(gc, at, kmerPre, previous, blockID, index, data, params)
		if temp != nil {
			return temp
		}
	}

	return Main_rule(kmerPre, previous, blockID, index, data)
}

// 0: too much gc - add at; 1: too much at - add gc; -1: no need for balance
func JudgeGC(gc, at int, params Params) int {
	gc_content := float64(gc) / float64(at+gc)
	if gc_content > params.GCmax {
		return 0
	} else if gc_content < params.GCmin {
		return 1
	}
	return -1
}

func UpdateBalance(balanceC []rune, kmerPre []rune, previous int, blockID int, index int, data map[string]Kmer, params Params) []rune {
	sA := string(append(kmerPre, 'A'))
	sC := string(append(kmerPre, 'C'))
	sG := string(append(kmerPre, 'G'))
	sT := string(append(kmerPre, 'T'))
	sum := data[sA].E_rate + data[sC].E_rate + data[sG].E_rate + data[sT].E_rate
	stemp0 := string(append(kmerPre, balanceC[0]))
	stemp1 := string(append(kmerPre, balanceC[1]))
	temp_rate0 := data[stemp0].E_rate / sum
	temp_rate1 := data[stemp1].E_rate / sum
	if temp_rate0 < params.Balance_Bound && temp_rate1 < params.Balance_Bound {
		return balanceC
	} else {
		return nil
	}
}

func Block2DNA(input Block, data map[string]Kmer, params Params) string {
	Pattern := []rune(InitPattern(params.PreviousNuc))
	runes := make([]rune, params.MaxDepth)
	previous := 0
	gc := 0
	at := 0
	index := 0
	for i := 0; i < params.HashLenSet[0]; i++ {
		var bit int
		if len(input.Hash) == 0 {
			bit = 0
		} else {
			bit = input.Hash[0][i]
		}
		thisC := GenNextC(gc, at, Pattern, previous, input.BlockID, index, data, params)[bit]
		if thisC == 'G' || thisC == 'C' {
			gc++
		} else {
			at++
		}
		runes[index] = thisC
		Pattern = append(Pattern[1:], thisC)
		previous = ((previous << 1) + bit) & Uint40Mask
		index += 1
	}

	for k := 1; k < len(params.HashLenSet); k++ {
		for i := 0; i < params.HashLenSet[k]; i++ {
			var bit int
			if len(input.Hash) == 0 {
				bit = 0
			} else {
				bit = input.Hash[k][i]
			}
			thisC := GenNextC(gc, at, Pattern, previous, input.BlockID, index, data, params)[bit]
			if thisC == 'G' || thisC == 'C' {
				gc++
			} else {
				at++
			}
			runes[index] = thisC
			Pattern = append(Pattern[1:], thisC)
			previous = ((previous << 1) + bit) & Uint40Mask
			index += 1
		}
		for i := 0; i < params.PayloadLenSet[k]; i++ {
			var bit int
			if len(input.Payload) == 0 {
				bit = 0
			} else {
				bit = input.Payload[k][i]
			}
			thisC := GenNextC(gc, at, Pattern, previous, input.BlockID, index, data, params)[bit]
			if thisC == 'G' || thisC == 'C' {
				gc++
			} else {
				at++
			}
			runes[index] = thisC
			Pattern = append(Pattern[1:], thisC)
			previous = ((previous << 1) + bit) & Uint40Mask
			index += 1
		}
	}

	for i := 0; i < params.PayloadLenSet[0]; i++ {
		var bit int
		if len(input.Payload) == 0 {
			bit = 0
		} else {
			bit = input.Payload[0][i]
		}
		thisC := GenNextC(gc, at, Pattern, previous, input.BlockID, index, data, params)[bit]
		if thisC == 'G' || thisC == 'C' {
			gc++
		} else {
			at++
		}
		runes[index] = thisC
		Pattern = append(Pattern[1:], thisC)
		previous = ((previous << 1) + bit) & Uint40Mask
		index += 1
	}

	return string(runes)
}

func Gen3Slice(bitstream []int) []int {
	if len(bitstream) > 11 {
		log.Fatal("Invalid Bitstream Length!")
	}
	sum := 0
	for i := 0; i < len(bitstream); i++ {
		sum += (bitstream[i] << i)
	}
	res := make([]int, 7)
	for i := 0; i < 7; i++ {
		res[i] = sum % 3
		sum = sum / 3
	}

	return res
}

func RecoverBitsFrom3Slice(threestream []int) ([]int, bool) {
	res := make([]int, 11)
	if len(threestream) > 7 {
		log.Fatal("Invalid Threestream Length!")
	}
	sum := 0
	temp := 1
	for i := 0; i < len(threestream); i++ {
		sum += threestream[i] * temp
		temp *= 3
	}
	for i := 0; i < 11; i++ {
		res[i] = sum & 1
		sum = sum >> 1
	}
	if sum > 0 {
		return res, false
	}

	return res, true
}

func BlockBits(input Block, params Params) []int {
	res := make([]int, 0)
	for i := 0; i < params.HashLenSet[0]; i++ {
		var bit int
		if len(input.Hash) == 0 {
			bit = 0
		} else {
			bit = input.Hash[0][i]
		}
		res = append(res, bit)
	}

	for k := 1; k < len(params.HashLenSet); k++ {
		for i := 0; i < params.HashLenSet[k]; i++ {
			var bit int
			if len(input.Hash) == 0 {
				bit = 0
			} else {
				bit = input.Hash[k][i]
			}
			res = append(res, bit)
		}
		for i := 0; i < params.PayloadLenSet[k]; i++ {
			var bit int
			if len(input.Payload) == 0 {
				bit = 0
			} else {
				bit = input.Payload[k][i]
			}
			res = append(res, bit)
		}
	}

	for i := 0; i < params.PayloadLenSet[0]; i++ {
		var bit int
		if len(input.Payload) == 0 {
			bit = 0
		} else {
			bit = input.Payload[0][i]
		}
		res = append(res, bit)
	}

	return res
}

func BitsToInt(bits []int) int {
	res := 0
	for i := 0; i < len(bits); i++ {
		res += (bits[i] << i)
	}
	return res
}

// 55 + 22 + this package (11 bits to 14 bits): 88 bits (represented as 91 bits)
// previous 0: bits in temp_previous -> 14 bits + 22 previous bits;
func GenPreviousSet(bitstream []int, temp_previous int, i int) []int {
	if i < 3 {
		previous := make([]int, 1)
		previous[0] = BitsToInt(bitstream[:11*i])
		previous[0] = previous[0] << 14
		previous[0] += temp_previous
		return previous
	} else if i < 8 {
		previous := make([]int, 2)
		previous[0] = BitsToInt(bitstream[11*(i-2) : 11*i])
		previous[0] = previous[0] << 14
		previous[0] += temp_previous
		previous[1] = BitsToInt(bitstream[0 : 11*(i-2)])
		return previous
	} else {
		previous := make([]int, 2)
		previous[0] = BitsToInt(bitstream[11*(i-2) : 11*i])
		previous[0] = previous[0] << 14
		previous[0] += temp_previous
		previous[1] = BitsToInt(bitstream[11*(i-7) : 11*(i-2)])
		return previous
	}
}

func Block2DNA_Three(input Block, params Params) string {
	bitstream := BlockBits(input, params)
	Pattern := []rune(InitPattern(params.PreviousNuc))
	runes := make([]rune, params.MaxDepth)
	gc := 0
	at := 0
	index := 0
	setlen := params.Seven3num
	if params.Res3num > 0 {
		setlen += 1
	}
	threeset := make([][]int, setlen)
	for i := 0; i < params.Seven3num; i++ {
		threeset[i] = make([]int, 7)
	}
	if params.Res3num > 0 {
		threeset[params.Seven3num] = make([]int, params.Res3num)
	}

	for i := 0; i < setlen; i++ {
		maxbound := i*11 + 11
		if maxbound > params.MaxBits {
			maxbound = params.MaxBits
		}
		// fmt.Println(bitstream[i*11:maxbound], params)
		// fmt.Println(Gen3Slice(bitstream[i*11:maxbound], params), len(threeset[i]))
		// fmt.Println(i, len(threeset[i]))
		threeset[i] = Gen3Slice(bitstream[i*11 : maxbound])[:len(threeset[i])]
	}

	// fmt.Println(threeset)

	for i := 0; i < setlen; i++ {
		temp_previous := 0
		// if i == 14 && input.BlockID == 1 {
		// 	fmt.Println(threeset[i])
		// 	fmt.Println(bitstream[154:])
		// 	fmt.Println(Gen3Slice(bitstream[154:]))
		// }
		for j := 0; j < len(threeset[i]); j++ {
			thisthree := threeset[i][j]
			previous := GenPreviousSet(bitstream, temp_previous, i)
			thisC := GenNextC_Three(gc, at, Pattern, previous, input.BlockID, index, params)[thisthree]
			// if i == 14 && input.BlockID == 1 {
			// 	fmt.Println(thisthree, temp_previous, previous, string(thisC), string(GenNextC_Three(gc, at, Pattern, previous, input.BlockID, index, params)), bitstream[11*i:11*i+11])
			// }
			if thisC == 'G' || thisC == 'C' {
				gc++
			} else {
				at++
			}
			temp_previous = temp_previous << 2
			temp_previous += thisthree
			runes[index] = thisC
			Pattern = append(Pattern[1:], thisC)
			index += 1
		}
	}

	return string(runes)
}

func TempPreviousToBits(temp_previous int, max_threenum int) ([]int, bool) {
	threes := make([]int, max_threenum)
	val := temp_previous
	for i := 0; i < max_threenum; i++ {
		threes[max_threenum-i-1] = val & 3
		val = val >> 2
	}
	for len(threes) < 7 {
		threes = append(threes, 0)
	}
	return RecoverBitsFrom3Slice(threes)
}

func Encode(Data []Block, params Params) []string {
	data, _ := Readjson()
	blocknum := len(Data)
	dna := make([]string, blocknum)
	if params.Option == Gungnir_Trit_Params {
		for i := 0; i < blocknum; i++ {
			dna[i] = Block2DNA_Three(Data[i], params)
		}
		return dna
	}
	for i := 0; i < blocknum; i++ {
		dna[i] = Block2DNA(Data[i], data, params)
	}
	return dna
}

// for random error generation
func Int2Nuc(a int) rune {
	if a == 0 {
		return 'A'
	} else if a == 1 {
		return 'C'
	} else if a == 2 {
		return 'T'
	} else if a == 3 {
		return 'G'
	} else {
		return 'N'
	}
}

func Nuc2Int(a rune) int {
	if a == 'A' {
		return 0
	} else if a == 'C' {
		return 1
	} else if a == 'T' {
		return 2
	} else if a == 'G' {
		return 3
	} else {
		return -1
	}
}

func GenRevString(s string) string {
	temp := []rune(s)
	res := make([]rune, len(temp))
	for i := 0; i < len(temp); i++ {
		if temp[len(temp)-1-i] == 'A' {
			res[i] = 'T'
		} else if temp[len(temp)-1-i] == 'T' {
			res[i] = 'A'
		} else if temp[len(temp)-1-i] == 'C' {
			res[i] = 'G'
		} else if temp[len(temp)-1-i] == 'G' {
			res[i] = 'C'
		}

	}
	return string(res)

}

func Calgc(dna []string) []float64 {
	res := make([]float64, len(dna))
	for i := 0; i < len(dna); i++ {
		at := 0
		gc := 0
		temp := []rune(dna[i])
		for j := 0; j < len(temp); j++ {
			if temp[j] == 'A' || temp[j] == 'T' {
				at++
			} else {
				gc++
			}
		}
		res[i] = float64(gc) / float64(at+gc)
	}
	return res
}

func calerrorrate_site(index int, str []rune, dataset []map[string]Kmer) float64 {
	err := 0.0
	count := 0
	for i := 0; i < KmerSize; i++ {
		start := index - i
		end := start + KmerSize
		if start >= 0 && end <= len(str) {
			count += 1
			thisKmer := string(str[start:end])
			err += dataset[i][thisKmer].E_rate
		}

	}
	if count > 0 {
		err = err / float64(count)
	}
	return err
}

func Calexpecterror(dna []string) {
	_, dataset, ee := ReadjsonAll()
	sim_e := 0.0
	count := 0
	pattern := InitPattern(KmerSize - 1)

	for i := 0; i < len(dna); i++ {
		str := pattern + dna[i]
		temp := []rune(str)

		for j := KmerSize - 1; j < len(temp); j++ {
			sim_e += calerrorrate_site(j, temp, dataset)
			count += 1
		}

	}
	res := sim_e / float64(count)
	fmt.Println("error rate: ", res, "expected error rate: ", ee, "reduced", (ee-res)/ee)

}

func CalError(dna []string) []float64 {
	res := make([]float64, len(dna))
	_, dataset, _ := ReadjsonAll()
	pattern := InitPattern(KmerSize - 1)
	for i := 0; i < len(dna); i++ {
		str := pattern + dna[i]
		temp := []rune(str)
		res[i] = 0.0
		for j := KmerSize - 1; j < len(temp); j++ {
			res[i] += calerrorrate_site(j, temp, dataset)
		}
		res[i] = res[i] / float64(len(dna[i]))
	}
	return res
}
