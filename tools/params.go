// Copyright 2025 The University of Hong Kong, Department of Computer Science
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the
//     names of its contributors may be used to endorse or promote products
//     derived from this software without specific prior written permission.
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
	"log"
	"math"
)

type ParamsRaw struct {
	HashLen    int
	PayloadLen int
}

type Params struct {
	Option                  int
	MaxDepth                int
	HashLen                 int
	PayloadLen              int
	HashLenSet              []int
	PayloadLenSet           []int
	CheckPoint              []int
	HashMask0               int
	HashMask1               int
	Necessary_Decoding_Ints int
	PreviousNuc             int
	PatternMask             int
	GCrightshift            int
	Depthrightshift         int
	Indexrightshift         int
	Addressrightshift       int
	GCmax                   float64
	GCmin                   float64
	Balance_Bound           float64
	MaxBits                 int
	Seven3num               int
	Res3num                 int
	MaxHashPackage          int
}

func GenParamsRaw(hashlen, payloadlen int) ParamsRaw {
	var p ParamsRaw
	p.HashLen = hashlen
	p.PayloadLen = payloadlen
	return p
}

func (paramsraw *ParamsRaw) Compile(option int) Params {
	var params Params
	if option == Gungnir_Trit_Params {
		params.Option = Gungnir_Trit_Params
		params.HashLen = paramsraw.HashLen
		params.PayloadLen = paramsraw.PayloadLen
		params.MaxBits = paramsraw.HashLen + paramsraw.PayloadLen
		elevenbits_num := params.MaxBits / 11
		resbits := params.MaxBits - elevenbits_num*11
		params.Seven3num = elevenbits_num
		params.Res3num = 0
		if resbits > 0 {
			temp := 1
			for i := 0; i < 7; i++ {
				if temp < (1 << resbits) {
					params.Res3num += 1
					temp *= 3
				}
			}
		}
		params.MaxDepth = params.Seven3num*7 + params.Res3num
		params.MaxHashPackage = 2
		Maxhash0 := params.MaxHashPackage * 11
		if paramsraw.HashLen <= Maxhash0 {
			params.HashLenSet = make([]int, 1)
			params.HashLenSet[0] = paramsraw.HashLen
			params.PayloadLenSet = make([]int, 1)
			params.PayloadLenSet[0] = paramsraw.PayloadLen
			params.CheckPoint = make([]int, 1)
			params.CheckPoint[0] = params.MaxDepth
			params.HashMask0 = (1 << params.HashLenSet[0]) - 1
			params.HashMask1 = 0
		} else {
			params.HashLenSet = make([]int, params.Seven3num-params.MaxHashPackage+1)
			params.HashLenSet[0] = Maxhash0
			hashres := params.HashLen - Maxhash0
			hash1 := hashres / (params.Seven3num - params.MaxHashPackage)
			// fmt.Println(params.HashLen, Maxhash0, hashres, params.Seven3num-2)
			if hash1*(params.Seven3num-params.MaxHashPackage) != hashres {
				log.Fatal("Invalid Hash Length!")
			}
			params.PayloadLenSet = make([]int, params.Seven3num-params.MaxHashPackage+1)
			params.PayloadLenSet[0] = resbits
			for i := 1; i < params.Seven3num-params.MaxHashPackage+1; i++ {
				params.HashLenSet[i] = hash1
				params.PayloadLenSet[i] = 11 - hash1
			}
			params.CheckPoint = make([]int, params.Seven3num-params.MaxHashPackage+1)
			params.CheckPoint[0] = params.MaxDepth
			for i := 1; i < params.Seven3num-params.MaxHashPackage+1; i++ {
				params.CheckPoint[i] = 7*params.MaxHashPackage + 7*i
			}
			params.HashMask0 = (1 << params.HashLenSet[0]) - 1
			params.HashMask1 = (1 << hash1) - 1
		}

		// 5 * 11-bit ~ 35 * three
		params.Necessary_Decoding_Ints = (params.MaxBits + 54) / 55
		params.PreviousNuc = HomopolymerLen
		for i := 0; i < len(Exclude_motifs); i++ {
			if params.PreviousNuc < len(Exclude_motifs[i])-1 {
				params.PreviousNuc = len(Exclude_motifs[i]) - 1
			}
		}
		for i := 0; i < len(Exclude_long_motifs); i++ {
			if params.PreviousNuc < len(Exclude_long_motifs[i])-1 {
				params.PreviousNuc = len(Exclude_long_motifs[i]) - 1
			}
		}
		params.PatternMask = (1 << (2 * params.PreviousNuc)) - 1
		params.GCrightshift = 2 * params.PreviousNuc
		params.Depthrightshift = params.GCrightshift + 8
		params.Indexrightshift = params.Depthrightshift + 8
		params.Addressrightshift = params.Indexrightshift + 8
		params.GCmax = 0.575
		params.GCmin = 0.425
		params.Balance_Bound = 1.0
	} else {
		params.Option = option
		params.MaxDepth = paramsraw.HashLen + paramsraw.PayloadLen
		params.HashLen = paramsraw.HashLen
		params.PayloadLen = paramsraw.PayloadLen
		Maxhash0 := int(math.Ceil(float64(params.MaxDepth)*0.15)) + 5
		if paramsraw.HashLen <= Maxhash0 {
			params.HashLenSet = make([]int, 1)
			params.HashLenSet[0] = paramsraw.HashLen
			params.PayloadLenSet = make([]int, 1)
			params.PayloadLenSet[0] = paramsraw.PayloadLen
			params.CheckPoint = make([]int, 1)
			params.CheckPoint[0] = params.MaxDepth
			params.HashMask0 = (1 << params.HashLenSet[0]) - 1
			params.HashMask1 = 0
		} else {
			extrahashnum := (paramsraw.HashLen - Maxhash0 + 4) / 5
			r := paramsraw.HashLen - Maxhash0 - 5*(extrahashnum-1) // 1~5
			params.HashLenSet = make([]int, extrahashnum+1)
			params.HashLenSet[0] = Maxhash0
			params.HashLenSet[1] = r
			for i := 2; i < len(params.HashLenSet); i++ {
				params.HashLenSet[i] = 5
			}
			params.PayloadLenSet = make([]int, extrahashnum+1)
			payload_average := (paramsraw.PayloadLen + extrahashnum) / (extrahashnum + 1)
			if paramsraw.PayloadLen < paramsraw.HashLen {
				payload_average = (paramsraw.PayloadLen + extrahashnum - 1) / (extrahashnum)
			}
			payload_r := paramsraw.PayloadLen - payload_average*extrahashnum
			p_index := 1
			for payload_r < 0 {
				payload_r += 1
				params.PayloadLenSet[p_index] -= 1
				p_index += 1
			}
			params.PayloadLenSet[0] = payload_r
			for i := 1; i < len(params.PayloadLenSet); i++ {
				params.PayloadLenSet[i] += payload_average
			}
			params.CheckPoint = make([]int, extrahashnum+1)
			params.CheckPoint[0] = params.MaxDepth
			temp := params.HashLenSet[0]
			for i := 1; i < extrahashnum+1; i++ {
				temp += params.HashLenSet[i]
				temp += params.PayloadLenSet[i]
				params.CheckPoint[i] = temp
			}
			params.HashMask0 = (1 << params.HashLenSet[0]) - 1
			params.HashMask1 = (1 << params.HashLenSet[1]) - 1
		}
		params.Necessary_Decoding_Ints = (params.MaxDepth + 63) / 64
		params.PreviousNuc = HomopolymerLen
		for i := 0; i < len(Exclude_motifs); i++ {
			if params.PreviousNuc < len(Exclude_motifs[i])-1 {
				params.PreviousNuc = len(Exclude_motifs[i]) - 1
			}
		}
		for i := 0; i < len(Exclude_long_motifs); i++ {
			if params.PreviousNuc < len(Exclude_long_motifs[i])-1 {
				params.PreviousNuc = len(Exclude_long_motifs[i]) - 1
			}
		}
		if option == Gungnir_ONT_Params {
			if params.PreviousNuc < KmerSize-1 {
				params.PreviousNuc = KmerSize - 1
			}
		}
		params.PatternMask = (1 << (2 * params.PreviousNuc)) - 1
		params.GCrightshift = 2 * params.PreviousNuc
		params.Depthrightshift = params.GCrightshift + 8
		params.Indexrightshift = params.Depthrightshift + 9
		params.Addressrightshift = params.Indexrightshift + 9
		params.GCmax = 0.6
		params.GCmin = 0.4
		params.Balance_Bound = 1.0
		if option == Gungnir_ONT_Params {
			params.GCmax = 0.575
			params.GCmin = 0.425
			params.Balance_Bound = 0.29
		}

	}
	return params

}

const Gungnir_Default_Params = 0
const Gungnir_ONT_Params = 1
const Gungnir_Trit_Params = 2
const Primer = "TCGAAGTCAGCGTGTATTGTATG"
const KmerSize = 7
const Uint40Mask = (1 << 40) - 1
const Maxhypo_firstround = 100000
const Maxhypo_secondround = 1000000
const Maxhypo_thirdround = 5000000
const Maxhypo_forthround = 20000000
const Maxhypo_ultra = 200000000
const Maxhypo_simple = 100
const Uint5Mask = (1 << 5) - 1
const Uint8Mask = (1 << 8) - 1
const Uint9Mask = (1 << 9) - 1
const Uint11Mask = (1 << 11) - 1
const Uint14Mask = (1 << 14) - 1
const Uint16Mask = (1 << 16) - 1
const ShardCountDefaultLog = 12 // num: 1<<12 = 4096
const ShardMask = (1 << ShardCountDefaultLog) - 1
