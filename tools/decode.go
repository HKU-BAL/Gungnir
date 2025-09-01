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
	"sort"
	"sync"

	"github.com/spaolacci/murmur3"
)

type Hypothesis struct {
	Info    int // lask 2k-2 bit ~ Pattern; 8 bit ~ GC; 8 bit ~ Depth; 8 bit ~ Index; 15 bit ~ sub, del, ins; 8 bit ~ last error index
	Penalty int
	Bits    []int
}

func (hypo *Hypothesis) Pattern(params Params) int {
	return hypo.Info & params.PatternMask
}

func (hypo *Hypothesis) GC(params Params) int {
	return (hypo.Info >> params.GCrightshift) & Uint8Mask
}

func (hypo *Hypothesis) Depth(params Params) int {
	return (hypo.Info >> params.Depthrightshift) & Uint8Mask
}

func (hypo *Hypothesis) Index(params Params) int {
	return (hypo.Info >> params.Indexrightshift) & Uint8Mask
}

func (hypo *Hypothesis) GenStrandID(params Params) int {
	return (hypo.Info >> params.Addressrightshift) & Uint16Mask
}

func (hypo *Hypothesis) BuildInfo(pattern, gc, depth, index, strandID int, params Params) int {
	res := 0
	res += (pattern & params.PatternMask)
	res += (gc & Uint8Mask) << params.GCrightshift
	res += (depth & Uint8Mask) << params.Depthrightshift
	res += (index & Uint8Mask) << params.Indexrightshift
	res += (strandID & Uint16Mask) << params.Addressrightshift
	return res
}

type HypoSet struct {
	hypos     []Hypothesis
	hypoMap   map[int]int
	valid_num int
}

type PlistWithLock struct {
	shards     []PlistShard
	shardCount int
}

type PlistShard struct {
	hypos []Hypothesis
}

type IDtobeDecode struct {
	IDset    []bool
	MaxID    int
	dataLock sync.Mutex
}

func (set *IDtobeDecode) Init(maxid int) {
	set.MaxID = maxid
	set.IDset = make([]bool, maxid)
	for i := 0; i < set.MaxID; i++ {
		set.IDset[i] = false
	}
}

func (set *IDtobeDecode) InitWithtempset(maxid int, tempset []int) {
	set.MaxID = maxid
	set.IDset = make([]bool, maxid)
	for i := 0; i < set.MaxID; i++ {
		set.IDset[i] = true
	}
	for i := 0; i < len(tempset); i++ {
		set.IDset[tempset[i]] = false
	}
}

func (set *IDtobeDecode) Analysis() {
	count := 0
	for i := 0; i < len(set.IDset); i++ {
		if set.IDset[i] {
			count++
		}
	}
	fmt.Println(count, "in set have been decoded!")
}

func (plist *PlistWithLock) Init(data map[string]Kmer, set *IDtobeDecode, params Params) {
	plist.shardCount = 1 << ShardCountDefaultLog
	plist.shards = make([]PlistShard, plist.shardCount)
	for i := 0; i < plist.shardCount; i++ {
		plist.shards[i].hypos = make([]Hypothesis, 0)
	}
	for i := 0; i < set.MaxID; i++ {
		if !set.IDset[i] {
			root := Genroothypo(data, i, params)
			key_id := root.CalKey(params) & ShardMask
			plist.shards[key_id].hypos = append(plist.shards[key_id].hypos, root)
		}
	}
}

type HypoShard struct {
	dataLock  sync.Mutex
	hypos     []Hypothesis
	hypoMap   map[int]int
	valid_num int
}

type HypoLock struct {
	shards     []HypoShard
	shardCount int
}

func (hypolock *HypoLock) Init() {
	hypolock.shardCount = 1 << ShardCountDefaultLog
	hypolock.shards = make([]HypoShard, hypolock.shardCount)
	for i := 0; i < hypolock.shardCount; i++ {
		hypolock.shards[i].hypos = make([]Hypothesis, 0)
		hypolock.shards[i].hypoMap = make(map[int]int)
		hypolock.shards[i].valid_num = 0
	}
}

func (hypo *Hypothesis) CalKey(params Params) int {
	index := hypo.Index(params)
	depth := hypo.Depth(params)
	Address := hypo.GenStrandID(params)
	data := make([]byte, 0)
	for i := 0; i < params.Necessary_Decoding_Ints; i++ {
		temp := IntToBytes(hypo.Bits[i])
		data = append(data, temp...)
	}
	hashpart := murmur3.Sum32(data)
	key := (index << 56) + (depth << 48) + (Address << 32) + int(hashpart)
	return key
}

func (hypo *Hypothesis) UpdateBits(thisbit int, depth int) {
	index := depth - 1
	valindex := index / 64
	bitindex := index & 63
	hypo.Bits[valindex] += (thisbit << bitindex)
}

func (hypo *Hypothesis) Selectbits(index int) int {
	valindex := index / 64
	bitindex := index & 63
	return (hypo.Bits[valindex] >> bitindex) & 1
}

func (hypo *Hypothesis) GenPrevious(params Params) int {
	res := 0
	depth := hypo.Depth(params)
	if depth < 40 {
		for i := 0; i < depth; i++ {
			res += (hypo.Selectbits(depth-1-i) << i)
		}
	} else {
		for i := 0; i < 40; i++ {
			res += (hypo.Selectbits(depth-1-i) << i)
		}
	}
	return res
}

func (hypo *Hypothesis) GenHash0(params Params) []int {
	res := make([]int, params.HashLenSet[0])

	for i := 0; i < len(res); i++ {
		res[i] = hypo.Selectbits(i)
	}
	return res
}

func (hypo *Hypothesis) GenHashIndex(hashindex int, params Params) []int {
	res := make([]int, params.HashLenSet[hashindex])
	previousLen := params.HashLenSet[0]
	for i := 1; i < hashindex; i++ {
		previousLen += params.HashLenSet[i]
		previousLen += params.PayloadLenSet[i]
	}

	for i := 0; i < len(res); i++ {
		res[i] = hypo.Selectbits(previousLen + i)
	}
	return res
}

func (hypo *Hypothesis) GenPayload0(params Params) []int {
	res := make([]int, params.PayloadLenSet[0])

	previousLen := params.HashLenSet[0]
	for i := 1; i < len(params.PayloadLenSet); i++ {
		previousLen += params.HashLenSet[i]
		previousLen += params.PayloadLenSet[i]
	}

	for i := 0; i < len(res); i++ {
		res[i] = hypo.Selectbits(previousLen + i)
	}
	return res
}

func (hypo *Hypothesis) GenPayloadIndex(payloadindex int, params Params) []int {
	res := make([]int, params.PayloadLenSet[payloadindex])

	previousLen := params.HashLenSet[0]
	for i := 1; i < payloadindex; i++ {
		previousLen += params.HashLenSet[i]
		previousLen += params.PayloadLenSet[i]
	}
	previousLen += params.HashLenSet[payloadindex]

	for i := 0; i < len(res); i++ {
		res[i] = hypo.Selectbits(previousLen + i)
	}
	return res
}

func (hypo *Hypothesis) GenPayload(params Params) [][]int {
	res := make([][]int, len(params.PayloadLenSet))

	res[0] = hypo.GenPayload0(params)

	for i := 1; i < len(res); i++ {
		res[i] = hypo.GenPayloadIndex(i, params)
	}
	return res
}

func (hypo *Hypothesis) CheckHash(params Params) bool {
	depth := hypo.Depth(params)
	for i := 1; i < len(params.CheckPoint); i++ {
		if depth == params.CheckPoint[i] {
			hash0 := hypo.GenHash0(params)
			thishash := hypo.GenHashIndex(i, params)
			thispayload := hypo.GenPayloadIndex(i, params)
			if depth == params.MaxDepth {
				return JudgeHashIndex(hypo.GenStrandID(params), thispayload, hash0, thishash, i, params) && hypo.CheckHash0(params)
			}
			return JudgeHashIndex(hypo.GenStrandID(params), thispayload, hash0, thishash, i, params)
		}
	}

	return hypo.CheckHash0(params)
}

func (hypo *Hypothesis) CheckHash0(params Params) bool {
	depth := hypo.Depth(params)
	if depth == params.MaxDepth {
		hash0 := hypo.GenHash0(params)
		payload := hypo.GenPayload(params)
		fullpayload := ReshapePayload(payload, params)
		return JudgeHash0(hypo.GenStrandID(params), fullpayload, hash0, params)
	}
	return true
}

func Pattern2runes(Pattern int, params Params) []rune {
	res := make([]rune, params.PreviousNuc)
	temp := Pattern
	for i := 0; i < params.PreviousNuc; i++ {
		res[params.PreviousNuc-1-i] = Int2Nuc(temp & 3)
		temp = temp >> 2
	}

	return res
}

func Runes2Pattern(p []rune, params Params) int {
	res := 0
	temp := 0
	for i := 0; i < params.PreviousNuc; i++ {
		res += Nuc2Int(p[params.PreviousNuc-1-i]) << temp
		temp += 2
	}

	return res
}

func Genroothypo(data map[string]Kmer, strandID int, params Params) Hypothesis {
	var hypo Hypothesis
	pattern := Runes2Pattern([]rune(InitPattern(params.PreviousNuc)), params)
	hypo.Penalty = 0.0
	hypo.Info = hypo.BuildInfo(pattern, 0, 0, Uint8Mask, strandID, params)
	hypo.Bits = make([]int, params.Necessary_Decoding_Ints)
	for i := 0; i < params.Necessary_Decoding_Ints; i++ {
		hypo.Bits[i] = 0
	}
	return hypo
}

func (hypo *Hypothesis) Addchild(DataC []rune, data map[string]Kmer, EDmax int, params Params) []Hypothesis {
	res := make([]Hypothesis, 0)
	if hypo.Depth(params) == params.MaxDepth {
		return []Hypothesis{*hypo}
	}

	previous := hypo.GenPrevious(params)
	strandID := hypo.GenStrandID(params)

	gc := hypo.GC(params)
	depth := hypo.Depth(params)
	index := hypo.Index(params)
	pattern := hypo.Pattern(params)

	potentialC := GenNextC(gc, depth-gc, Pattern2runes(pattern, params), previous, strandID, depth, data, params)

	// sub or ok ~ childX
	for i := 0; i < 2; i++ {
		if index == len(DataC)-1 {
			break
		}
		var childX Hypothesis

		x_depth := depth + 1
		x_gc := gc
		x_index := (index + 1) & Uint8Mask

		childX.Bits = make([]int, params.Necessary_Decoding_Ints)
		for j := 0; j < params.Necessary_Decoding_Ints; j++ {
			childX.Bits[j] = hypo.Bits[j]
		}

		childX.UpdateBits(i, x_depth)
		ThisC := potentialC[i] // updated

		x_pattern := ((pattern << 2) + Nuc2Int(ThisC)) & params.PatternMask

		if potentialC[i] == 'C' || potentialC[i] == 'G' {
			x_gc += 1
		}

		if ThisC == DataC[x_index] {
			childX.Penalty = hypo.Penalty
		} else {
			childX.Penalty = hypo.Penalty + 1.0 // - math.Log(data[string(append(previous, childX[i].ThisC))].X_rate)
		}

		if x_depth == params.MaxDepth {
			childX.Penalty += (len(DataC) - 1 - x_index)
			x_index = len(DataC) - 1
		}

		childX.Info = childX.BuildInfo(x_pattern, x_gc, x_depth, x_index, strandID, params)
		if childX.Penalty < EDmax+1 && childX.CheckHash(params) {
			res = append(res, childX)
		}
	}

	// del re-insert 2 possible characters
	for i := 0; i < 2; i++ {
		var childD Hypothesis

		d_depth := depth + 1
		d_index := index
		d_gc := gc

		childD.Bits = make([]int, params.Necessary_Decoding_Ints)
		for j := 0; j < params.Necessary_Decoding_Ints; j++ {
			childD.Bits[j] = hypo.Bits[j]
		}

		childD.UpdateBits(i, d_depth)
		ThisC := potentialC[i]

		d_pattern := ((pattern << 2) + Nuc2Int(ThisC)) & params.PatternMask
		childD.Penalty = hypo.Penalty + 1.0 // - math.Log(data[string(append(previous, childD[i].ThisC))].D_rate)

		if potentialC[i] == 'C' || potentialC[i] == 'G' {
			d_gc += 1
		}

		if d_depth == params.MaxDepth {
			childD.Penalty += len(DataC) - 1 - d_index
			d_index = len(DataC) - 1
		}

		childD.Info = childD.BuildInfo(d_pattern, d_gc, d_depth, d_index, strandID, params)
		if childD.Penalty < EDmax+1 && childD.CheckHash(params) {
			res = append(res, childD)
		}
	}

	// ins ignore this character
	if index != len(DataC)-1 {
		var childI Hypothesis

		i_index := index + 1

		childI.Bits = make([]int, params.Necessary_Decoding_Ints)
		for i := 0; i < params.Necessary_Decoding_Ints; i++ {
			childI.Bits[i] = hypo.Bits[i]
		}

		childI.Penalty = hypo.Penalty + 1.0 // - math.Log(data[string(append(hypo.Previousinfo, hypo.ThisC))].I_rate)

		childI.Info = childI.BuildInfo(pattern, gc, depth, i_index, strandID, params)
		if childI.Penalty < EDmax+1 {
			res = append(res, childI)
		}
	}

	return res
}

func (hypo *Hypothesis) Traceback(params Params) (b Block) {
	if hypo.Depth(params) != params.MaxDepth {
		return
	}
	b.BlockID = hypo.GenStrandID(params)
	b.Payload = hypo.GenPayload(params)
	b.Hash = make([][]int, len(params.HashLenSet))
	b.Hash[0] = hypo.GenHash0(params)
	for i := 1; i < len(b.Hash); i++ {
		b.Hash[i] = hypo.GenHashIndex(i, params)
	}

	return
}

func Updateinfo_Multithread(consensus []rune, Maxhypo int, data map[string]Kmer, EDmax int, threads_num int, set *IDtobeDecode, params Params) (b Block, suc bool) {
	var hypotree_plist PlistWithLock
	hypotree_plist.Init(data, set, params)

	for i := 0; i < len(consensus)+EDmax; i++ {
		var hypotree_qlist HypoLock
		hypotree_qlist.Init()

		var wg sync.WaitGroup
		wg.Add(hypotree_plist.shardCount)

		ch := make(chan struct{}, threads_num)

		for j := 0; j < hypotree_plist.shardCount; j++ {

			index := j
			ch <- struct{}{}
			go func() {

				num_thistask := len(hypotree_plist.shards[index].hypos)

				this_hypo := make([]Hypothesis, 0)

				for k := 0; k < num_thistask; k++ {
					temp := hypotree_plist.shards[index].hypos[k].Addchild(consensus, data, EDmax, params)
					if temp != nil {
						this_hypo = append(this_hypo, temp...)
					}
				}

				keyset := make([]int, len(this_hypo))
				for k := 0; k < len(this_hypo); k++ {
					keyset[k] = this_hypo[k].CalKey(params)
				}

				for k := 0; k < len(this_hypo); k++ {
					targetshard := keyset[k] & ShardMask
					hypotree_qlist.shards[targetshard].dataLock.Lock()

					if id, exists := hypotree_qlist.shards[targetshard].hypoMap[keyset[k]]; exists {
						if this_hypo[k].Penalty < hypotree_qlist.shards[targetshard].hypos[id].Penalty {
							hypotree_qlist.shards[targetshard].hypos[id] = this_hypo[k]
						}
					} else {
						hypotree_qlist.shards[targetshard].hypoMap[keyset[k]] = hypotree_qlist.shards[targetshard].valid_num
						hypotree_qlist.shards[targetshard].hypos = append(hypotree_qlist.shards[targetshard].hypos, this_hypo[k])
						hypotree_qlist.shards[targetshard].valid_num += 1
					}

					hypotree_qlist.shards[targetshard].dataLock.Unlock()
				}

				<-ch
				wg.Done()
			}()

		}

		wg.Wait()

		for j := 0; j < hypotree_plist.shardCount; j++ {
			hypotree_plist.shards[j].hypos = nil
		}

		totalvalidnum := 0
		for k := 0; k < hypotree_qlist.shardCount; k++ {
			totalvalidnum += hypotree_qlist.shards[k].valid_num
		}

		if totalvalidnum > Maxhypo {

			uppbound := EDmax + 2
			sumlist := make([]int, uppbound)
			for k := 0; k < uppbound; k++ {
				sumlist[k] = 0
			}

			temp_sumlist := make([][]int, hypotree_qlist.shardCount)

			var wg1 sync.WaitGroup
			wg1.Add(hypotree_qlist.shardCount)

			ch1 := make(chan struct{}, threads_num)

			for j := 0; j < hypotree_qlist.shardCount; j++ {

				index := j
				ch1 <- struct{}{}
				go func() {

					temp_sumlist[index] = make([]int, uppbound)
					for k := 0; k < uppbound; k++ {
						temp_sumlist[index][k] = 0
					}
					for k := 0; k < hypotree_qlist.shards[index].valid_num; k++ {
						val := hypotree_qlist.shards[index].hypos[k].Penalty
						temp_sumlist[index][val] += 1
					}

					<-ch1
					wg1.Done()
				}()

			}

			wg1.Wait()

			for j := 0; j < hypotree_qlist.shardCount; j++ {
				for k := 0; k < uppbound; k++ {
					sumlist[k] += temp_sumlist[j][k]
				}
			}

			count := 0
			valid_index := 0
			for k := 0; k < uppbound; k++ {
				count += sumlist[k]
				if count > Maxhypo {
					// count -= sumlist[k]
					valid_index = k - 1
					break
				}
			}

			maxpenal := valid_index

			var wg2 sync.WaitGroup
			wg2.Add(hypotree_qlist.shardCount)

			ch2 := make(chan struct{}, threads_num)

			for j := 0; j < hypotree_qlist.shardCount; j++ {

				index := j
				ch2 <- struct{}{}
				go func() {
					count_index := 0
					for k := 0; k <= valid_index; k++ {
						count_index += temp_sumlist[index][k]
					}
					hypotree_plist.shards[index].hypos = make([]Hypothesis, count_index)
					iterator := 0

					for k := 0; k < hypotree_qlist.shards[index].valid_num; k++ {
						newbits := make([]int, params.Necessary_Decoding_Ints)
						for l := 0; l < params.Necessary_Decoding_Ints; l++ {
							newbits[l] = hypotree_qlist.shards[index].hypos[k].Bits[l]
						}
						if hypotree_qlist.shards[index].hypos[k].Penalty <= maxpenal {
							hypotree_plist.shards[index].hypos[iterator] = Hypothesis{Info: hypotree_qlist.shards[index].hypos[k].Info,
								Penalty: hypotree_qlist.shards[index].hypos[k].Penalty,
								Bits:    newbits,
							}
							iterator += 1
						}
					}
					hypotree_qlist.shards[index].hypoMap = nil
					hypotree_qlist.shards[index].hypos = nil

					<-ch2
					wg2.Done()
				}()

			}

			wg2.Wait()

		} else {
			for j := 0; j < hypotree_plist.shardCount; j++ {
				hypotree_plist.shards[j].hypos = hypotree_qlist.shards[j].hypos
			}
		}

	}

	final_hypos := make([]Hypothesis, 0)
	for i := 0; i < hypotree_plist.shardCount; i++ {
		final_hypos = append(final_hypos, hypotree_plist.shards[i].hypos...)
	}

	sort.Slice(final_hypos, func(a, b int) bool {
		return final_hypos[a].Penalty < final_hypos[b].Penalty
	})

	suc = false

	for i := 0; i < len(final_hypos); i++ {
		if final_hypos[i].Depth(params) == params.MaxDepth {
			strandID := final_hypos[i].GenStrandID(params)
			set.dataLock.Lock()
			set.IDset[strandID] = true
			set.dataLock.Unlock()
			b = final_hypos[i].Traceback(params)
			suc = true
			return
		}
	}

	return
}

func Updateinfo_Singlethread(consensus []rune, Maxhypo int, data map[string]Kmer, EDmax int, set *IDtobeDecode, params Params) (b Block, suc bool) {

	hypotree_plist := make([]Hypothesis, 0)
	set.dataLock.Lock()
	for i := 0; i < set.MaxID; i++ {
		if !set.IDset[i] {
			root := Genroothypo(data, i, params)
			hypotree_plist = append(hypotree_plist, root)
		}
	}
	set.dataLock.Unlock()

	for i := 0; i < len(consensus)+EDmax; i++ {
		hypotree_qlist := &HypoSet{
			hypos:     make([]Hypothesis, 5*len(hypotree_plist)),
			hypoMap:   make(map[int]int),
			valid_num: 0,
		}

		task := len(hypotree_plist)

		for j := 0; j < task; j++ {
			temp := hypotree_plist[j].Addchild(consensus, data, EDmax, params)
			if temp != nil {
				for k := 0; k < len(temp); k++ {
					key := temp[k].CalKey(params)
					if id, exists := hypotree_qlist.hypoMap[key]; exists {
						if temp[k].Penalty < hypotree_qlist.hypos[id].Penalty {
							hypotree_qlist.hypos[id] = temp[k]
						}
					} else {
						hypotree_qlist.hypoMap[key] = hypotree_qlist.valid_num
						hypotree_qlist.hypos[hypotree_qlist.valid_num] = temp[k]
						hypotree_qlist.valid_num += 1
					}
				}
			}

		}

		hypotree_plist = nil

		if hypotree_qlist.valid_num > Maxhypo {

			uppbound := EDmax + 2
			sumlist := make([]int, uppbound)
			for k := 0; k < uppbound; k++ {
				sumlist[k] = 0
			}

			for k := 0; k < hypotree_qlist.valid_num; k++ {
				val := hypotree_qlist.hypos[k].Penalty
				sumlist[val] += 1
			}

			count := 0
			valid_index := 0
			for k := 0; k < uppbound; k++ {
				count += sumlist[k]
				if count > Maxhypo {
					count -= sumlist[k]
					valid_index = k - 1
					if k == 0 && hypotree_qlist.valid_num > 0 {
						valid_index = 0
						count = sumlist[0]
					}
					break
				}
			}

			hypotree_plist = make([]Hypothesis, count)
			iterator := 0
			maxpenal := valid_index

			for k := 0; k < hypotree_qlist.valid_num; k++ {
				newbits := make([]int, params.Necessary_Decoding_Ints)
				for l := 0; l < params.Necessary_Decoding_Ints; l++ {
					newbits[l] = hypotree_qlist.hypos[k].Bits[l]
				}
				if hypotree_qlist.hypos[k].Penalty <= maxpenal {
					hypotree_plist[iterator] = Hypothesis{
						Info: hypotree_qlist.hypos[k].Info, Penalty: hypotree_qlist.hypos[k].Penalty,
						Bits: newbits,
					}
					iterator += 1
				}
			}
		} else {
			hypotree_plist = hypotree_qlist.hypos[:hypotree_qlist.valid_num]
		}

	}

	sort.Slice(hypotree_plist, func(a, b int) bool {
		return hypotree_plist[a].Penalty < hypotree_plist[b].Penalty
	})

	suc = false

	for i := 0; i < len(hypotree_plist); i++ {
		if hypotree_plist[i].Depth(params) == params.MaxDepth {
			strandID := hypotree_plist[i].GenStrandID(params)
			set.dataLock.Lock()
			set.IDset[strandID] = true
			set.dataLock.Unlock()
			b = hypotree_plist[i].Traceback(params)
			suc = true
			return
		}
	}

	return
}

func Decode_Multithread(Consensus []string, Maxhypo int, threads_num int, EDmax int, set *IDtobeDecode, params Params) ([]Block, []bool) {
	data, _ := Readjson()
	blocknum := len(Consensus)
	bit_stream := make([]Block, blocknum)
	decode_res := make([]bool, blocknum)

	// fmt.Println("Decoding... Threads number: " + strconv.Itoa(threads_num))

	for i := 0; i < blocknum; i++ {
		var temp Block
		temp, decode_res[i] = Updateinfo_Multithread([]rune(Consensus[i]), Maxhypo, data, EDmax, threads_num, set, params)
		if decode_res[i] {
			bit_stream[i] = temp
		}

	}
	return bit_stream, decode_res
}

func Decode_Parallel(Consensus []string, Maxhypo int, threads_num int, EDmax int, set *IDtobeDecode, params Params) ([]Block, []bool) {
	data, _ := Readjson()
	blocknum := len(Consensus)
	bit_stream := make([]Block, blocknum)
	decode_res := make([]bool, blocknum)

	var wg sync.WaitGroup
	wg.Add(blocknum)

	ch := make(chan struct{}, threads_num)

	for i := 0; i < blocknum; i++ {

		index := i
		ch <- struct{}{}
		go func() {

			var temp Block
			temp, decode_res[index] = Updateinfo_Singlethread([]rune(Consensus[index]), Maxhypo, data, EDmax, set, params)
			if decode_res[index] {
				bit_stream[index] = temp
			}

			<-ch
			wg.Done()
		}()

	}

	wg.Wait()
	return bit_stream, decode_res
}

func Decode_Mix(Consensus []string, Maxhypo int, threads_num1 int, threads_num2 int, EDmax int, set *IDtobeDecode, params Params) ([]Block, []bool) {
	data, _ := Readjson()
	blocknum := len(Consensus)
	bit_stream := make([]Block, blocknum)
	decode_res := make([]bool, blocknum)

	var wg sync.WaitGroup
	wg.Add(blocknum)

	ch := make(chan struct{}, threads_num1)

	for i := 0; i < blocknum; i++ {

		index := i
		ch <- struct{}{}
		go func() {

			var temp Block
			temp, decode_res[index] = Updateinfo_Multithread([]rune(Consensus[index]), Maxhypo, data, EDmax, threads_num2, set, params)
			if decode_res[index] {
				bit_stream[index] = temp
			}

			<-ch
			wg.Done()
		}()

	}

	wg.Wait()
	return bit_stream, decode_res
}
