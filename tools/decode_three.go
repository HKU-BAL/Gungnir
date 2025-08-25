package tools

import (
	"fmt"
	"sort"
	"strconv"
	"sync"

	"github.com/spaolacci/murmur3"
)

type Hypothesis_Three struct {
	Info         int
	Penalty_Info int // temp_previous 14 bits; Penalty 16 bits
	Bits         []int
}

func (hypo *Hypothesis_Three) Pattern(params Params) int {
	return hypo.Info & params.PatternMask
}

func (hypo *Hypothesis_Three) GC(params Params) int {
	return (hypo.Info >> params.GCrightshift) & Uint8Mask
}

func (hypo *Hypothesis_Three) Depth(params Params) int {
	return (hypo.Info >> params.Depthrightshift) & Uint8Mask
}

func (hypo *Hypothesis_Three) Index(params Params) int {
	return (hypo.Info >> params.Indexrightshift) & Uint8Mask
}

func (hypo *Hypothesis_Three) GenStrandID(params Params) int {
	return (hypo.Info >> params.Addressrightshift) & Uint16Mask
}

func (hypo *Hypothesis_Three) BuildInfo(pattern, gc, depth, index, strandID int, params Params) int {
	res := 0
	res += (pattern & params.PatternMask)
	res += (gc & Uint8Mask) << params.GCrightshift
	res += (depth & Uint8Mask) << params.Depthrightshift
	res += (index & Uint8Mask) << params.Indexrightshift
	res += (strandID & Uint16Mask) << params.Addressrightshift
	return res
}

func (hypo *Hypothesis_Three) Penalty() int {
	return hypo.Penalty_Info & Uint16Mask
}

func (hypo *Hypothesis_Three) TempPrevious() int {
	return (hypo.Penalty_Info >> 16) & Uint14Mask
}

func (hypo *Hypothesis_Three) BuildPenaltyInfo(penalty, temp_previous int, params Params) int {
	res := 0
	res += (penalty & Uint16Mask)
	res += (temp_previous & Uint14Mask) << 16
	return res
}

type HypoSet_Three struct {
	hypos     []Hypothesis_Three
	hypoMap   map[int]int
	valid_num int
}

type PlistWithLock_Three struct {
	shards     []PlistShard_Three
	shardCount int
}

type PlistShard_Three struct {
	hypos []Hypothesis_Three
}

func Genroothypo_Three(strandID int, params Params) Hypothesis_Three {
	var hypo Hypothesis_Three
	pattern := Runes2Pattern([]rune(InitPattern(params.PreviousNuc)), params)
	hypo.Penalty_Info = 0
	hypo.Info = hypo.BuildInfo(pattern, 0, 0, Uint8Mask, strandID, params)
	hypo.Bits = make([]int, params.Necessary_Decoding_Ints)
	for i := 0; i < params.Necessary_Decoding_Ints; i++ {
		hypo.Bits[i] = 0
	}
	return hypo
}

func (hypo *Hypothesis_Three) CalKey(params Params) int {
	index := hypo.Index(params)
	depth := hypo.Depth(params)
	Address := hypo.GenStrandID(params)
	data := make([]byte, 0)
	for i := 0; i < params.Necessary_Decoding_Ints; i++ {
		temp := IntToBytes(hypo.Bits[i])
		data = append(data, temp...)
	}
	temp_pre := hypo.TempPrevious()
	temp_pre_bytes := IntToBytes(temp_pre)
	data = append(data, temp_pre_bytes...)
	hashpart := murmur3.Sum32(data)
	key := (index << 56) + (depth << 48) + (Address << 32) + int(hashpart)
	return key
}

func (plist *PlistWithLock_Three) Init(set *IDtobeDecode, params Params) {
	plist.shardCount = 1 << ShardCountDefaultLog
	plist.shards = make([]PlistShard_Three, plist.shardCount)
	for i := 0; i < plist.shardCount; i++ {
		plist.shards[i].hypos = make([]Hypothesis_Three, 0)
	}
	for i := 0; i < set.MaxID; i++ {
		if !set.IDset[i] {
			root := Genroothypo_Three(i, params)
			key_id := root.CalKey(params) & ShardMask
			plist.shards[key_id].hypos = append(plist.shards[key_id].hypos, root)
		}
	}
}

type HypoShard_Three struct {
	dataLock  sync.Mutex
	hypos     []Hypothesis_Three
	hypoMap   map[int]int
	valid_num int
}

type HypoLock_Three struct {
	shards     []HypoShard_Three
	shardCount int
}

func (hypolock *HypoLock_Three) Init() {
	hypolock.shardCount = 1 << ShardCountDefaultLog
	hypolock.shards = make([]HypoShard_Three, hypolock.shardCount)
	for i := 0; i < hypolock.shardCount; i++ {
		hypolock.shards[i].hypos = make([]Hypothesis_Three, 0)
		hypolock.shards[i].hypoMap = make(map[int]int)
		hypolock.shards[i].valid_num = 0
	}
}

func (hypo *Hypothesis_Three) UpdateBits(thisbit int, index int) {
	valindex := index / 55
	bitindex := index % 55
	hypo.Bits[valindex] += (thisbit << bitindex)
}

func (hypo *Hypothesis_Three) Selectbits(index int) int {
	valindex := index / 55
	bitindex := index % 55
	return (hypo.Bits[valindex] >> bitindex) & 1
}

func (hypo *Hypothesis_Three) GenPrevious(depth int, params Params) []int {
	bit_stream := make([]int, params.MaxBits)
	bitnum := depth * 11 / 7
	for i := 0; i < bitnum; i++ {
		bit_stream[i] = hypo.Selectbits(i)
	}
	return GenPreviousSet(bit_stream, hypo.TempPrevious(), depth/7)
}

func (hypo *Hypothesis_Three) GenHash0(params Params) []int {
	res := make([]int, params.HashLenSet[0])

	for i := 0; i < len(res); i++ {
		res[i] = hypo.Selectbits(i)
	}
	return res
}

func (hypo *Hypothesis_Three) GenHashIndex(hashindex int, params Params) []int {
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

func (hypo *Hypothesis_Three) GenPayload0(params Params) []int {
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

func (hypo *Hypothesis_Three) GenPayloadIndex(payloadindex int, params Params) []int {
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

func (hypo *Hypothesis_Three) GenPayload(params Params) [][]int {
	res := make([][]int, len(params.PayloadLenSet))

	res[0] = hypo.GenPayload0(params)

	for i := 1; i < len(res); i++ {
		res[i] = hypo.GenPayloadIndex(i, params)
	}
	return res
}

func (hypo *Hypothesis_Three) CheckHash(params Params) bool {
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

func (hypo *Hypothesis_Three) CheckHash0(params Params) bool {
	depth := hypo.Depth(params)
	if depth == params.MaxDepth {
		hash0 := hypo.GenHash0(params)
		payload := hypo.GenPayload(params)
		fullpayload := ReshapePayload(payload, params)
		return JudgeHash0(hypo.GenStrandID(params), fullpayload, hash0, params)
	}
	return true
}

func (hypo *Hypothesis_Three) Addchild(DataC []rune, params Params) []Hypothesis_Three {
	res := make([]Hypothesis_Three, 0)
	if hypo.Depth(params) == params.MaxDepth {
		return []Hypothesis_Three{*hypo}
	}

	gc := hypo.GC(params)
	depth := hypo.Depth(params)
	index := hypo.Index(params)
	pattern := hypo.Pattern(params)

	previous := hypo.GenPrevious(depth, params)
	strandID := hypo.GenStrandID(params)

	penalty := hypo.Penalty()
	temp_previous := hypo.TempPrevious()

	potentialC := GenNextC_Three(gc, depth-gc, Pattern2runes(pattern, params), previous, strandID, depth, params)

	// To be done: set temp_previous to 0 and update bits if depth of child = 7, 8, ...

	// sub or ok ~ childX
	for i := 0; i < 3; i++ {
		if index == len(DataC)-1 {
			break
		}
		var childX Hypothesis_Three

		x_depth := depth + 1
		x_gc := gc
		x_index := (index + 1) & Uint8Mask

		childX.Bits = make([]int, params.Necessary_Decoding_Ints)
		for j := 0; j < params.Necessary_Decoding_Ints; j++ {
			childX.Bits[j] = hypo.Bits[j]
		}

		// childX.UpdateBits(i, x_index)
		ThisC := potentialC[i]

		x_pattern := ((pattern << 2) + Nuc2Int(ThisC)) & params.PatternMask

		if potentialC[i] == 'C' || potentialC[i] == 'G' {
			x_gc += 1
		}

		x_penalty := penalty
		x_temp_previous := temp_previous << 2
		x_temp_previous += i
		x_temp_previous = x_temp_previous & Uint14Mask

		if ThisC != DataC[x_index] {
			x_penalty += 1
		}

		if x_depth == params.MaxDepth {
			x_penalty += (len(DataC) - 1 - x_index)
			x_index = len(DataC) - 1
		}

		childX.Info = childX.BuildInfo(x_pattern, x_gc, x_depth, x_index, strandID, params)

		if x_depth%7 == 0 {
			temp_bits, valid := TempPreviousToBits(x_temp_previous, 7)
			if !valid {
				continue
			}
			// if x_penalty == 0 && x_depth == 7 {
			// 	fmt.Println(temp_bits)
			// }
			package_index := (x_depth / 7) - 1
			for j := 0; j < 11; j++ {
				childX.UpdateBits(temp_bits[j], package_index*11+j)
			}
			x_temp_previous = 0
		} else if x_depth == params.MaxDepth {
			temp_bits, valid := TempPreviousToBits(x_temp_previous, params.Res3num)
			if !valid {
				continue
			}
			temp_bits = temp_bits[:params.MaxBits%11]
			package_index := x_depth / 7
			for j := 0; j < len(temp_bits); j++ {
				childX.UpdateBits(temp_bits[j], package_index*11+j)
			}

			x_temp_previous = 0
		}

		childX.Penalty_Info = childX.BuildPenaltyInfo(x_penalty, x_temp_previous, params)

		if x_penalty < params.Penalty_upperbound && childX.CheckHash(params) {
			res = append(res, childX)
		}

		// if hypo.Depth(params) == 99 && hypo.Penalty() == 0 && string(DataC) == "CACGATCCCGCTAGGGATTATCGCACTGTCGACCCTTAGCGTGGTAAGTTTGACTGTCACCTCCTAACGAAAGAAACGAAACGTGGTAGGTTCACTCATA" {
		// 	fmt.Println(childX.GenStrandID(params), childX.Depth(params), childX.TempPrevious(), previous, childX.Penalty())
		// 	temp_bits := TempPreviousToBits(temp_fake, params.Res3num)[:params.MaxBits%11]
		// 	for j := 0; j < len(temp_bits); j++ {
		// 		fmt.Println(temp_fake, temp_bits)
		// 	}
		// }
	}

	// del re-insert 3 possible characters
	for i := 0; i < 3; i++ {
		var childD Hypothesis_Three

		d_depth := depth + 1
		d_index := index
		d_gc := gc

		childD.Bits = make([]int, params.Necessary_Decoding_Ints)
		for j := 0; j < params.Necessary_Decoding_Ints; j++ {
			childD.Bits[j] = hypo.Bits[j]
		}

		// childD.UpdateBits(i, d_depth)
		ThisC := potentialC[i]

		d_pattern := ((pattern << 2) + Nuc2Int(ThisC)) & params.PatternMask

		d_penalty := penalty + 1
		d_temp_previous := temp_previous << 2
		d_temp_previous += i
		d_temp_previous = d_temp_previous & Uint14Mask

		if potentialC[i] == 'C' || potentialC[i] == 'G' {
			d_gc += 1
		}

		if d_depth == params.MaxDepth {
			d_penalty += len(DataC) - 1 - d_index
			d_index = len(DataC) - 1
		}

		childD.Info = childD.BuildInfo(d_pattern, d_gc, d_depth, d_index, strandID, params)

		if d_depth%7 == 0 {
			temp_bits, valid := TempPreviousToBits(d_temp_previous, 7)
			if !valid {
				continue
			}
			package_index := (d_depth / 7) - 1
			for j := 0; j < 11; j++ {
				childD.UpdateBits(temp_bits[j], package_index*11+j)
			}
			if d_depth != params.MaxDepth {
				d_temp_previous = 0
			}
		} else if d_depth == params.MaxDepth {
			temp_bits, valid := TempPreviousToBits(d_temp_previous, params.Res3num)
			if !valid {
				continue
			}
			temp_bits = temp_bits[:params.MaxBits%11]
			package_index := d_depth / 7
			for j := 0; j < len(temp_bits); j++ {
				childD.UpdateBits(temp_bits[j], package_index*11+j)
			}
		}

		childD.Penalty_Info = childD.BuildPenaltyInfo(d_penalty, d_temp_previous, params)

		if d_penalty < params.Penalty_upperbound && childD.CheckHash(params) {
			res = append(res, childD)
		}
	}

	// ins ignore this character
	if index != len(DataC)-1 {
		var childI Hypothesis_Three

		i_index := index + 1

		childI.Bits = make([]int, params.Necessary_Decoding_Ints)
		for i := 0; i < params.Necessary_Decoding_Ints; i++ {
			childI.Bits[i] = hypo.Bits[i]
		}

		childI.Penalty_Info = hypo.Penalty_Info + 1

		childI.Info = childI.BuildInfo(pattern, gc, depth, i_index, strandID, params)
		if childI.Penalty() < params.Penalty_upperbound {
			res = append(res, childI)
		}
	}

	// k := 98

	// if hypo.Depth(params) == k && hypo.Penalty() == 0 && string(DataC) == "CACGATCCCGCTAGGGATTATCGCACTGTCGACCCTTAGCGTGGTAAGTTTGACTGTCACCTCCTAACGAAAGAAACGAAACGTGGTAGGTTCACTCATA" {
	// 	fmt.Println(k, string(potentialC), hypo.TempPrevious())
	// 	for i := 0; i < len(res); i++ {
	// 		fmt.Println(res[i].GenStrandID(params), res[i].Depth(params), res[i].TempPrevious(), previous)
	// 	}
	// }

	return res
}

func (hypo *Hypothesis_Three) Traceback(params Params) (b Block) {
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

func Updateinfo_Three_Multithread(consensus []rune, Maxhypo int, threads_num int, set *IDtobeDecode, params Params) (b Block, suc bool) {
	var hypotree_plist PlistWithLock_Three
	hypotree_plist.Init(set, params)

	for i := 0; i < len(consensus)+params.MaxError; i++ {
		var hypotree_qlist HypoLock_Three
		hypotree_qlist.Init()

		var wg sync.WaitGroup
		wg.Add(hypotree_plist.shardCount)

		ch := make(chan struct{}, threads_num)

		for j := 0; j < hypotree_plist.shardCount; j++ {

			index := j
			ch <- struct{}{}
			go func() {

				num_thistask := len(hypotree_plist.shards[index].hypos)

				this_hypo := make([]Hypothesis_Three, 0)

				for k := 0; k < num_thistask; k++ {
					temp := hypotree_plist.shards[index].hypos[k].Addchild(consensus, params)
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
						if this_hypo[k].Penalty() < hypotree_qlist.shards[targetshard].hypos[id].Penalty() {
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

			uppbound := params.Penalty_upperbound + 1
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
						val := hypotree_qlist.shards[index].hypos[k].Penalty()
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
					hypotree_plist.shards[index].hypos = make([]Hypothesis_Three, count_index)
					iterator := 0

					for k := 0; k < hypotree_qlist.shards[index].valid_num; k++ {
						newbits := make([]int, params.Necessary_Decoding_Ints)
						for l := 0; l < params.Necessary_Decoding_Ints; l++ {
							newbits[l] = hypotree_qlist.shards[index].hypos[k].Bits[l]
						}
						if hypotree_qlist.shards[index].hypos[k].Penalty() <= maxpenal {
							hypotree_plist.shards[index].hypos[iterator] = Hypothesis_Three{Info: hypotree_qlist.shards[index].hypos[k].Info,
								Penalty_Info: hypotree_qlist.shards[index].hypos[k].Penalty_Info,
								Bits:         newbits,
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

	final_hypos := make([]Hypothesis_Three, 0)
	for i := 0; i < hypotree_plist.shardCount; i++ {
		final_hypos = append(final_hypos, hypotree_plist.shards[i].hypos...)
	}

	sort.Slice(final_hypos, func(a, b int) bool {
		return final_hypos[a].Penalty() < final_hypos[b].Penalty()
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

func Updateinfo_Three_Singlethread(consensus []rune, Maxhypo int, set *IDtobeDecode, params Params) (b Block, suc bool) {

	hypotree_plist := make([]Hypothesis_Three, 0)
	set.dataLock.Lock()
	for i := 0; i < set.MaxID; i++ {
		if !set.IDset[i] {
			root := Genroothypo_Three(i, params)
			hypotree_plist = append(hypotree_plist, root)
		}
	}
	set.dataLock.Unlock()

	for i := 0; i < len(consensus)+params.MaxError; i++ {
		hypotree_qlist := &HypoSet_Three{
			hypos:     make([]Hypothesis_Three, 7*len(hypotree_plist)),
			hypoMap:   make(map[int]int),
			valid_num: 0,
		}

		task := len(hypotree_plist)

		for j := 0; j < task; j++ {
			temp := hypotree_plist[j].Addchild(consensus, params)
			if temp != nil {
				for k := 0; k < len(temp); k++ {
					key := temp[k].CalKey(params)
					if id, exists := hypotree_qlist.hypoMap[key]; exists {
						if temp[k].Penalty() < hypotree_qlist.hypos[id].Penalty() {
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

			uppbound := params.Penalty_upperbound + 1
			sumlist := make([]int, uppbound)
			for k := 0; k < uppbound; k++ {
				sumlist[k] = 0
			}

			for k := 0; k < hypotree_qlist.valid_num; k++ {
				val := hypotree_qlist.hypos[k].Penalty()
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

			hypotree_plist = make([]Hypothesis_Three, count)
			iterator := 0
			maxpenal := valid_index

			for k := 0; k < hypotree_qlist.valid_num; k++ {
				newbits := make([]int, params.Necessary_Decoding_Ints)
				for l := 0; l < params.Necessary_Decoding_Ints; l++ {
					newbits[l] = hypotree_qlist.hypos[k].Bits[l]
				}
				if hypotree_qlist.hypos[k].Penalty() <= maxpenal {
					hypotree_plist[iterator] = Hypothesis_Three{
						Info: hypotree_qlist.hypos[k].Info, Penalty_Info: hypotree_qlist.hypos[k].Penalty_Info,
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
		return hypotree_plist[a].Penalty() < hypotree_plist[b].Penalty()
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

func Decode_Three_Multithread(Consensus []string, Maxhypo int, threads_num int, maxid int, set *IDtobeDecode, params Params) ([]Block, []bool) {
	blocknum := len(Consensus)
	bit_stream := make([]Block, blocknum)
	decode_res := make([]bool, blocknum)

	// fmt.Println("Decoding... Threads number: " + strconv.Itoa(threads_num))

	for i := 0; i < blocknum; i++ {
		var temp Block
		temp, decode_res[i] = Updateinfo_Three_Multithread([]rune(Consensus[i]), Maxhypo, threads_num, set, params)
		if decode_res[i] {
			bit_stream[i] = temp
		}

	}
	return bit_stream, decode_res
}

func Decode_Parallel_Three_Raw(Consensus []string, threads_num int, maxid int, params Params) ([]Block, []bool) {
	set := &IDtobeDecode{}
	set.Init(maxid)
	blocknum := len(Consensus)
	bit_stream := make([]Block, blocknum)
	decode_res := make([]bool, blocknum)

	fmt.Println("Decoding... Threads number: " + strconv.Itoa(threads_num))

	var wg sync.WaitGroup
	wg.Add(blocknum)

	ch := make(chan struct{}, threads_num)

	for i := 0; i < blocknum; i++ {

		index := i
		ch <- struct{}{}
		go func() {

			var temp Block
			temp, decode_res[index] = Updateinfo_Three_Singlethread([]rune(Consensus[index]), Maxhypo_simple, set, params)
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

func Decode_Three_Parallel(Consensus []string, Maxhypo int, threads_num int, set *IDtobeDecode, params Params) ([]Block, []bool) {
	blocknum := len(Consensus)
	bit_stream := make([]Block, blocknum)
	decode_res := make([]bool, blocknum)

	fmt.Println("Decoding... Threads number: " + strconv.Itoa(threads_num))

	var wg sync.WaitGroup
	wg.Add(blocknum)

	ch := make(chan struct{}, threads_num)

	for i := 0; i < blocknum; i++ {

		index := i
		ch <- struct{}{}
		go func() {

			var temp Block
			temp, decode_res[index] = Updateinfo_Three_Singlethread([]rune(Consensus[index]), Maxhypo, set, params)
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

func Decode_Three_Mix(Consensus []string, Maxhypo int, threads_num1 int, threads_num2 int, set *IDtobeDecode, params Params) ([]Block, []bool) {
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
			temp, decode_res[index] = Updateinfo_Three_Multithread([]rune(Consensus[index]), Maxhypo, threads_num2, set, params)
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
