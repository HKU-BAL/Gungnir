package tools

import (
	"encoding/binary"

	"github.com/spaolacci/murmur3"
)

func IntToBytes(n int) []byte {
	b := make([]byte, 8)
	binary.LittleEndian.PutUint64(b, uint64(n))
	return b
}

func PadBytes(b []byte, desiredLength int) []byte {
	paddedBytes := make([]byte, desiredLength-len(b))
	b = append(b, paddedBytes...)
	return b
}

func Ran_hash(u uint64) uint64 {
	v := u*3935559000370003845 + 2691343689449507681
	v ^= v >> 21
	v ^= v << 37
	v ^= v >> 4
	v *= 4768777513237032717
	v ^= v << 20
	v ^= v >> 41
	v ^= v << 5
	return v
}

func Ran_hash_multi(u []uint64) uint64 {
	res := uint64(0)
	for i := 0; i < len(u); i++ {
		temp := Ran_hash(u[i])
		res = res ^ temp
	}
	return res
}

func HashScale(blockID, index int, previous int) int {
	val := uint64((blockID << 48) + (index << 40) + previous)
	return int(Ran_hash(val) % 6)
}

func ReshapePayload(payload [][]int, params Params) []int {
	new_payload := make([]int, params.PayloadLen)
	count := 0
	for i := 0; i < len(payload); i++ {
		for j := 0; j < len(payload[i]); j++ {
			new_payload[count] = payload[i][j]
			count += 1
		}
	}
	return new_payload
}

func Hash0Payload(strandID int, fullpayload []int, params Params) []int {
	data := Bits2Bytes(fullpayload)
	stranddata := IntToBytes(strandID)
	data = append(data, stranddata...)
	val := int(murmur3.Sum64(data)) & params.HashMask0
	hash := make([]int, params.HashLenSet[0])
	for i := 0; i < len(hash); i++ {
		hash[i] = val & 1
		val = val >> 1
	}
	return hash
}

func HashIndex(strandID int, payloadIndex []int, hash0 []int, index int, params Params) []int {
	hash := make([]int, params.HashLenSet[index])
	data := Bits2Bytes(payloadIndex)
	stranddata := IntToBytes(strandID)
	data = append(data, stranddata...)
	hash0data := Bits2Bytes(hash0)
	data = append(data, hash0data...)
	val := int(murmur3.Sum64(data))
	if index == 1 {
		val = val & params.HashMask1
	} else {
		val = val & Uint5Mask
	}
	for i := 0; i < len(hash); i++ {
		hash[i] = val & 1
		val = val >> 1
	}
	return hash
}

func JudgeHash0(strandID int, fullpayload []int, hash0 []int, params Params) bool {
	new_hash := Hash0Payload(strandID, fullpayload, params)
	if len(new_hash) != len(hash0) {
		return false
	}
	for i := 0; i < len(new_hash); i++ {
		if hash0[i] != new_hash[i] {
			return false
		}
	}
	return true
}

func JudgeHashIndex(strandID int, payloadIndex []int, hash0 []int, hashthis []int, index int, params Params) bool {
	new_hash := HashIndex(strandID, payloadIndex, hash0, index, params)
	if len(new_hash) != len(hashthis) {
		return false
	}
	for i := 0; i < len(new_hash); i++ {
		if hashthis[i] != new_hash[i] {
			return false
		}
	}
	return true
}
