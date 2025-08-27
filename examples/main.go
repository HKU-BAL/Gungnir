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

package main

import (
	"Gungnir/tools"
	"flag"
	"fmt"
	"math"
)

var Gungnir_density = []float64{0.5, 0.6, 0.7, 0.8, 0.9}
var Trit_density = []float64{0.87, 0.99, 1.11, 1.23, 1.35, 1.45}

func numinset(a float64, set []float64) bool {
	for i := 0; i < len(set); i++ {
		if a == set[i] {
			return true
		}
	}
	return false
}

func main() {

	information_density_ := flag.Float64("density", 0.8, "Bits/Base")
	Subrate_ := flag.Float64("sub", 0.01, "Substitution Error Rate")
	Insrate_ := flag.Float64("ins", 0.01, "Insertion Error Rate")
	Delrate_ := flag.Float64("del", 0.01, "Deletion Error Rate")
	Option_ := flag.String("option", "Gungnir", "Gungnir, Gungnir-ONT or Gungnir-Trit")
	Action_ := flag.String("action", "Encode", "Encode, AddNoise, Decode or Reconstruction")
	Input_ := flag.String("input", "../files/The Ugly Duckling", "File to be encoded")
	Output_ := flag.String("output", "../newfile", "Path for output")
	MaxSeqNum_ := flag.Int("seqnum", -1, "Maximum number of sequences allowed to be generated")
	EDmax_ := flag.Int("EDmax", 100, "Maximum Edit distance allowed")
	DecodeOption_ := flag.Bool("DecodeEDmax", true, "Whether using advancing EDmax for decoding (ignore EDmax if true)")
	threads_num1_ := flag.Int("thread1", 1, "Sequences processed in parallel")
	threads_num2_ := flag.Int("thread2", 1, "Threads for each Sequences")

	flag.Parse()

	option := *Option_
	var config, PayloadLen, HashLen int
	if option == "Gungnir" {
		config = tools.Gungnir_Default_Params
		if !numinset(*information_density_, Gungnir_density) {
			fmt.Println("Invalid Information Density!")
			return
		}
		PayloadLen = int(math.Ceil(100 * *information_density_))
		HashLen = int(100 - PayloadLen)
	} else if option == "Gungnir-ONT" {
		config = tools.Gungnir_ONT_Params
		if !numinset(*information_density_, Gungnir_density) {
			fmt.Println("Invalid Information Density!")
			return
		}
		PayloadLen = int(math.Ceil(100 * *information_density_))
		HashLen = int(100 - PayloadLen)
	} else if option == "Gungnir-Trit" {
		config = tools.Gungnir_Trit_Params
		if !numinset(*information_density_, Trit_density) {
			fmt.Println("Invalid Information Density!")
			return
		}
		PayloadLen = int(math.Ceil(100 * *information_density_))
		HashLen = int(157 - PayloadLen)
	} else {
		fmt.Println("Invalid Option!")
		return
	}

	paramsRaw := tools.GenParamsRaw(HashLen, PayloadLen)
	params := paramsRaw.Compile(config)

	sub := *Subrate_
	ins := *Insrate_
	del := *Delrate_
	err := sub + ins + del
	action := *Action_
	maxseq := *MaxSeqNum_
	edmax := *EDmax_
	thread1 := *threads_num1_
	thread2 := *threads_num2_
	input := *Input_
	output := *Output_

	if action == "Encode" {
		tools.EncodeFile(input, output, params, maxseq)
	} else if action == "AddNoise" {
		tools.AddNoise(output, del, ins, err)
	} else if action == "Decode" {
		if *DecodeOption_ {
			tools.DecodeWithEDmax(output, thread1, thread2, params)
		} else {
			tools.DecodeWithFixEDmax(output, thread1, thread2, edmax, params)
		}
	} else if action == "Reconstruction" {
		tools.ReconstructFile(output, params)
	} else {
		fmt.Println("Invalid action!")
	}

}
