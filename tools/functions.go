package tools

import (
	"fmt"
	"os"
	"runtime"
)

func Genfilename(filepath string) (string, string, string) {
	var Origin = "/Origin"
	var Error = "/Add_Error"
	var Decode = "/Decoded"

	Origin = filepath + Origin
	Error = filepath + Error
	Decode = filepath + Decode
	return Origin, Error, Decode
}

func EncodeFile(inputfile string, outputpath string, params Params, maximumseq int) {

	os.MkdirAll(outputpath, 0755)

	Origin_Name, _, _ := Genfilename(outputpath)

	file := Readfile(inputfile)

	b := File2Bits(file)

	d := BitsintoBlocks(b, params)

	if maximumseq > 0 && maximumseq < len(d) {
		d = d[:maximumseq]
	}

	fmt.Println("Strand Num: ", len(d))

	ori_seqs := Encode(d, params)

	GenFasta(ori_seqs, Origin_Name)
}

// Sub = Total - Del - Ins
func AddNoise(filepath string, Delrate, Insrate, ErrorRate float64) {
	Origin_Name, Error_Name, _ := Genfilename(filepath)

	ori_seqs := ReadFasta(Origin_Name)
	Subrate := ErrorRate - Delrate - Insrate

	seqs := AddError(ori_seqs, Subrate, Delrate, Insrate, ErrorRate) // Fixed Error Rate

	GenFasta(seqs, Error_Name)

}

func AnalysisFailure(set *IDtobeDecode, filepath string, params Params) (int, int) {

	Origin_Name, _, Decode_Name := Genfilename(filepath)
	seqs := ReadFasta(Origin_Name)
	dec_seqs := ReadFasta(Decode_Name)

	data, _ := Readjson()
	var zero_block Block
	var zero_seq string

	if params.Option != Gungnir_Trit_Params {
		zero_seq = Block2DNA(zero_block, data, params)
	} else {
		zero_seq = Block2DNA_Three(zero_block, params)
	}

	fail := 0
	tempset := make([]int, 0)

	for i := 0; i < len(seqs); i++ {
		if Hamming_Distrance([]rune(dec_seqs[i]), []rune(zero_seq)) == 0 {
			tempset = append(tempset, i)
			fail++
		}
	}
	set.InitWithtempset(len(seqs), tempset)
	fmt.Println("Total Sequences:", len(seqs), " Failure Sequnces:", fail)
	return fail, len(seqs)
}

func AnalysisAll(filepath string, params Params) (float64, int, int, int) {
	if params.Option != Gungnir_Trit_Params {
		Origin_Name, _, Decode_Name := Genfilename(filepath)
		seqs := ReadFasta(Origin_Name)
		dec_seqs := ReadFasta(Decode_Name)

		data, _ := Readjson()
		var zero_block Block
		zero_seq := Block2DNA(zero_block, data, params)

		res := make([]int, 0)
		count := 0
		fail := 0
		for i := 0; i < len(seqs); i++ {
			val := Hamming_Distrance([]rune(dec_seqs[i]), []rune(seqs[i]))
			res = append(res, val)
		}

		// tempset := make([]int, 0)
		for i := 0; i < len(res); i++ {
			if res[i] != 0 {
				if dec_seqs[i] == zero_seq {
					// tempset = append(tempset, i)
					fail++
				}
				count++
			}
		}
		suc_rate := 1.0 - (float64(count) / float64(len(seqs)))

		fmt.Println("Total Sequences:", len(seqs), " Failure Sequnces:", fail, " Mistaken Sequences:", count-fail)
		return suc_rate, fail, count - fail, len(seqs)
	} else {
		Origin_Name, _, Decode_Name := Genfilename(filepath)
		seqs := ReadFasta(Origin_Name)
		dec_seqs := ReadFasta(Decode_Name)
		// err_seqs := ReadFasta(Error_Name)
		res := make([]int, 0)
		count := 0
		fail := 0
		for i := 0; i < len(seqs); i++ {
			val := Hamming_Distrance([]rune(dec_seqs[i]), []rune(seqs[i]))
			res = append(res, val)
		}
		var zero_block Block
		zero_seq := Block2DNA_Three(zero_block, params)
		// tempset := make([]int, 0)
		for i := 0; i < len(res); i++ {
			if res[i] != 0 {
				if dec_seqs[i] == zero_seq {
					// tempset = append(tempset, i)
					fail++
				}
				count++
			}
		}
		suc_rate := 1.0 - (float64(count) / float64(len(res)))
		fmt.Println("Total Sequences:", len(seqs), " Failure Sequnces:", fail, " Mistaken Sequences:", count-fail)
		return suc_rate, fail, count - fail, len(res)
	}

}

// threads_num1: How many sequences in parallel; thread_num2: How many threads for each sequence
func FirstDecode(filepath string, threads_num1, threads_num2 int, Hmax int, EDmax int, params Params) (int, int) {

	var deco []Block
	var deco_suc []bool

	set := &IDtobeDecode{}

	_, Error_Name, Decode_Name := Genfilename(filepath)

	seqs := ReadFasta(Error_Name)

	set.Init(len(seqs))

	if params.Option != Gungnir_Trit_Params {
		if threads_num2 == 1 {
			deco, deco_suc = Decode_Parallel(seqs, Hmax, threads_num1, EDmax, set, params)
		} else if threads_num1 == 1 {
			deco, deco_suc = Decode_Multithread(seqs, Hmax, threads_num2, EDmax, set, params)
		} else {
			deco, deco_suc = Decode_Mix(seqs, Hmax, threads_num1, threads_num2, EDmax, set, params)
		}
	} else {
		if threads_num2 == 1 {
			deco, deco_suc = Decode_Three_Parallel(seqs, Hmax, threads_num1, EDmax, set, params)
		} else if threads_num1 == 1 {
			deco, deco_suc = Decode_Three_Multithread(seqs, Hmax, threads_num2, EDmax, set, params)
		} else {
			deco, deco_suc = Decode_Three_Mix(seqs, Hmax, threads_num1, threads_num2, EDmax, set, params)
		}
	}

	res := make([]Block, len(seqs))

	for i := 0; i < len(seqs); i++ {
		if deco_suc[i] && deco[i].BlockID < len(seqs) {
			if len(res[deco[i].BlockID].Hash) == 0 {
				res[deco[i].BlockID] = deco[i]
			}
		}

	}

	dec_seqs_temp := Encode(res, params)

	GenFasta(dec_seqs_temp, Decode_Name)
	SaveBoolsToFile(deco_suc, filepath+"/whetheroutput")

	fmt.Println("Decoded finished, Hmax:", Hmax)

	return AnalysisFailure(set, filepath, params)
}

func DealWithError(filepath string, threads_num1, threads_num2 int, Hmax int, EDmax int, params Params) (int, int) {
	_, Error_Name, Decode_Name := Genfilename(filepath)

	data, _ := Readjson()

	dec_seqs := ReadFasta(Decode_Name)
	err_seqs := ReadFasta(Error_Name)

	set := &IDtobeDecode{}

	AnalysisFailure(set, filepath, params)

	deco_suc, _ := LoadBoolsFromFile(filepath + "/whetheroutput")

	tobefix := make([]string, 0)
	tobefixID := make([]int, 0)
	fixed := 0
	for i := 0; i < len(err_seqs); i++ {
		if !deco_suc[i] {
			tobefix = append(tobefix, err_seqs[i])
			tobefixID = append(tobefixID, i)
		}
	}
	fmt.Println(len(tobefix), " sequences with no output! Hmax:", Hmax)

	var new_set []Block
	var new_decres []bool

	if params.Option != Gungnir_Trit_Params {
		if threads_num2 == 1 {
			new_set, new_decres = Decode_Parallel(tobefix, Hmax, threads_num1, EDmax, set, params)
		} else if threads_num1 == 1 {
			new_set, new_decres = Decode_Multithread(tobefix, Hmax, threads_num2, EDmax, set, params)
		} else {
			new_set, new_decres = Decode_Mix(tobefix, Hmax, threads_num1, threads_num2, EDmax, set, params)
		}
	} else {
		if threads_num2 == 1 {
			new_set, new_decres = Decode_Three_Parallel(tobefix, Hmax, threads_num1, EDmax, set, params)
		} else if threads_num1 == 1 {
			new_set, new_decres = Decode_Three_Multithread(tobefix, Hmax, threads_num2, EDmax, set, params)
		} else {
			new_set, new_decres = Decode_Three_Mix(tobefix, Hmax, threads_num1, threads_num2, EDmax, set, params)
		}
	}

	for i := 0; i < len(new_set); i++ {
		if new_decres[i] {
			index := new_set[i].BlockID
			if index >= len(dec_seqs) {
				continue
			}
			if params.Option != Gungnir_Trit_Params {
				deco_suc[tobefixID[i]] = true
				dec_seqs[index] = Block2DNA(new_set[i], data, params)
			} else {
				deco_suc[tobefixID[i]] = true
				dec_seqs[index] = Block2DNA_Three(new_set[i], params)
			}

			fixed += 1
		}
	}
	fmt.Println(len(tobefix), " no ouput, ", fixed, " fixed!")
	SaveBoolsToFile(deco_suc, filepath+"/whetheroutput")
	GenFasta(dec_seqs, Decode_Name)
	return AnalysisFailure(set, filepath, params)
}

func DecodeWithEDmax(filepath string, threads_num1, threads_num2 int, params Params) {

	var fail, old_fail int
	var totalseq int
	fisrt_decode := true

	Continue := true
	EDmax := 3

	for Continue {
		if fisrt_decode {
			fail, totalseq = FirstDecode(filepath, threads_num1, threads_num2, Maxhypo_firstround, EDmax, params)
			fisrt_decode = false
			old_fail = totalseq
			fmt.Println("First round Finished at Edit Distance upperbound: ", EDmax)
		} else {
			fail, _ = DealWithError(filepath, threads_num1, threads_num2, Maxhypo_firstround, EDmax, params)
			fmt.Println("First round Finished at Edit Distance upperbound: ", EDmax)
		}

		if fail < old_fail {
			if float64(fail)/float64(totalseq) < 0.5 {
				old_fail = fail
				fmt.Println("Second round begin!")
				fail, _ = DealWithError(filepath, threads_num1, threads_num2, Maxhypo_secondround, EDmax, params)
				fmt.Println("Second round Finished at Edit Distance upperbound: ", EDmax)
			}

			if float64(fail)/float64(totalseq) < 0.2 && fail < old_fail {
				old_fail = fail
				fmt.Println("Third round begin!")
				fail, _ = DealWithError(filepath, threads_num1, threads_num2, Maxhypo_thirdround, EDmax, params)
				fmt.Println("Third round Finished at Edit Distance upperbound: ", EDmax)
			}

			if float64(fail)/float64(totalseq) < 0.1 && fail < old_fail {
				old_fail = fail
				fmt.Println("Forth round begin!")
				fail, _ = DealWithError(filepath, threads_num1, threads_num2, Maxhypo_forthround, EDmax, params)
				fmt.Println("Forth round Finished at Edit Distance upperbound: ", EDmax)
			}

			if float64(fail)/float64(totalseq) < 0.01 && fail < old_fail {
				fmt.Println("Fifth round begin!")
				fail, _ = DealWithError(filepath, threads_num1, threads_num2, Maxhypo_ultra, EDmax, params)
				fmt.Println("Fifth round Finished at Edit Distance upperbound: ", EDmax)
			}

		} else {
			fmt.Println("Fail to decode at Edit Distance upperbound: ", EDmax)
		}

		EDmax += 1
		old_fail = fail

		if fail == 0 {
			Continue = false
		}

		if EDmax > 20 {
			Continue = false
		}

	}

	suc_rate, fail, mistake, total := AnalysisAll(filepath, params)
	precision := float64(total-fail-mistake) / float64(total-fail)
	recall := float64(total-fail) / float64(total)
	fmt.Println("Data recovery: ", suc_rate, " Precision: ", precision, " Recall: ", recall)
}

func DecodeWithFixEDmax(filepath string, threads_num1, threads_num2 int, EDmax int, params Params) {
	fmt.Println("First round begin!")
	FirstDecode(filepath, threads_num1, threads_num2, Maxhypo_firstround, EDmax, params)
	fmt.Println("First round finish!")

	fmt.Println("Second round begin!")
	DealWithError(filepath, threads_num1, threads_num2, Maxhypo_secondround, EDmax, params)
	fmt.Println("Second round finish!")

	fmt.Println("Third round begin!")
	DealWithError(filepath, threads_num1, threads_num2, Maxhypo_thirdround, EDmax, params)
	fmt.Println("Third round finish!")

	fmt.Println("Forth round begin!")
	DealWithError(filepath, threads_num1, threads_num2, Maxhypo_forthround, EDmax, params)
	fmt.Println("Forth round finish!")

	suc_rate, fail, mistake, total := AnalysisAll(filepath, params)
	precision := float64(total-fail-mistake) / float64(total-fail)
	recall := float64(total-fail) / float64(total)
	fmt.Println("Data recovery: ", suc_rate, " Precision: ", precision, " Recall: ", recall)
}

func ReconstructFile(filepath string, params Params) {
	_, _, Decode_Name := Genfilename(filepath)
	numCores := runtime.NumCPU() * 2 / 3
	set := &IDtobeDecode{}
	dec_seqs := ReadFasta(Decode_Name)
	set.Init(len(dec_seqs))
	var datablock []Block
	if params.Option != Gungnir_Trit_Params {
		datablock, _ = Decode_Parallel(dec_seqs, Maxhypo_simple, numCores, 100, set, params)
	} else {
		datablock, _ = Decode_Three_Parallel(dec_seqs, Maxhypo_simple, numCores, 100, set, params)
	}

	bits := BlocksintoBits(datablock, params)

	file := Bits2File(bits)
	WriteStringToFile(file, filepath+"/output")
}
