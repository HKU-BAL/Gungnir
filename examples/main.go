package main

import (
	"dnastorage/tools"
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"strconv"
	"time"
)

const coderate = 0.5

var PayloadLen = int(math.Ceil(100 * coderate))
var HashLen = int(100 - PayloadLen)

var Subrate, Delrate, Insrate float64

var ErrorRate = 0.2

var ErrorNum = 3

func UpdateError() {
	Subrate = ErrorRate / 3
	Insrate = ErrorRate / 3
	Delrate = ErrorRate / 3
}

func Genfilename(params tools.Params) (string, string, string) {
	var Origin = "Time_Origin_"
	var Error = "Time_Add_Error_"
	var Decode = "Time_Decoded_"
	var option string
	if params.Option == tools.Gungnir_Default_Params {
		option = "Default_"
	} else if params.Option == tools.Gungnir_ONT_Params {
		option = "ONT_"
	} else {
		option = "Three_"
	}
	errorrate := ErrorRate

	crate := fmt.Sprintf("%.2f", coderate) + "_"

	sub := fmt.Sprintf("%.4f", Subrate) + "_"
	del := fmt.Sprintf("%.4f", Delrate) + "_"
	ins := fmt.Sprintf("%.4f", Insrate) + "_"
	total := fmt.Sprintf("%.4f", errorrate)
	Origin = Origin + option + crate + sub + del + ins + total
	Error = Error + option + crate + sub + del + ins + total
	Decode = Decode + option + crate + sub + del + ins + total
	return Origin, Error, Decode
}

func GC(params tools.Params) {
	hamlet := tools.Readfile("../file/The Ugly Duckling")
	fmt.Println("Length of txt (23 per block)", len(hamlet))

	b := tools.File2Bits(hamlet)

	d := tools.BitsintoBlocks(b, params)

	e := tools.Encode(d, params)
	fmt.Println("GC", tools.Calgc(e))
	tools.Calexpecterror(e)
}

func Gen_Test(params tools.Params) {
	hamlet := tools.Readfile("../file/The Ugly Duckling")
	b := tools.File2Bits(hamlet)

	d := tools.BitsintoBlocks(b, params)
	fmt.Println("Block Num: ", len(d))

	set := &tools.IDtobeDecode{}
	set.Init(len(d))

	ori_seqs := tools.Encode(d, params)

	temp_seqs := make([]string, len(ori_seqs))
	for i := 0; i < len(temp_seqs); i++ {
		temp_seqs[i] = tools.GenRevString(ori_seqs[i])
	}

	tools.Calexpecterror(ori_seqs)
	tools.Calexpecterror(temp_seqs)

	gc := tools.Calgc(ori_seqs)

	count := 0
	for i := 0; i < len(gc); i++ {
		if gc[i] < 0.4 || gc[i] > 0.6 {
			fmt.Println(i, gc[i])
			count += 1
		}
	}

	fmt.Println("GC mismatch num: ", count)

	sum := 0.0
	for i := 0; i < len(gc); i++ {
		sum += gc[i]
	}
	sum = sum / float64(len(gc))

	fmt.Println(sum)

}

func Gen(set *tools.IDtobeDecode, params tools.Params) {

	Origin_Name, Error_Name, _ := Genfilename(params)

	hamlet := tools.Readfile("../file/The Ugly Duckling")

	b := tools.File2Bits(hamlet)

	d := tools.BitsintoBlocks(b, params)

	fmt.Println("Block Num: ", len(d))
	set.Init(len(d))

	ori_seqs := tools.Encode(d, params)

	tools.GenFasta(ori_seqs, Origin_Name)

	seqs := tools.AddError(ori_seqs, Subrate, Delrate, Insrate, ErrorRate) // Fixed Error Rate

	tools.GenFasta(seqs, Error_Name)

}

func FirstDecode(set *tools.IDtobeDecode, threads_num int, params tools.Params) (float64, int, int, int) {
	var deco []tools.Block
	var deco_suc []bool

	_, Error_Name, Decode_Name := Genfilename(params)

	seqs := tools.ReadFasta(Error_Name)

	if params.Option != 2 {
		deco, deco_suc = tools.Decode_Parallel(seqs, tools.Maxhypo_firstround, threads_num, set, params) // fisrt round
	} else {
		deco, deco_suc = tools.Decode_Three_Parallel(seqs, tools.Maxhypo_firstround, threads_num, set, params) // fisrt round
	}

	res := make([]tools.Block, len(seqs))

	for i := 0; i < len(seqs); i++ {
		if deco_suc[i] && deco[i].BlockID < len(seqs) {
			if len(res[deco[i].BlockID].Hash) == 0 {
				res[deco[i].BlockID] = deco[i]
			}
		}

	}

	dec_seqs_temp := tools.Encode(res, params) // fisrt round

	tools.GenFasta(dec_seqs_temp, Decode_Name) // fisrt round

	fmt.Println("First Round finished!")

	return Analysis(set, params)
}

func Analysis(set *tools.IDtobeDecode, params tools.Params) (float64, int, int, int) {
	if params.Option != tools.Gungnir_Three_Params {
		Origin_Name, Error_Name, Decode_Name := Genfilename(params)
		seqs := tools.ReadFasta(Origin_Name)
		dec_seqs := tools.ReadFasta(Decode_Name)
		err_seqs := tools.ReadFasta(Error_Name)

		data, _ := tools.Readjson()
		var zero_block tools.Block
		zero_seq := tools.Block2DNA(zero_block, data, params)

		// for i := 0; i < len(seqs); i++ {
		// 	dec_seqs[i] = seqs[i]
		// }
		// dec_seqs[919] = zero_seq

		// tools.GenFasta(dec_seqs, Decode_Name)

		res := make([]int, 0)
		count := 0
		fail := 0
		for i := 0; i < len(seqs); i++ {
			val := tools.Hamming_Distrance([]rune(dec_seqs[i]), []rune(seqs[i]))
			res = append(res, val)
		}
		fmt.Println("Total Blocks", len(res))

		tempset := make([]int, 0)
		for i := 0; i < len(res); i++ {
			if res[i] != 0 {
				if dec_seqs[i] == zero_seq {
					fmt.Println(tools.EditingDistance(seqs[i], err_seqs[i]))
					fmt.Println(seqs[i])
					fmt.Println(err_seqs[i])

					tempset = append(tempset, i)
					fail++
				}
				count++
			}
		}
		fmt.Println("Error Num", count-fail)
		fmt.Println("Failure Num", fail)
		set.InitWithtempset(len(seqs), tempset)
		suc_rate := 1.0 - (float64(count) / float64(len(res)))
		return suc_rate, fail, count - fail, len(res)
	} else {
		Origin_Name, _, Decode_Name := Genfilename(params)
		seqs := tools.ReadFasta(Origin_Name)
		dec_seqs := tools.ReadFasta(Decode_Name)
		// err_seqs := tools.ReadFasta(Error_Name)
		res := make([]int, 0)
		count := 0
		fail := 0
		for i := 0; i < len(seqs); i++ {
			val := tools.Hamming_Distrance([]rune(dec_seqs[i]), []rune(seqs[i]))
			res = append(res, val)
		}
		fmt.Println("Total Blocks", len(res))
		var zero_block tools.Block
		zero_seq := tools.Block2DNA_Three(zero_block, params)
		tempset := make([]int, 0)
		for i := 0; i < len(res); i++ {
			if res[i] != 0 {
				if dec_seqs[i] == zero_seq {
					tempset = append(tempset, i)
					fail++
				}
				count++
			}
		}
		fmt.Println("Error Num", count-fail)
		fmt.Println("Failure Num", fail)
		set.InitWithtempset(len(seqs), tempset)
		suc_rate := 1.0 - (float64(count) / float64(len(res)))
		return suc_rate, fail, count - fail, len(res)
	}

}

// func TestDecodeError(threads_num int, set *tools.IDtobeDecode, params tools.Params) {
// 	_, _, Decode_Name := Genfilename(params)

// 	// data, _ := tools.Readjson()

// 	dec_seqs := tools.ReadFasta(Decode_Name)
// 	dec_test, test_res := tools.Decode_Parallel_Raw(dec_seqs, threads_num, len(dec_seqs), params) // test error
// 	errornum := 0
// 	for i := 0; i < len(dec_test); i++ {
// 		if !test_res[i] {
// 			// fmt.Println(i, " is Error!")
// 			errornum++
// 		}
// 	}
// 	fmt.Println(errornum, " error has been found!")
// }

func DealWithError_Parallel(hypoMax int, threads_num int, set *tools.IDtobeDecode, params tools.Params) (float64, int, int, int) {
	_, Error_Name, Decode_Name := Genfilename(params)

	data, _ := tools.Readjson()

	dec_seqs := tools.ReadFasta(Decode_Name)
	err_seqs := tools.ReadFasta(Error_Name)

	var dec_test []tools.Block
	var test_res []bool

	if params.Option != tools.Gungnir_Three_Params {
		dec_test, test_res = tools.Decode_Parallel_Raw(dec_seqs, threads_num, len(dec_seqs), params) // test error
	} else {
		dec_test, test_res = tools.Decode_Parallel_Three_Raw(dec_seqs, threads_num, len(dec_seqs), params) // test error
	}

	tobefix := make([]string, 0)
	errornum := 0
	fixed := 0
	for i := 0; i < len(dec_test); i++ {
		if !test_res[i] {
			errornum += 1
			tobefix = append(tobefix, err_seqs[i])
		}
	}

	fmt.Println(errornum, " error has been found!")

	var new_set []tools.Block
	var new_decres []bool

	if params.Option != tools.Gungnir_Three_Params {
		new_set, new_decres = tools.Decode_Parallel(tobefix, hypoMax, threads_num, set, params)
	} else {
		new_set, new_decres = tools.Decode_Three_Parallel(tobefix, hypoMax, threads_num, set, params)
	}

	for i := 0; i < len(new_set); i++ {
		if new_decres[i] {
			index := new_set[i].BlockID
			if index >= len(dec_test) {
				continue
			}
			if params.Option != tools.Gungnir_Three_Params {
				dec_seqs[index] = tools.Block2DNA(new_set[i], data, params)
			} else {
				dec_seqs[index] = tools.Block2DNA_Three(new_set[i], params)
			}

			fixed += 1
		}
	}
	fmt.Println(errornum, " error, ", fixed, " fixed!")
	tools.GenFasta(dec_seqs, Decode_Name)
	return Analysis(set, params)
}

func DealWithError_Mix(hypoMax int, threads_num1 int, threads_num2 int, set *tools.IDtobeDecode, params tools.Params) (float64, int, int, int) {
	_, Error_Name, Decode_Name := Genfilename(params)

	data, _ := tools.Readjson()

	dec_seqs := tools.ReadFasta(Decode_Name)
	err_seqs := tools.ReadFasta(Error_Name)

	var dec_test []tools.Block
	var test_res []bool

	if params.Option != tools.Gungnir_Three_Params {
		dec_test, test_res = tools.Decode_Parallel_Raw(dec_seqs, threads_num1, len(dec_seqs), params) // test error
	} else {
		dec_test, test_res = tools.Decode_Parallel_Three_Raw(dec_seqs, threads_num1, len(dec_seqs), params) // test error
	}

	tobefix := make([]string, 0)
	errornum := 0
	fixed := 0
	for i := 0; i < len(dec_test); i++ {
		if !test_res[i] {
			errornum += 1
			tobefix = append(tobefix, err_seqs[i])
		}
	}

	fmt.Println(errornum, " error has been found!")

	var new_set []tools.Block
	var new_decres []bool

	if params.Option != tools.Gungnir_Three_Params {
		new_set, new_decres = tools.Decode_Mix(tobefix, hypoMax, threads_num1, threads_num2, set, params)
	} else {
		new_set, new_decres = tools.Decode_Three_Mix(tobefix, hypoMax, threads_num1, threads_num2, set, params)
	}

	for i := 0; i < len(new_set); i++ {
		if new_decres[i] {
			index := new_set[i].BlockID
			if index >= len(dec_test) {
				continue
			}
			if params.Option != tools.Gungnir_Three_Params {
				dec_seqs[index] = tools.Block2DNA(new_set[i], data, params)
			} else {
				dec_seqs[index] = tools.Block2DNA_Three(new_set[i], params)
			}

			fixed += 1
		}
	}
	fmt.Println(errornum, " error, ", fixed, " fixed!")
	tools.GenFasta(dec_seqs, Decode_Name)
	return Analysis(set, params)
}

func DealWithError(hypoMax int, threads_num int, set *tools.IDtobeDecode, params tools.Params) (float64, int, int, int) {
	_, Error_Name, Decode_Name := Genfilename(params)

	data, _ := tools.Readjson()

	dec_seqs := tools.ReadFasta(Decode_Name)
	err_seqs := tools.ReadFasta(Error_Name)

	var dec_test []tools.Block
	var test_res []bool

	if params.Option != tools.Gungnir_Three_Params {
		dec_test, test_res = tools.Decode_Parallel_Raw(dec_seqs, threads_num, len(dec_seqs), params) // test error
	} else {
		dec_test, test_res = tools.Decode_Parallel_Three_Raw(dec_seqs, threads_num, len(dec_seqs), params) // test error
	}

	errornum := 0
	fixed := 0
	for i := 0; i < len(dec_test); i++ {
		if !test_res[i] {
			errornum += 1
		}
	}

	fmt.Println(errornum, " error has been found!")

	for i := 0; i < len(dec_test); i++ {
		if !test_res[i] {
			fmt.Println(strconv.Itoa(i) + " error! Try to fix...")
			if params.Option != tools.Gungnir_Three_Params {
				dec_test[i], test_res[i] = tools.Updateinfo_Multithread([]rune(err_seqs[i]), hypoMax, data, threads_num, set, params)
			} else {
				dec_test[i], test_res[i] = tools.Updateinfo_Three_Multithread([]rune(err_seqs[i]), hypoMax, threads_num, set, params)
			}

			if test_res[i] {
				fmt.Println(strconv.Itoa(i) + " has been fixed!")
				if params.Option != tools.Gungnir_Three_Params {
					dec_seqs[i] = tools.Block2DNA(dec_test[i], data, params)
				} else {
					dec_seqs[i] = tools.Block2DNA_Three(dec_test[i], params)
				}

				fixed += 1
				tools.GenFasta(dec_seqs, Decode_Name)
			} else {
				fmt.Println(strconv.Itoa(i) + " fail to fix!")
			}
		}
	}
	fmt.Println(errornum, " error, ", fixed, " fixed!")
	return Analysis(set, params)
}

func numinset(a int, set []int) bool {
	for i := 0; i < len(set); i++ {
		if a == set[i] {
			return true
		}
	}
	return false
}

func DealWithError_MultiServer(hypoMax int, threads_num int, set *tools.IDtobeDecode, params tools.Params) (float64, int, int, int) {
	_, Error_Name, Decode_Name := Genfilename(params)

	data, _ := tools.Readjson()

	dec_seqs := tools.ReadFasta(Decode_Name)
	err_seqs := tools.ReadFasta(Error_Name)

	tried_set := make([]int, 0)

	for {
		var dec_test []tools.Block
		var test_res []bool

		if params.Option != tools.Gungnir_Three_Params {
			dec_test, test_res = tools.Decode_Parallel_Raw(dec_seqs, threads_num, len(dec_seqs), params) // test error
		} else {
			dec_test, test_res = tools.Decode_Parallel_Three_Raw(dec_seqs, threads_num, len(dec_seqs), params) // test error
		}

		set_tobefix := make([]int, 0)
		for i := 0; i < len(dec_test); i++ {
			if !test_res[i] && !numinset(i, tried_set) {
				set_tobefix = append(set_tobefix, i)
			}
		}

		if len(set_tobefix) == 0 {
			break
		}

		r_seed := rand.New(rand.NewSource(time.Now().UnixNano()))
		rng := r_seed.Intn(len(set_tobefix))
		dec_id := set_tobefix[rng]
		tried_set = append(tried_set, dec_id)

		fmt.Println(strconv.Itoa(dec_id) + " failure! Try to fix...")
		if params.Option != tools.Gungnir_Three_Params {
			dec_test[dec_id], test_res[dec_id] = tools.Updateinfo_Multithread([]rune(err_seqs[dec_id]), hypoMax, data, threads_num, set, params)
		} else {
			dec_test[dec_id], test_res[dec_id] = tools.Updateinfo_Three_Multithread([]rune(err_seqs[dec_id]), hypoMax, threads_num, set, params)
		}

		if test_res[dec_id] {
			fmt.Println(strconv.Itoa(dec_id) + " has been fixed!")
			dec_seqs = tools.ReadFasta(Decode_Name)
			if params.Option != tools.Gungnir_Three_Params {
				dec_seqs[dec_id] = tools.Block2DNA(dec_test[dec_id], data, params)
			} else {
				dec_seqs[dec_id] = tools.Block2DNA_Three(dec_test[dec_id], params)
			}

			tools.GenFasta(dec_seqs, Decode_Name)
		} else {
			fmt.Println(strconv.Itoa(dec_id) + " fail to fix!")
		}

		Analysis(set, params)

	}

	return Analysis(set, params)
}

func Pipeline(params tools.Params) float64 {
	set := &tools.IDtobeDecode{}
	now := time.Now()
	numCores := runtime.NumCPU() * 2 / 3

	Gen(set, params)
	suc := 0.0
	var fail int
	old_fail := set.MaxID
	fisrt_decode := true

	Continue := true

	for Continue {
		if fisrt_decode {
			suc, fail, _, _ = FirstDecode(set, numCores, params)
			fmt.Printf("First0 round takes (%s)\n", time.Since(now))
			fisrt_decode = false
			now = time.Now()
		} else {
			now = time.Now()
			suc, fail, _, _ = DealWithError_Parallel(tools.Maxhypo_firstround, numCores, set, params)
			fmt.Printf("First0 round takes (%s)\n", time.Since(now))
			now = time.Now()
		}

		if fail < old_fail {
			fmt.Println("Try to continue at Edit Distance upperbound: ", params.MaxError)
			// old_fail = fail
			// fmt.Println("Try to continue at Edit Distance upperbound: ", params.MaxError)
			// _, fail, _, _ = DealWithError_Parallel(tools.Maxhypo_firstround, numCores, set, params)
			// fmt.Printf("First1 round takes (%s)\n", time.Since(now))

			// now = time.Now()
			// if fail < old_fail {
			// 	old_fail = fail
			// 	fmt.Println("First2 round begin!")
			// 	_, fail, _, _ = DealWithError_Parallel(tools.Maxhypo_firstround, numCores, set, params)
			// 	fmt.Printf("First2 round takes (%s)\n", time.Since(now))
			// 	now = time.Now()
			// }

			// if fail < old_fail {
			// 	fmt.Println("First3 round begin!")
			// 	suc, fail, _, _ = DealWithError_Parallel(tools.Maxhypo_firstround, numCores, set, params)
			// 	fmt.Printf("First3 round takes (%s)\n", time.Since(now))
			// }

			// old_fail = fail

			if float64(fail)/float64(set.MaxID) < 0.5 {
				now = time.Now()
				old_fail = fail
				fmt.Println("Second1 round begin!")
				suc, fail, _, _ = DealWithError_Parallel(tools.Maxhypo_secondround, numCores, set, params)
				fmt.Printf("Second1 round takes (%s)\n", time.Since(now))

				// if fail < old_fail {
				// 	now = time.Now()
				// 	fmt.Println("Second2 round begin!")
				// 	suc, fail, _, _ = DealWithError_Parallel(tools.Maxhypo_secondround, numCores, set, params)
				// 	fmt.Printf("Second2 round takes (%s)\n", time.Since(now))
				// 	old_fail = fail
				// }
			}

			if float64(fail)/float64(set.MaxID) < 0.2 && fail < old_fail {
				now = time.Now()
				old_fail = fail

				fmt.Println("Third round begin!")
				// thread_num := numCores
				// if thread_num > 64 {
				// 	thread_num = 64
				// }
				suc, fail, _, _ = DealWithError_Mix(tools.Maxhypo_thirdround, 22, 4, set, params)
				// DealWithError_Parallel(tools.Maxhypo_thirdround, numCores, set, params)
				fmt.Printf("Third round takes (%s)\n", time.Since(now))
			}

			if float64(fail)/float64(set.MaxID) < 0.1 && fail < old_fail {
				now = time.Now()
				old_fail = fail
				fmt.Println("Forth round begin!")

				thread1 := numCores / 8
				thread2 := 8

				if thread1 < 8 {
					thread1 = 8
					thread2 = numCores / 8
				}

				suc, fail, _, _ = DealWithError_Mix(tools.Maxhypo_forthround, thread1, thread2, set, params)
				fmt.Printf("Forth round takes (%s)\n", time.Since(now))
			}

			if float64(fail)/float64(set.MaxID) < 0.01 && fail < old_fail {
				now = time.Now()
				suc, fail, _, _ = DealWithError(tools.Maxhypo_ultra, numCores, set, params)
				fmt.Printf("Fifth round takes (%s)\n", time.Since(now))
			}

		} else {
			fmt.Println("Fail to decode at Edit Distance upperbound: ", params.MaxError)
		}

		params.MaxError += 1
		params.Penalty_upperbound += 1
		old_fail = fail

		if fail == 0 {
			Continue = false
		}

		if params.MaxError > 20 {
			Continue = false
		}

	}

	return suc
}

func NewRound_Multiserver(Maxhypo int, params tools.Params) {
	set := &tools.IDtobeDecode{}
	now := time.Now()
	numCores := runtime.NumCPU() * 2 / 3
	Analysis(set, params)
	DealWithError_MultiServer(Maxhypo, numCores, set, params)
	fmt.Printf("Fifth round takes (%s)\n", time.Since(now))

}

func main() {
	paramsRaw := tools.GenParamsRaw(HashLen, PayloadLen)

	UpdateError()
	params := paramsRaw.Compile(tools.Gungnir_Default_Params, ErrorNum)
	params.MaxError = 20
	params.Penalty_upperbound = 21
	set := &tools.IDtobeDecode{}
	Analysis(set, params)
	DealWithError(200000000, 90, set, params)
	DealWithError(300000000, 90, set, params)
	// suc_rate := Pipeline(params)
	// fmt.Printf("Default setting! CodeRate: %.2f; Payload: %d; Hash: %d; SubRate: %.5f; InsRate: %.5f; DelRate: %.5f; ErrorRate: %.5f; Successful Rate: %.5f\n",
	// 	coderate, PayloadLen, HashLen, Subrate, Insrate, Delrate, ErrorRate, suc_rate)

}
