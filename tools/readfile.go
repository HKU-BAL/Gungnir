package tools

import (
	"bufio"
	"compress/gzip"
	"encoding/json"
	"fmt"
	"io"
	"math/big"
	"os"
	"sort"
	"strconv"
	"strings"
)

type Kmer_json struct {
	X     int `json:"X"`
	I     int `json:"I"`
	D     int `json:"D"`
	Total int `json:"Total"`
}

type Kmer struct {
	X_rate float64
	I_rate float64
	D_rate float64
	E_rate float64
}

type Kmer_rank struct {
	Input  string
	E_rate float64
}

func ReadjsonAll() (map[string]Kmer, []map[string]Kmer, float64) {
	fulljson, errorjson := Readjson()
	JsonSet := make([]map[string]Kmer, KmerSize)
	for i := 0; i < KmerSize; i++ {
		JsonSet[i], _ = ReadjsonIndex(i)
	}

	return fulljson, JsonSet, errorjson / KmerSize
}

func Readjson() (map[string]Kmer, float64) {
	res, errorrate := ReadjsonIndex(0)
	for i := 1; i < KmerSize; i++ {
		etable, etemp := ReadjsonIndex(i)
		errorrate += etemp
		for key, value := range etable {
			temp := res[key]
			temp.D_rate += value.D_rate
			temp.E_rate += value.E_rate
			temp.I_rate += value.I_rate
			temp.X_rate += value.X_rate
			res[key] = temp
		}
	}
	return res, errorrate
}

func ReadjsonIndex(Index int) (map[string]Kmer, float64) {
	filepath := "../error_pattern/hg002_ONT_" + strconv.Itoa(KmerSize) + "_" + strconv.Itoa(Index) + ".json"
	file, _ := os.Open(filepath)
	defer file.Close()
	byteValue, _ := io.ReadAll(file)
	var data map[string]Kmer_json
	json.Unmarshal(byteValue, &data)

	Xsum := big.NewInt(0)
	Isum := big.NewInt(0)
	Dsum := big.NewInt(0)
	Totalsum := big.NewInt(0)

	res := make(map[string]Kmer)
	for key, value := range data {
		var temp Kmer
		temp.X_rate = float64(value.X) / float64(value.Total)
		temp.I_rate = float64(value.I) / float64(value.Total)
		temp.D_rate = float64(value.D) / float64(value.Total)
		temp.E_rate = temp.X_rate + temp.I_rate + temp.D_rate
		res[key] = temp

		Xsum.Add(Xsum, big.NewInt(int64(value.X)))
		Isum.Add(Isum, big.NewInt(int64(value.I)))
		Dsum.Add(Dsum, big.NewInt(int64(value.D)))
		Totalsum.Add(Totalsum, big.NewInt(int64(value.Total)))
	}

	Xfloat := new(big.Float).SetInt(Xsum)
	Ifloat := new(big.Float).SetInt(Isum)
	Dfloat := new(big.Float).SetInt(Dsum)
	TotalFloat := new(big.Float).SetInt(Totalsum)

	resultX := new(big.Float).Quo(Xfloat, TotalFloat)
	resultI := new(big.Float).Quo(Ifloat, TotalFloat)
	resultD := new(big.Float).Quo(Dfloat, TotalFloat)
	xrate, _ := resultX.Float64()
	irate, _ := resultI.Float64()
	drate, _ := resultD.Float64()

	error_rate := xrate + irate + drate

	return res, error_rate
}

func RankKmer(data map[string]Kmer) []Kmer_rank {
	var res []Kmer_rank
	for key, value := range data {
		var temp Kmer_rank
		temp.Input = key
		temp.E_rate = value.E_rate
		res = append(res, temp)
	}
	sort.Slice(res, func(i, j int) bool {
		return res[i].E_rate < res[j].E_rate
	})
	return res
}

func Readfile(filepath string) string {
	data, _ := os.ReadFile(filepath)
	content := string(data)
	return content
}

func WriteStringToFile(content string, filepath string) error {
	file, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0644)
	if err != nil {
		return fmt.Errorf("Can not open file: %v", err)
	}
	defer file.Close()

	_, err = file.WriteString(content)
	if err != nil {
		return fmt.Errorf("Fail to write: %v", err)
	}

	return nil
}

func File2bytes(content string) []byte {
	bytes := []byte(content)
	if len(bytes) == 0 {
		return nil
	}
	return bytes
}

func Bytes2File(bytes []byte) string {
	content := string(bytes)
	return content
}

func Bytes2Bits(bytes []byte) []int {
	bits := make([]int, 8*len(bytes))
	for i := 0; i < len(bytes); i++ {
		val := int(bytes[i])
		for j := 0; j < 8; j++ {
			bits[8*i+j] = val & 1
			val = val >> 1
		}
	}
	return bits
}

func Bits2Bytes(bits []int) []byte {
	bytesLen := (len(bits)-1)/8 + 1
	bytes := make([]byte, bytesLen)
	for i := 0; i < bytesLen; i++ {
		val := 0
		for j := 0; j < 8 && 8*i+j < len(bits); j++ {
			val += (bits[8*i+j] << j)
		}
		bytes[i] = byte(val)
	}

	return bytes
}

func File2Bits(content string) []int {
	bytes := File2bytes(content)
	return Bytes2Bits(bytes)
}

func Bits2File(data []int) string {
	bytes := Bits2Bytes(data)
	return Bytes2File(bytes)
}

// 0-less 1-more
func Rank(symbol int, length int) {
	data, _ := Readjson()

	rank := RankKmer(data)
	for i := 0; i < length; i++ {
		if symbol == 0 {
			fmt.Println(rank[i].Input, rank[i].E_rate)
		} else {
			fmt.Println(rank[len(rank)-i-1].Input, rank[len(rank)-i-1].E_rate)
		}
	}
}

func ReadFastqgz(filepath string) ([]string, [][]int) {
	file, _ := os.Open(filepath)
	defer file.Close()
	reader, _ := gzip.NewReader(file)
	defer reader.Close()
	scanner := bufio.NewScanner(reader)
	linenum := 0
	sequence := make([]string, 0)
	qual := make([][]int, 0)
	for scanner.Scan() {
		if linenum&3 == 1 {
			sequence = append(sequence, scanner.Text())
		}
		if linenum&3 == 3 {
			temp := scanner.Text()
			qual_temp := make([]int, len(temp))
			for i := 0; i < len(temp); i++ {
				qual_temp[i] = int(temp[i]) - 33
			}
			qual = append(qual, qual_temp)
		}
		if linenum > 10 {
			break
		}
		linenum += 1
	}

	return sequence, qual

}

func GenFasta(sequence []string, filepath string) {
	f, _ := os.Create(filepath)
	w := bufio.NewWriter(f)

	for i := 0; i < len(sequence); i++ {
		w.WriteString(">index_" + strconv.Itoa(i) + "\n")
		w.WriteString(sequence[i] + "\n")
	}

	w.Flush()
	f.Close()
}

func ReadFasta(filepath string) []string {
	file, _ := os.Open(filepath)
	defer file.Close()

	scanner := bufio.NewScanner(file)

	res := make([]string, 0)

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			continue
		} else {
			if len(line) > 50 {
				res = append(res, line)
			}

		}
	}

	return res
}

func SaveBoolsToFile(deco_suc []bool, filepath string) error {
	file, err := os.Create(filepath)
	if err != nil {
		return fmt.Errorf("Fail to create a file: %w", err)
	}
	defer file.Close()

	writer := bufio.NewWriter(file)
	for _, val := range deco_suc {
		line := "false"
		if val {
			line = "true"
		}
		_, err := writer.WriteString(line + "\n")
		if err != nil {
			return fmt.Errorf("Fail to write a file: %w", err)
		}
	}
	return writer.Flush()
}

func LoadBoolsFromFile(filepath string) ([]bool, error) {
	file, err := os.Open(filepath)
	if err != nil {
		return nil, fmt.Errorf("Fail to open a file: %w", err)
	}
	defer file.Close()

	var deco_suc []bool
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		switch line {
		case "true":
			deco_suc = append(deco_suc, true)
		case "false":
			deco_suc = append(deco_suc, false)
		default:
			return nil, fmt.Errorf("Invalid bool value: %s", line)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("Fail to read a file: %w", err)
	}
	return deco_suc, nil
}
