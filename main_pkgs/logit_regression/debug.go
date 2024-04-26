package main

import (
	mathFunc "lattigov4_dev/src_pkgs/mathematical_func"

	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// 打印密文的调试信息
func ciphertextDebug(params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []float64, decryptor rlwe.Decryptor, encoder ckks.Encoder) (valuesTest []float64) {
	tmp := encoder.Decode(decryptor.DecryptNew(ciphertext), params.LogSlots())

	valuesTest = make([]float64, len(tmp))
	for i := range tmp {
		valuesTest[i] = real(tmp[i])
	}
	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])
	fmt.Println()
	//precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, params.LogSlots(), 0)
	//fmt.Println(precStats.String())
	return
}

// 对数据进行预测并计算准确率
func predictDebug(data [][]float64, weightsWant []float64, weightTest []float64, doPrint bool) (correctWant int, correctTest int) {
	correctWant, correctTest = 0, 0
	for i := 0; i < len(data); i++ {
		sumWant, sumTest := weightsWant[0], weightTest[0]
		for j := 1; j < len(data[i]); j++ {
			sumWant += data[i][j] * weightsWant[j]
			sumTest += data[i][j] * weightTest[j]
		}
		sumWant, sumTest = mathFunc.Sigmoid(sumWant), mathFunc.Sigmoid(sumTest)
		if (sumWant > 0.5) == (data[i][0] > 0) {
			correctWant++
		}
		if (sumTest > 0.5) == (data[i][0] > 0) {
			correctTest++
		}
		if doPrint {
			fmt.Println("PredictTest:", (sumTest > 0.5), "\tReal:", (data[i][0] > 0), "\tPredictWant:", (sumWant > 0.5), "\tData:", data[i][0:10])
		}
	}
	return correctWant, correctTest
}
