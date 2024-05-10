package main

import (
	dataProc "lattigov4_dev/src_pkgs/data_process"
	mathFunc "lattigov4_dev/src_pkgs/mathematical_func"
	typeConv "lattigov4_dev/src_pkgs/type_convert"

	"fmt"
	"math/rand"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func main() {
	// 数据预处理
	rawData := dataProc.ReadAllFromCSV("lbw.txt")
	normData := dataProc.NormalizeMaxMin(rawData)
	relabeledData := dataProc.ModifyLabel(normData, len(normData[0])-1, 0)
	testData, trainData := dataProc.SplitData(relabeledData, 0.2)
	invertedData := dataProc.InvertFeatures(relabeledData[len(testData):], 0)
	paddedData := dataProc.Padding(invertedData, 4, 10)

	// 生成一个权重向量，并随机初始化(范围在0-1之间)
	weights := make([]float64, len(paddedData[0]))
	for i := 0; i < len(weights); i++ {
		weights[i] = rand.Float64()
	}

	// ckks密钥对象生成
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: 14,
		Q: []uint64{0x10000140001, 0x1fff90001, 0x200080001,
			0x1fff60001, 0x200100001, 0x1fff00001,
			0x1ffef0001, 0x1ffe60001, 0x2001d0001,
			0x2002e0001}, // 40 + 9 x 33

		P:            []uint64{0x1ffffe0001, 0x1ffffc0001}, // 37, 37
		RingType:     ring.ConjugateInvariant,
		LogSlots:     14,
		DefaultScale: 1 << 33,
	})
	if err != nil {
		panic(err)
	}
	encoder := ckks.NewEncoder(params)
	kgen := ckks.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPair()
	rlk := kgen.GenRelinearizationKey(sk, 1)
	rtk := kgen.GenRotationKeysForRotations([]int{1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, -1, -2, -4, -8}, false, sk)
	decryptor := ckks.NewDecryptor(params, sk)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: rtk})

	// 编码与加密
	dataCiphertext, _ := dataProc.EncMatrixF64(paddedData, params.MaxLevel(), params.DefaultScale(), params.LogSlots(), params, pk)
	//printDebug(params, dataCiphertext, paddedData[0], decryptor, encoder)
	weightsCiphertext, _ := dataProc.EncRowVectorF64(weights, 1024, params.MaxLevel(), params.DefaultScale(), params.LogSlots(), params, pk)
	//printDebug(params, weightsCiphertext, weights, decryptor, encoder)
	_, firstColPlainText := dataProc.EncRowVectorF64([]float64{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 1024, 8, params.DefaultScale(), params.LogSlots(), params, nil)

	// sigmoid函数近似
	approxSigmoidReversed := ckks.Approximate(mathFunc.SigmoidReversed, -8.0, 8.0, 14)
	// approxSigmoid := ckks.NewPoly([]complex128{0.5, -1.53048, 2.3533056, -1.3511295})

	// 迭代7次
	for t := 0; t < 7; t++ {
		// 样本与权重相乘并求和
		// 密文===============================================================================================================
		operateCiphertext := evaluator.MulRelinNew(dataCiphertext, weightsCiphertext)
		evaluator.Rescale(operateCiphertext, params.DefaultScale(), operateCiphertext)
		for i, r := 0, 1; i < 4; i, r = i+1, r<<1 {
			evaluator.Add(operateCiphertext, evaluator.RotateNew(operateCiphertext, r), operateCiphertext)
		}
		evaluator.Mul(operateCiphertext, firstColPlainText, operateCiphertext)
		evaluator.Rescale(operateCiphertext, params.DefaultScale(), operateCiphertext)
		for i, r := 0, 1; i < 4; i, r = i+1, r<<1 {
			evaluator.Add(operateCiphertext, evaluator.RotateNew(operateCiphertext, -r), operateCiphertext)
		}
		// 明文===============================================================================================================
		operateData := make([][]float64, len(paddedData))
		for i := 0; i < len(paddedData); i++ {
			operateData[i] = make([]float64, len(paddedData[i]))
			copy(operateData[i], paddedData[i])
		}
		for i := 0; i < len(operateData); i++ {
			sum := 0.0
			for j := 0; j < len(operateData[i]); j++ {
				sum += operateData[i][j] * weights[j]
			}
			for j := 0; j < len(operateData[i]); j++ {
				operateData[i][j] = sum
			}
		}
		//printDebug(params, dataC, operateData[0], decryptor, encoder)

		// 对现在样本每个值计算sigmoid函数值并乘上原始样本
		// 密文===============================================================================================================
		evaluator.MultByConst(operateCiphertext, 1.0/8.0, operateCiphertext)
		evaluator.Rescale(operateCiphertext, params.DefaultScale(), operateCiphertext)
		if operateCiphertext, err = evaluator.EvaluatePoly(operateCiphertext, approxSigmoidReversed, params.DefaultScale()); err != nil {
			panic(err)
		}
		evaluator.MulRelin(operateCiphertext, evaluator.DropLevelNew(dataCiphertext, 7), operateCiphertext)
		evaluator.Rescale(operateCiphertext, params.DefaultScale(), operateCiphertext)
		// 明文===============================================================================================================
		for i := 0; i < len(operateData); i++ {
			for j := 0; j < len(operateData[i]); j++ {
				operateData[i][j] = mathFunc.Sigmoid(-operateData[i][j]) * paddedData[i][j]
			}
		}
		//printDebug(params, dataC, operateData[0], decryptor, encoder)

		// 对每一列的平均值然后作为梯度乘上学习率来更新权重
		// 密文===============================================================================================================
		learningRate := 10 / (float64(t+1) * float64(len(trainData)))
		for i, r := 4, (1 << 4); i < 4+10; i, r = i+1, r<<1 {
			evaluator.Add(operateCiphertext, evaluator.RotateNew(operateCiphertext, r), operateCiphertext)
		}
		evaluator.MultByConst(operateCiphertext, learningRate, operateCiphertext)
		evaluator.Rescale(operateCiphertext, params.DefaultScale(), operateCiphertext)
		decodedCiphertext := dataProc.DecCiphertext(operateCiphertext, params.LogSlots(), params, sk)
		operateCiphertext, _ = dataProc.EncRowVectorF64(typeConv.Complex128ToFloat64(decodedCiphertext), 1, params.MaxLevel(), params.DefaultScale(), params.LogSlots(), params, pk)
		evaluator.Add(weightsCiphertext, operateCiphertext, weightsCiphertext)
		// 明文===============================================================================================================
		for i := 0; i < len(weights); i++ {
			sum := 0.0
			for j := 0; j < len(operateData); j++ {
				sum += operateData[j][i]
			}
			weights[i] += learningRate * sum
		}
		fmt.Println("Iteration:", t+1)
		ciphertextDebug(params, weightsCiphertext, weights, decryptor, encoder)
	}

	// 输出权重
	plainWeightsDecoded := dataProc.DecCiphertext(weightsCiphertext, params.LogSlots(), params, sk)
	plainWeightsFloat := typeConv.Complex128ToFloat64(plainWeightsDecoded[0:len(testData[0])])
	fmt.Println()
	fmt.Println("WeightsTest:", plainWeightsFloat[0:len(testData[0])])
	fmt.Println("WeightsWant:", weights[0:len(testData[0])])
	fmt.Println()

	// 对测试集进行预测, 并计算准确率
	// 分别输出预测值和真实值(均重新转换为0和1)
	correctWant, correctTest := predictDebug(testData, weights, plainWeightsFloat, true)
	fmt.Println()
	fmt.Println("On test data:")
	fmt.Println("AccuracyTest:", float64(correctTest)/float64(len(testData)))
	fmt.Println("AccuracyWant:", float64(correctWant)/float64(len(testData)))
	fmt.Println()
	// 对训练集进行预测, 并计算准确率
	correctWant, correctTest = predictDebug(trainData, weights, plainWeightsFloat, false)
	fmt.Println()
	fmt.Println("On train data:")
	fmt.Println("AccuracyTest:", float64(correctTest)/float64(len(trainData)))
	fmt.Println("AccuracyWant:", float64(correctWant)/float64(len(trainData)))
	fmt.Println()
}
