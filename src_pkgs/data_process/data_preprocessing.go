package data_process

import (
	"log"
	"math/rand"
)

// 对样本的每个特征进行线性归一化(Max-Min)
func NormalizeMaxMin(data [][]int) [][]float64 {
	minVals := make([]float64, len(data[0]))
	maxVals := make([]float64, len(data[0]))
	for i := 0; i < len(data[0]); i++ {
		minVals[i] = float64(data[0][i])
		maxVals[i] = float64(data[0][i])
	}
	for _, record := range data {
		for i, v := range record {
			if float64(v) < minVals[i] {
				minVals[i] = float64(v)
			}
			if float64(v) > maxVals[i] {
				maxVals[i] = float64(v)
			}
		}
	}
	floatData := make([][]float64, len(data))
	for i := 0; i < len(data); i++ {
		floatData[i] = make([]float64, len(data[i]))
		for j := 0; j < len(data[i]); j++ {
			floatData[i][j] = (float64(data[i][j]) - minVals[j]) / (maxVals[j] - minVals[j])
		}
	}
	return floatData
}

// 将用0将样本的行数和列数分别填充到2的幂次方
func Padding(data [][]float64, logColSize int, logRowSize int) [][]float64 {
	col := len(data[0])
	row := len(data)
	if col > (1<<logColSize) || row > (1<<logRowSize) {
		log.Fatal("The number of rows and columns are greater than 2^logRows and 2^logCols")
	}
	newRowSize := 1 << logRowSize
	newColSize := 1 << logColSize
	newData := make([][]float64, newRowSize)
	for i := 0; i < newRowSize; i++ {
		newData[i] = make([]float64, newColSize)
		if i < row {
			copy(newData[i], data[i])
		}
	}
	return newData
}

// 对样本指定标签列从{0,1}转换为{-1,1}并将其位置重新调整到指定列
func ModifyLabel(data [][]float64, labelCol int, destCol int) [][]float64 {
	for i := 0; i < len(data); i++ {
		data[i][destCol], data[i][labelCol] = data[i][labelCol], data[i][destCol]
		if data[i][destCol] == 0 {
			data[i][destCol] = -1
		}
	}
	return data
}

// 根据指定标签列，将负样本其他特征值列取反
func InvertFeatures(data [][]float64, labelCol int) [][]float64 {
	for i := 0; i < len(data); i++ {
		if data[i][labelCol] == -1 {
			for j := 0; j < len(data[i]); j++ {
				if j != labelCol {
					data[i][j] = -data[i][j]
				}
			}
		}
	}
	return data
}

// 随机划分测试集和训练集
func SplitData(data [][]float64, percent float64) (testData [][]float64, trainData [][]float64) {
	rand.Shuffle(len(data), func(i, j int) {
		data[i], data[j] = data[j], data[i]
	})
	testData = make([][]float64, int(float64(len(data))*percent))
	for i := 0; i < len(testData); i++ {
		testData[i] = make([]float64, len(data[i]))
		copy(testData[i], data[i])
	}
	trainData = make([][]float64, len(data)-len(testData))
	for i := 0; i < len(trainData); i++ {
		trainData[i] = make([]float64, len(data[i+len(testData)]))
		copy(trainData[i], data[i+len(testData)])
	}
	return testData, trainData
}
