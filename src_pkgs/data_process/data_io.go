package data_process

import (
	"encoding/csv"
	"log"
	"os"
	"strconv"
)

// 从CSV文件中读取所有样本
func ReadAllFromCSV(path string) [][]int {
	file, err := os.Open(path)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	reader := csv.NewReader(file)
	reader.Comma = ','
	records, err := reader.ReadAll()
	if err != nil {
		log.Fatal(err)
	}
	intRecords := make([][]int, len(records[1:]))
	for i, record := range records[1:] {
		intRecords[i] = make([]int, len(record))
		for j, v := range record {
			intVal, err := strconv.Atoi(v)
			if err != nil {
				log.Fatal(err)
			}
			intRecords[i][j] = intVal
		}
	}
	return intRecords
}
