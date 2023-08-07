package auxiliary_io

/*
 * This part of package auxiliary_io aims to implement IO of csv file. It has been tested that library "encoding/csv"
 * has worse performance than simply using library "bufio", since it does not support loacting and reading a specific
 * row (neither can a Scanner of bufio, but the latter is faster). Although seeking some specific offset of bytes in the file
 * is available with *os.File object, correctly locating the starting byte of some line in a csv file is not so easy.
 * So here we use a Scanner of bufio to scan the file until we reach the target line.
 */

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"strings"
)

// Retrieve one CSV block using the rowIdx and ColTdx as location. This routine will retrieve a rowSize x colSize matrix with (rowIdx,colIdx) as the first entry.
func GetOneCSVSubblock(csvfilepath string, rowIdx int, colIdx int, rowSize int, colSize int) (subblock [][]string, err error) {
	var csvfile *os.File
	// var line []string
	csvfile, err = os.Open(csvfilepath)
	if err != nil {
		log.Fatalln("Couldn't open the csv file", err)
	}
	defer csvfile.Close()
	scanner := bufio.NewScanner(csvfile)
	scanner.Split(bufio.ScanLines)
	for i := 0; i < rowIdx; i++ {
		scanner.Scan()
		// test := scanner.Text()
		// fmt.Printf("%s", test)
	}
	subblock = make([][]string, rowSize)
	for i := 0; i < rowSize; i++ {
		subblock[i] = make([]string, colSize)
		scanner.Scan()
		line_str := scanner.Text()
		fields := strings.Split(line_str, ",")
		copy(subblock[i], fields[colIdx:colIdx+colSize])
	}
	return
}

// Export a float64 matrix to a csv file.
func ExportToCSV(data [][]float64, filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)

	for _, row := range data {
		strRow := make([]string, len(row))
		for i, val := range row {
			strRow[i] = fmt.Sprintf("%f", val)
		}
		if err := writer.Write(strRow); err != nil {
			return err
		}
	}

	writer.Flush()

	return nil
}
