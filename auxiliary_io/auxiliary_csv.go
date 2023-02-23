package auxiliary_io

/*
 * This part of package auxiliary_io is mean to implement IO of csv file. It has been tested that library "encoding/csv"
 * has worse performance than simply using library "bufio", since it does not support loacting and reading a specific
 * row (neither can a Scanner of bufio, but the latter is faster). Although seeking some specific offset of bytes in the file
 * is available with *os.File object, correctly locating the starting byte of some line in a csv file is not so easy.
 * So here we use a Scanner of bufio to scan the file until we reach the target line.
 */

import (
	"bufio"
	"log"
	"os"
	"strings"
)

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

	/* version 2 (Parse err occurs in reader.Read(), basically the function cant get the right fields for lines)
	scanner := bufio.NewScanner(csvfile)
	scanner.Split(bufio.ScanLines)
	for i := 0; i < rowIdx; i++ {
		scanner.Scan()
		test := scanner.Text()
		fmt.Printf("%s", test)
	}
	_, err = csvfile.Seek(int64(scanner.Bytes()[0]), 0)
	if err != nil {
		panic(err)
	}
	reader := csv.NewReader(csvfile)
	subblock = make([][]string, rowSize)
	for i := 0; i < rowSize; i++ {
		subblock[i] = make([]string, colSize)
		line, err = reader.Read()
		if err != nil {
			return nil, err
		}
		copy(subblock[i], line[colIdx:colIdx+colSize])
	}
	*/

	/* version 1 (No signifficant err, but performs badly when the CSV file is big)
	for i := 0; i < rowIdx; i++ {
		_, err = reader.Read()
	}
	if err != nil {
		return nil, err
	}

	subblock = make([][]string, rowSize)
	for i := 0; i < rowSize; i++ {
		subblock[i] = make([]string, colSize)
		line, err = reader.Read()
		if err != nil {
			return nil, err
		}
		copy(subblock[i], line[colIdx:colIdx+colSize])
	}
	*/

	return
}
