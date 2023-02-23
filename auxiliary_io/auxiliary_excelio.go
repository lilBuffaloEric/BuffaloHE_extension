package auxiliary_io

/*
 * This part of package auxiliary_io works Extremly slow in debugg version... This should be deprecated.
 * This is orginally intended to do some IO of excel files using library excelize. However, it does not
 * perform well.
 */
import (
	"strconv"

	"github.com/xuri/excelize/v2"
)

type Subblock_str struct {
	Row    int
	Col    int
	Values [][]string
}

type Subblock_f64 struct {
	Row    int
	Col    int
	Values [][]float64
}

func GetSubblock(file *excelize.File, sheetname string, ltop string, rbot string) (sb *Subblock_str, err error) {
	// var rows [][]string
	var lcol, rcol int
	var lrow, rrow int
	var rowsize, colsize int
	var cell string
	// rows, err = file.GetRows(sheetname)
	if err != nil {
		return sb, err
	}
	lcol, lrow, err = excelize.CellNameToCoordinates(ltop)
	if err != nil {
		return sb, err
	}
	rcol, rrow, err = excelize.CellNameToCoordinates(rbot)
	if err != nil {
		return sb, err
	}
	rowsize = rrow - lrow + 1
	colsize = rcol - lcol + 1
	SB := Subblock_str{Row: rowsize, Col: colsize, Values: make([][]string, rowsize)}
	sb = &SB
	for i := lrow; i <= rrow; i++ {
		k := i - lrow
		sb.Values[k] = make([]string, colsize)
		for j := lcol; j <= rcol; j++ {
			l := j - lcol
			cell, _ = excelize.CoordinatesToCellName(j, i)
			sb.Values[k][l], err = file.GetCellValue(sheetname, cell)
		}
	}
	return

}

func SwitchBlock_str2f64(sbIn Subblock_str) (sbOut Subblock_f64, err error) {
	sbOut.Col = sbIn.Col
	sbOut.Row = sbIn.Row
	sbOut.Values = make([][]float64, sbIn.Row)
	// float, err := strconv.ParseFloat(string, 64)
	for i := 0; i < sbIn.Row; i++ {
		sbOut.Values[i] = make([]float64, sbIn.Col)
		for j := 0; j < sbIn.Col; j++ {
			sbOut.Values[i][j], err = strconv.ParseFloat(sbIn.Values[i][j], 64)
			if err != nil {
				return
			}
		}
	}
	return
}
