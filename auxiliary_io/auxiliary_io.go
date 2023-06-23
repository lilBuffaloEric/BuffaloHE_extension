package auxiliary_io

/*
 * This part of package auxiliary_io is mean to print data structures in matrix or vector form, and transform
 * the type of data in the structures. Notice that most of the function is designed for type float64 and string,
 * More flexible interface for any type of data (such as complex, int, ect...) is needed.
 */

import (
	"fmt"
	"strconv"
)

func Print_matrix_f64(arr []float64, row_num int, col_num int) {
	if row_num*col_num > len(arr) {
		panic("row_num * col_num is bigger than the input slice's length")
	}
	print_row_size := 4
	if print_row_size > row_num {
		print_row_size = row_num
	}
	print_col_size := 4
	if print_col_size > col_num {
		print_col_size = col_num
	}
	for i := 0; i < print_row_size; i++ {
		fmt.Printf("[")
		for j := 0; j < print_col_size; j++ {
			fmt.Printf("%.8f ", arr[col_num*i+j])
		}
		fmt.Printf(".........")
		for j := col_num - print_col_size; j < col_num; j++ {
			fmt.Printf("%.8f ", arr[col_num*i+j])
		}
		fmt.Printf("]\n")
	}
	fmt.Print(".........\n")
	for i := row_num - print_row_size; i < row_num; i++ {
		fmt.Printf("[")
		for j := 0; j < print_col_size; j++ {
			fmt.Printf("%.8f ", arr[col_num*i+j])
		}
		fmt.Printf(".........")
		for j := col_num - print_col_size; j < col_num; j++ {
			fmt.Printf("%.8f ", arr[col_num*i+j])
		}
		fmt.Printf("]\n")
	}
	fmt.Printf("\n")
}

func Print_vector_f64(arr []float64, size int) {
	print_ele_size := 4
	if print_ele_size > size {
		print_ele_size = size
	}

	fmt.Printf("[")
	for i := 0; i < print_ele_size; i++ {
		fmt.Printf("%.8f ", arr[i])
	}
	fmt.Printf(".........")
	for i := size - print_ele_size; i < size; i++ {
		fmt.Printf("%.8f ", arr[i])
	}
	fmt.Printf("]\n")
}

func Print_vector_f64_full(arr []float64, prec ...int) {
	fmt.Printf("[")
	if len(prec) > 0 {
		for i := 0; i < len(arr); i++ {
			fmt.Printf("%."+strconv.Itoa(prec[0])+"f ", arr[i])
		}
	} else {
		for i := 0; i < len(arr); i++ {
			fmt.Printf("%.8f ", arr[i])
		}
	}
	fmt.Printf("]\n")
}

func Print_matrix_f64_full(arr []float64, row_num int, col_num int) {
	if row_num*col_num > len(arr) {
		panic("row_num * col_num is bigger than the input slice's length")
	}
	for i := 0; i < row_num; i++ {
		fmt.Printf("[")
		for j := 0; j < col_num; j++ {
			fmt.Printf("%.8f ", arr[col_num*i+j])
		}
		fmt.Printf("]\n")
	}
	fmt.Printf("\n")
}

func Print_matrix_f64_full_2d(arr2d [][]float64, row_num int, col_num int, prec ...int) {
	fmt.Printf("\n")
	for i := 0; i < row_num; i++ {
		Print_vector_f64_full(arr2d[i], prec...)
	}
	fmt.Printf("\n")
}

func Print_matrix_f64_2d(arr2d [][]float64, row_num int, col_num int) {

	fmt.Printf("\n")
	print_row_size := 4
	if print_row_size > row_num {
		print_row_size = row_num
	}
	for i := 0; i < print_row_size; i++ {
		Print_vector_f64(arr2d[i], col_num)
	}
	fmt.Print(".........")
	for i := row_num - print_row_size; i < row_num; i++ {
		Print_vector_f64(arr2d[i], col_num)
	}
	fmt.Printf("\n")
}

func Switch1d_str2f64(str []string) (arrf64 []float64, err error) {
	arrf64 = make([]float64, len(str))
	for i := 0; i < len(str); i++ {
		arrf64[i], _ = strconv.ParseFloat(str[i], 64)
	}
	return
}

func Switch2d_str2f64(str [][]string) (mtrxf64 [][]float64, err error) {
	mtrxf64 = make([][]float64, len(str))
	for i := 0; i < len(str); i++ {
		mtrxf64[i] = make([]float64, len(str[i]))
		mtrxf64[i], err = Switch1d_str2f64(str[i])
		if err != nil {
			return
		}
	}
	return
}

func Switch3d_str2f64(str [][][]string) (t3df64 [][][]float64, err error) {
	t3df64 = make([][][]float64, len(str))
	for i := 0; i < len(str); i++ {
		t3df64[i] = make([][]float64, len(str[i]))
		t3df64[i], err = Switch2d_str2f64(str[i])
		if err != nil {
			return
		}
	}
	return
}

func Switch_int2uint64(Invec []int) (Outvec []uint64) {
	Outvec = make([]uint64, len(Invec))
	for i := 0; i < len(Invec); i++ {
		Outvec[i] = uint64(Invec[i])
	}
	return
}
