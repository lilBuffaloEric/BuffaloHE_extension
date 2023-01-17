package auxiliary_io

import (
	"fmt"
)

func Print_matrix_f64(arr []float64, row_num int, col_num int) {
	if row_num*col_num > len(arr) {
		panic("row_num * col_num is bigger than the input slice's length")
	}
	print_row_size := 4

	for i := 0; i < print_row_size; i++ {
		fmt.Print("[")
		for j := 0; j < col_num; j++ {
			fmt.Printf("%4f ", arr[col_num*i+j])
		}
		fmt.Print("]\n")
	}
	fmt.Print(".........")
	for i := row_num - print_row_size; i < row_num; i++ {
		fmt.Print("[")
		for j := 0; j < col_num; j++ {
			fmt.Printf("%4f ", arr[col_num*i+j])
		}
		fmt.Print("]\n")
	}
	fmt.Print("\n")
}

func Print_vector_f64(arr []float64, size int) {
	print_ele_size := 4
	fmt.Print("[")
	for i := 0; i < print_ele_size; i++ {
		fmt.Printf("%4f ", arr[i])
	}
	fmt.Print(".........")
	for i := size - print_ele_size; i < size; i++ {
		fmt.Printf("%4f ", arr[i])
	}
	fmt.Print("]\n")
}

func Print_vector_f64_full(arr []float64) {
	fmt.Print("[")
	for i := 0; i < len(arr); i++ {
		fmt.Printf("%4f ", arr[i])
	}
	fmt.Print("]\n")
}

func Print_matrix_f64_full(arr []float64, row_num int, col_num int) {
	if row_num*col_num > len(arr) {
		panic("row_num * col_num is bigger than the input slice's length")
	}
	for i := 0; i < row_num; i++ {
		fmt.Print("[")
		for j := 0; j < col_num; j++ {
			fmt.Printf("%4f ", arr[col_num*i+j])
		}
		fmt.Print("]\n")
	}
	fmt.Print("\n")
}

func Print_matrix_f64_full_2d(arr2d [][]float64, row_num int, col_num int) {
	fmt.Print("\n")
	for i := 0; i < row_num; i++ {
		Print_vector_f64_full(arr2d[i])
	}
	fmt.Print("\n")
}
