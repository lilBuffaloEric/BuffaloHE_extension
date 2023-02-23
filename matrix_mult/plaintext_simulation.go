package matrix_mult

import (
	"errors"
	"math"
)

func Sigma_permute(A [][]float64, d int) (rslt [][]float64) {
	i := 0
	j := 0
	rslt = make([][]float64, d)
	for i = 0; i < d; i++ {
		rslt[i] = make([]float64, d)
	}
	for i = 0; i < d; i++ {
		for j = 0; j < d; j++ {
			rslt[i][j] = A[i][(i+j)%d]
		}
	}
	return
}

func Tao_permute(A [][]float64, d int) (rslt [][]float64) {
	i := 0
	j := 0
	rslt = make([][]float64, d)
	for i = 0; i < d; i++ {
		rslt[i] = make([]float64, d)
	}
	for i = 0; i < d; i++ {
		for j = 0; j < d; j++ {
			rslt[i][j] = A[(i+j)%d][j]
		}
	}
	return
}

func Phi_permute(A [][]float64, d int) (rslt [][]float64) {
	i := 0
	j := 0
	rslt = make([][]float64, d)
	for i = 0; i < d; i++ {
		rslt[i] = make([]float64, d)
	}
	for i = 0; i < d; i++ {
		for j = 0; j < d; j++ {
			rslt[i][j] = A[i][(j+1)%d]
		}
	}
	return
}

func Psi_permute(A [][]float64, d int) (rslt [][]float64) {
	i := 0
	j := 0
	rslt = make([][]float64, d)
	for i = 0; i < d; i++ {
		rslt[i] = make([]float64, d)
	}
	for i = 0; i < d; i++ {
		for j = 0; j < d; j++ {
			rslt[i][j] = A[(i+1)%d][j]
		}
	}
	return
}

func Zero_padding(A [][]float64, row int, col int) (rslt [][]float64) {
	i := 0
	if row == col {
		return A
	}
	if row > col {
		rslt = make([][]float64, row)
		for i = 0; i < row; i++ {
			rslt[i] = make([]float64, row)
			copy(rslt[i], A[i])
		}
		return
	}
	if col > row {
		rslt = make([][]float64, col)
		for i = 0; i < col; i++ {
			rslt[i] = make([]float64, col)
			copy(rslt[i], A[i])
		}
		return
	}
	return nil
}

func Row_ordering(a []float64) (A [][]float64, err error) {
	n := len(a)
	d := int(math.Sqrt(float64(n)))
	if d*d != n {
		return nil, errors.New("err: input vector's length can't be squared into integer")
	}
	A = make([][]float64, d)
	for i := 0; i < d; i++ {
		A[i] = make([]float64, d)
		for j := 0; j < d; j++ {
			A[i][j] = a[i*d+j]
		}
	}
	return
}

func Row_ordering_multiple(a []float64, g int) (A [][]float64, err error) {
	if len(a)%g != 0 {
		return nil, errors.New("err: the length of the input vector is not a multiple of k")
	}
	n := len(a) / g
	d := int(math.Sqrt(float64(n)))
	if n%d != 0 {
		return nil, errors.New("err: input vector's length can't be squared into integer")
	}
	A = make([][]float64, d)
	for i := 0; i < g*d; i++ {
		A[i] = make([]float64, d)
	}
	for k := 0; k < g; k++ {
		for i := 0; i < d; i++ {
			for j := 0; j < d; j++ {
				A[k*d+i][j] = a[g*(d*i+j)+k]
			}
		}
	}
	return
}

func Row_orderingInv(A [][]float64) (a []float64, err error) {
	/*
		if len(A) != len(A[0]) {
			return nil, errors.New("err: input matrix is not a square one")
		}
	*/
	row := len(A)
	col := len(A[0])
	a = make([]float64, row*col)
	for i := 0; i < row; i++ {
		for j := 0; j < col; j++ {
			a[i*col+j] = A[i][j]
		}
	}
	return
}

func Row_orderingInvZeroPad(A [][]float64, d int) (a []float64, err error) {
	if d < len(A) || d < len(A[0]) {
		return nil, errors.New("square dimension smaller than original colsize and rowsize")
	}
	row := len(A)
	col := len(A[0])
	a = make([]float64, d*d)
	for i := 0; i < row; i++ {
		for j := 0; j < col; j++ {
			a[i*col+j] = A[i][j]
		}
	}
	return
}

func Row_orderingInv_multiple(A [][][]float64) (a []float64, err error) {
	g := len(A)
	if len(A[0]) != len(A[0][0]) {
		return nil, errors.New("err: input matrices is not square matrices")
	}
	d := len(A[0])
	a = make([]float64, d*d*g)
	for k := 0; k < g; k++ {
		for i := 0; i < d; i++ {
			for j := 0; j < d; j++ {
				a[g*(d*i+j)+k] = A[k][i][j]
			}
		}
	}
	return
}

func Col_ordering(a []float64) (A [][]float64, err error) {
	n := len(a)
	d := int(math.Sqrt(float64(n)))
	if d*d != n {
		return nil, errors.New("err: input vector's length can't be squared into integer")
	}
	A = make([][]float64, d)
	for i := 0; i < d; i++ {
		A[i] = make([]float64, d)
	}
	for i := 0; i < d; i++ {
		for j := 0; j < d; j++ {
			A[j][i] = a[i*d+j]
		}
	}
	return
}

func Col_orderingInv(A [][]float64) (a []float64, err error) {
	if len(A) != len(A[0]) {
		return nil, errors.New("err: input matrix is not a square one")
	}
	d := len(A)
	a = make([]float64, d*d)
	for i := 0; i < d; i++ {
		for j := 0; j < d; j++ {
			a[i*d+j] = A[j][i]
		}
	}
	return
}

func Hadmard_matrix_mult(A [][]float64, B [][]float64) (C [][]float64, err error) {
	if len(A) != len(B) {
		return nil, errors.New("err: A's row number does not match B's row number")
	}
	if len(A[0]) != len(B[0]) {
		return nil, errors.New("err: A's column number does not match B's column number")
	}
	row := len(A)
	col := len(A[0])
	C = make([][]float64, row)
	for i := 0; i < row; i++ {
		C[i] = make([]float64, col)
		for j := 0; j < col; j++ {
			C[i][j] = A[i][j] * B[i][j]
		}
	}
	return
}

func Hadmard_matrix_add(A [][]float64, B [][]float64) (C [][]float64, err error) {
	if len(A) != len(B) {
		return nil, errors.New("err: A's row number does not match B's row number")
	}
	if len(A[0]) != len(B[0]) {
		return nil, errors.New("err: A's column number does not match B's column number")
	}
	row := len(A)
	col := len(A[0])
	C = make([][]float64, row)
	for i := 0; i < row; i++ {
		C[i] = make([]float64, col)
		for j := 0; j < col; j++ {
			C[i][j] = A[i][j] + B[i][j]
		}
	}
	return
}

func SquareMatrix_product_permute_version(A [][]float64, B [][]float64) (C [][]float64, err error) {
	if len(A) != len(B) {
		return nil, errors.New("err: A's row number does not match B's row number")
	}
	if len(A[0]) != len(B[0]) {
		return nil, errors.New("err: A's column number does not match B's column number")
	}
	row := len(A)
	col := len(A[0])
	if row != col {
		return nil, errors.New("err: inpute Matrices is not Square Matrices")
	}
	d := row
	C = make([][]float64, d)
	for i := 0; i < d; i++ {
		C[i] = make([]float64, d)
	}
	for k := 0; k < d; k++ {
		colshiftA := Sigma_permute(A, d)
		//auxio.Print_matrix_f64_full_2d(colshiftA, d, d)
		rowshiftB := Tao_permute(B, d)
		//auxio.Print_matrix_f64_full_2d(rowshiftB, d, d)
		for i := 0; i < k; i++ {
			colshiftA = Phi_permute(colshiftA, d)
			//auxio.Print_matrix_f64_full_2d(colshiftA, d, d)
			rowshiftB = Psi_permute(rowshiftB, d)
			//auxio.Print_matrix_f64_full_2d(rowshiftB, d, d)
		}
		var mult_rslt [][]float64
		mult_rslt, err = Hadmard_matrix_mult(colshiftA, rowshiftB)
		//auxio.Print_matrix_f64_full_2d(mult_rslt, d, d)
		if err != nil {
			return nil, err
		}
		C, err = Hadmard_matrix_add(C, mult_rslt)
		//auxio.Print_matrix_f64_full_2d(C, d, d)
		if err != nil {
			return nil, err
		}
	}
	return
}

func Gen_sigma_diagonalVecotrs(d int) (U map[int][]float64, err error) {
	if d <= 0 {
		return nil, errors.New("dimension d <= 0 ")
	}
	U = make(map[int][]float64, 2*d-1)
	for k := -d + 1; k < d; k++ {
		U[k] = make([]float64, d*d)
		for i := 0; i < d*d; i++ {
			if k >= 0 {
				if (i-d*k) >= 0 && (i-d*k) < (d-k) {
					U[k][i] = 1
				} else {
					U[k][i] = 0
				}
			} else {
				if (i-(d+k)*d) >= -k && (i-(d+k)*d) < d {
					U[k][i] = 1
				} else {
					U[k][i] = 0
				}
			}
		}
	}
	return
}

func Gen_tao_diagonalVectors(d int) (U map[int][]float64, err error) {
	if d <= 0 {
		return nil, errors.New("dimension d <= 0 ")
	}
	U = make(map[int][]float64, d)
	for k := 0; k < d; k++ {
		U[d*k] = make([]float64, d*d)
		for i := 0; i < d; i++ {
			U[d*k][k+d*i] = 1
		}
	}
	return
}

func Gen_colShift_diagonalVectors(d int, k int) (U map[int][]float64, err error) {
	if k < 0 || k >= d || d <= 0 { //  k < 1 || k >= d || d <= 0
		return nil, errors.New("dimension d <= 0 or k < 1 or k >= d")
	}
	U = make(map[int][]float64, 2)
	U[k] = make([]float64, d*d)
	U[k-d] = make([]float64, d*d)
	for i := 0; i < d*d; i++ {
		if (0 <= i%d) && (i%d < (d - k)) {
			U[k][i] = 1
		} else {
			U[k][i] = 0
		}
		if (d-k) <= (i%d) && (i%d) < d {
			U[k-d][i] = 1
		} else {
			U[k-d][i] = 0
		}
	}
	return
}

func Gen_rowShift_diagonalVectors(d int, k int) (U map[int][]float64, err error) {
	if k < 0 || k >= d || d <= 0 { //  k < 1 || k >= d || d <= 0
		return nil, errors.New("dimension d <= 0 or k < 1 or k >= d")
	}
	U = make(map[int][]float64, 1)
	U[k*d] = make([]float64, d*d)
	for i := 0; i < d*d; i++ {
		U[k*d][i] = 1
	}
	return
}

func Gen_transpose_diagonalVectors(d int) (U map[int][]float64, err error) {
	if d <= 0 {
		return nil, errors.New("dimension d <= 0 ")
	}
	U = make(map[int][]float64, 2*d-1)
	for i := -d + 1; i < d; i++ {
		U[(d-1)*i] = make([]float64, d*d)
		for l := 0; l < d*d; l++ {
			k := (l - i) % (d + 1)
			j := (l - i) / (d + 1)
			if i >= 0 {
				if k == 0 && 0 <= j && j < d-i {
					U[(d-1)*i][l] = 1
				} else {
					U[(d-1)*i][l] = 0
				}
			} else {
				if k == 0 && -i <= j && j < d {
					U[(d-1)*i][l] = 1
				} else {
					U[(d-1)*i][l] = 0
				}
			}
		}
	}
	return
}

func Gen_trans_C_tao_diagonalVectors(d int) (U map[int][]float64, err error) {
	var U_tao map[int][]float64
	var U_trans map[int][]float64
	var sum float64
	U_tao, err = Gen_tao_diagonalVectors(d)
	if err != nil {
		return nil, err
	}
	U_trans, err = Gen_transpose_diagonalVectors(d)
	if err != nil {
		return nil, err
	}
	M_tao := DiagonalVectors2Matrix(U_tao, d*d)
	M_trans := DiagonalVectors2Matrix(U_trans, d*d)
	M := make([][]float64, d*d)
	for i := 0; i < d*d; i++ {
		M[i] = make([]float64, d*d)
	}
	// Do Matrix Mult
	for i := 0; i < d*d; i++ {
		for j := 0; j < d*d; j++ {
			sum = 0
			for k := 0; k < d*d; k++ {
				// data[i]
				sum += M_tao[j][k] * M_trans[k][i]
			}
			M[j][i] = sum
		}
	}
	// Transform to Diagonal Form
	U, err = Matrix2DiagonalVectors(M)
	return
}

func plainDotMult(a []int, b []int) int {
	n := len(a)
	res := 0
	for i := 0; i < n; i++ {
		res += a[i] * b[i]
	}
	return res

}

func IsInMap(key int, M map[int]interface{}) int {
	_, ok := M[key]
	if ok {
		return 1
	} else {
		return 0
	}
}

func isAllzero(arr []float64) int {
	for i := 0; i < len(arr); i++ {
		if arr[i] != 0 {
			return 0
		}
	}
	return 1
}

func DiagonalVectors2Matrix(U map[int][]float64, n int) (M [][]float64) {
	M = make([][]float64, n)
	for i := 0; i < n; i++ {
		M[i] = make([]float64, n)
	}
	for key, data := range U {
		i := 0
		j := (key + n) % n
		for ; i < n; i++ {
			M[i][j] = data[i]
			j = (j + 1) % n
		}
	}
	return
}

func Matrix2DiagonalVectors(M [][]float64) (U map[int][]float64, err error) {
	if len(M[0]) != len(M) {
		return nil, errors.New("input matrix is not a square one")
	}
	U = make(map[int][]float64, len(M))
	for l := 0; l < len(M); l++ {
		U[l] = make([]float64, len(M))
		j := l
		for i := 0; i < len(M); i++ {
			U[l][i] = M[i][j]
			j = (j + 1) % len(M)
		}
		if isAllzero(U[l]) == 1 {
			delete(U, l)
		}
	}
	return
}
