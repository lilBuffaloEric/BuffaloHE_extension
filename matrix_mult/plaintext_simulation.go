package matrix_mult

import (
	"errors"
	"fmt"
	"math"
	"sort"
)

type Point struct {
	X int
	Y int
}
type NonZeroEntry struct {
	NormalCoordinate  Point
	DiagVecCoordinate Point
	Choice            map[Point]int
	Trytime           int
	Rest              int
}

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
			a[i*col+j+i*(d-col)] = A[i][j]
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
	/*
		for i := -d + 1; i < d; i++ {
			U[((d-1)*i+d*d)%(d*d)] = make([]float64, d*d)
			for l := 0; l < d*d; l++ {
				k := (l - i) % (d + 1)
				j := (l - i) / (d + 1)
				if i >= 0 {
					if k == 0 && 0 <= j && j < d-i {
						U[((d-1)*i+d*d)%(d*d)][l] = 1
					} else {
						U[((d-1)*i+d*d)%(d*d)][l] = 0
					}
				} else {
					if k == 0 && -i <= j && j < d {
						U[((d-1)*i+d*d)%(d*d)][l] = 1
					} else {
						U[((d-1)*i+d*d)%(d*d)][l] = 0
					}
				}
			}
		}
	*/

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

func GenSigmaDiagnalDecomposeMatrices(MatrixDimension int, TargetMaxDiagonalNo int) (Z1s []map[int][]float64, Z2s []map[int][]float64, err error) {
	maxNo := TargetMaxDiagonalNo
	n := MatrixDimension
	n2 := int(math.Ceil(float64(n-1) / 2.0))
	dpNum1 := int(math.Ceil((float64(n-1) / float64(maxNo)))) - 1
	dpNum2 := int(math.Ceil((float64(n2) / float64(maxNo)))) - 1
	// gen original Sigma Matrix in diagonal structure.
	Z, err := Gen_sigma_diagonalVecotrs(n)
	if err != nil {
		return nil, nil, err
	}
	// Seperate Sigma Matrix into two Matrices, one contains NonZero-Digonals with |[No.]_{n^2}| > \lceil (n-1)/2 \rceil , the other contains the rest.
	Z1 := make(map[int][]float64)
	Z2 := make(map[int][]float64)
	Zl := make(map[int][]float64)
	Zr := make(map[int][]float64)
	Z1s = make([]map[int][]float64, dpNum1+1)
	Z2s = make([]map[int][]float64, dpNum2+1)
	recorder1 := make(map[int]int) // record the key-value: (submatrixNo,nonZeroDiagNo)
	recorder2 := make(map[int]int)

	for key, value := range Z {
		if math.Abs(float64(key)) > math.Ceil(float64(n-1)/2.0) {
			Z1[key] = value
		} else {
			Z2[key] = value
		}
	}

	if (n % 2) == 0 {
		Z1[-int(math.Ceil(float64(n-1)/2.0))] = Z[-int(math.Ceil(float64(n-1)/2.0))]
		delete(Z2, -int(math.Ceil(float64(n-1)/2.0)))
	}

	for i := 0; i < n; i++ {
		if i <= n2 {
			if Z1[i-n] != nil {
				recorder1[i] = i - n
			}
			if Z2[i] != nil {
				recorder2[i] = i
			}

		} else {
			if Z1[i] != nil {
				recorder1[i] = i
			}
			if Z2[i-n] != nil {
				recorder2[i] = i - n
			}
		}
	}

	if dpNum1 == 0 {
		Z1s[0] = DeepIntFloat64MapCopy(Z) // Z
		Z2s = nil
		return
	}
	for i := 0; i < dpNum1; i++ {
		for subMtrxNo, key := range recorder1 { // go through each submatrix.
			value := Z1[key]
			if value == nil {
				continue
			}
			keyAbs := int(math.Abs(float64(key)))
			keySign := 1
			if key != 0 {
				keySign = key / keyAbs
			}
			var jl int
			var jr int
			var start int
			var end int
			var startl int
			var endl int
			var startr int
			var endr int
			if keyAbs >= maxNo {
				jl = keySign * maxNo
				jr = key - keySign*maxNo
			} else {
				jl = 0
				jr = key
			}
			if Zl[jl] == nil {
				Zl[jl] = make([]float64, n*n)
			}
			if Zr[jr] == nil {
				Zr[jr] = make([]float64, n*n)
			}

			if keySign == -1 {
				distance := -key + jr
				start = subMtrxNo*n + keyAbs
				end = (subMtrxNo + 1) * n
				startr = start - distance
				endr = end - distance
				startl = start
				endl = end
			} else {
				distance := key - jr
				start = subMtrxNo * n
				end = start + (n - key)
				startr = start + distance
				endr = end + distance
				startl = start
				endl = end
			}
			copy(Zr[jr][startr:endr], value[start:end])
			copy(Zl[jl][startl:endl], value[start:end])
			// update the NonZeroDiagNo. of subMtrxNo in recorder1.
			recorder1[subMtrxNo] = jr
		}
		// update the result in Z1s:
		// FIXME!!!!!!!!!!!!!!! Can't do Deep copy of map
		Z1s[i] = DeepIntFloat64MapCopy(Zl) //Zl
		if i == dpNum1-1 {
			Z1s[i+1] = DeepIntFloat64MapCopy(Zr) //Zr
		} else {
			Z1 = DeepIntFloat64MapCopy(Zr) //Zr
		}

		for k := range Zl {
			delete(Zl, k)
		}
		for k := range Zr {
			delete(Zr, k)
		}
	}

	if dpNum2 == 0 {
		Z2s[0] = DeepIntFloat64MapCopy(Z2) //Z2
		return
	}
	for i := 0; i < dpNum2; i++ {
		for subMtrxNo, key := range recorder2 {
			value := Z2[key]
			if value == nil {
				continue
			}
			keyAbs := int(math.Abs(float64(key)))
			keySign := 1
			if key != 0 {
				keySign = key / keyAbs
			}
			var jl int
			var jr int
			var start int
			var end int
			var startl int
			var endl int
			var startr int
			var endr int
			if keyAbs >= maxNo {
				jl = keySign * maxNo
				jr = key - keySign*maxNo
			} else {
				jl = 0
				jr = key
			}
			if Zl[jl] == nil {
				Zl[jl] = make([]float64, n*n)
			}
			if Zr[jr] == nil {
				Zr[jr] = make([]float64, n*n)
			}

			if keySign == -1 {
				distance := -key + jr
				start = subMtrxNo*n + keyAbs
				end = (subMtrxNo + 1) * n
				startr = start - distance
				endr = end - distance
				startl = start
				endl = end
			} else {
				distance := key - jr
				start = subMtrxNo * n
				end = start + (n - key)
				startr = start + distance
				endr = end + distance
				startl = start
				endl = end
			}
			copy(Zr[jr][startr:endr], value[start:end])
			copy(Zl[jl][startl:endl], value[start:end])
			// update the NonZeroDiagNo. of subMtrxNo in recorder1.
			recorder2[subMtrxNo] = jr
		}
		// update the result in Z1s:
		Z2s[i] = DeepIntFloat64MapCopy(Zl) //Zl
		if i == dpNum2-1 {
			Z2s[i+1] = DeepIntFloat64MapCopy(Zr) //Zr
		} else {
			Z2 = DeepIntFloat64MapCopy(Zr) //Zr
		}
		for k := range Zl {
			delete(Zl, k)
		}
		for k := range Zr {
			delete(Zr, k)
		}
	}
	return
}

func GenSigmaDiagnalDecomposeMatrices_Ver2(MatrixDimension int, TargetMaxDiagonalNo int) (Z1s []map[int][]float64, Z2s []map[int][]float64, err error) {
	maxNo := TargetMaxDiagonalNo
	n := MatrixDimension
	dpNum := int(math.Ceil((float64(n-1) / float64(maxNo)))) - 1
	Z, err := Gen_sigma_diagonalVecotrs(n)
	if err != nil {
		return nil, nil, err
	}
	Z1 := make(map[int][]float64)
	Z2 := make(map[int][]float64)
	Zl := make(map[int][]float64)
	Zr := make(map[int][]float64)
	Z1s = make([]map[int][]float64, dpNum+1)
	Z2s = make([]map[int][]float64, dpNum+1)
	recorder1 := make(map[int]int) // record the key-value: (submatrixNo,nonZeroDiagNo)
	recorder2 := make(map[int]int)
	for key, value := range Z {
		if key >= 0 {
			Z1[key] = value
		} else {
			Z2[key] = value
		}
	}
	for i := 0; i < n; i++ {
		recorder1[i] = i
		if i != 0 {
			recorder2[i] = i - n
		}
	}
	if dpNum == 0 {
		Z1s[0] = DeepIntFloat64MapCopy(Z) // Z
		Z2s = nil
		return
	}
	for i := 0; i < dpNum; i++ {
		for subMtrxNo, key := range recorder1 {
			value := Z1[key]
			if value == nil {
				continue
			}
			var jl int
			var jr int
			var start int
			var end int
			var startl int
			var endl int
			var startr int
			var endr int
			if key >= maxNo {
				jl = key - maxNo
				jr = maxNo
			} else {
				jl = key
				jr = 0
			}
			if Zl[jl] == nil {
				Zl[jl] = make([]float64, n*n)
			}
			if Zr[jr] == nil {
				Zr[jr] = make([]float64, n*n)
			}
			distance := key - jr
			start = subMtrxNo * n
			end = start + (n - key)
			startr = start + distance
			endr = end + distance
			startl = start
			endl = end
			copy(Zr[jr][startr:endr], value[start:end])
			copy(Zl[jl][startl:endl], value[start:end])
			recorder1[subMtrxNo] = jl
		}
		Z1s[dpNum+1-i-1] = DeepIntFloat64MapCopy(Zr)
		if i == dpNum-1 {
			Z1s[0] = DeepIntFloat64MapCopy(Zl)
		} else {
			Z1 = DeepIntFloat64MapCopy(Zl)
		}
		for k := range Zl {
			delete(Zl, k)
		}
		for k := range Zr {
			delete(Zr, k)
		}
	}
	for i := 0; i < dpNum; i++ {
		for subMtrxNo, key := range recorder2 {
			value := Z2[key]
			if value == nil {
				continue
			}
			keyAbs := int(math.Abs(float64(key)))
			keySign := -1
			if key != 0 {
				keySign = key / keyAbs
			}
			var jl int
			var jr int
			var start int
			var end int
			var startl int
			var endl int
			var startr int
			var endr int
			if keyAbs >= maxNo {
				jl = key - keySign*maxNo
				jr = keySign * maxNo
			} else {
				jl = key
				jr = 0
			}
			if Zl[jl] == nil {
				Zl[jl] = make([]float64, n*n)
			}
			if Zr[jr] == nil {
				Zr[jr] = make([]float64, n*n)
			}
			distance := -key + jr
			start = subMtrxNo*n + keyAbs
			end = (subMtrxNo + 1) * n
			startr = start - distance
			endr = end - distance
			startl = start
			endl = end
			copy(Zr[jr][startr:endr], value[start:end])
			copy(Zl[jl][startl:endl], value[start:end])
			recorder2[subMtrxNo] = jl
		}
		Z2s[dpNum+1-i-1] = DeepIntFloat64MapCopy(Zr)
		if i == dpNum-1 {
			Z2s[0] = DeepIntFloat64MapCopy(Zl)
		} else {
			Z2 = DeepIntFloat64MapCopy(Zl)
		}
		for k := range Zl {
			delete(Zl, k)
		}
		for k := range Zr {
			delete(Zr, k)
		}
	}
	return

}

func GenTauDiagonalDecomposeMatrices(MatrixDimension int, TargetMaxDiagonalNo int) (Ts []map[int][]float64, err error) {
	maxNo := TargetMaxDiagonalNo
	n := MatrixDimension
	if maxNo%n != 0 {
		return nil, errors.New("targetMaxDiagonalNois not divisible by MatrixDimension")
	}
	dpNum1 := int(math.Ceil((float64(n-1) / float64(maxNo/n)))) - 1
	// gen original Sigma Matrix in diagonal structure.
	T, err := Gen_tao_diagonalVectors(n)
	if err != nil {
		return nil, err
	}
	Tl := make(map[int][]float64)
	Tr := make(map[int][]float64)
	Ts = make([]map[int][]float64, dpNum1+1)
	if dpNum1 == 0 {
		Ts[0] = DeepIntFloat64MapCopy(T)
		return Ts, err
	}
	for t := 0; t < dpNum1; t++ {
		for i := 0; i <= maxNo/n; i++ {
			key := i * n
			value := T[key]
			if value == nil {
				continue
			}
			if Tl[0] == nil {
				Tl[0] = make([]float64, n*n)
			}
			if Tr[key] == nil {
				Tr[key] = make([]float64, n*n)
			}
			for idx, entry := range value {
				if entry == 0 {
					continue
				}
				Tl[0][idx] = 1
				Tr[key][idx] = 1
			}
		}
		for i := maxNo/n + 1; i < n; i++ {
			key := i * n
			value := T[key]
			fixkey := maxNo
			if value == nil {
				continue
			}
			if Tl[fixkey] == nil {
				Tl[fixkey] = make([]float64, n*n)
			}
			if Tr[key-fixkey] == nil {
				Tr[key-fixkey] = make([]float64, n*n)
			}
			for idx, entry := range value {
				if entry == 0 {
					continue
				}
				Tl[fixkey][idx] = 1
				Tr[key-fixkey][(key+idx-(key-fixkey)+n*n)%(n*n)] = 1
			}
		}
		Ts[t] = DeepIntFloat64MapCopy(Tl)
		if t == dpNum1-1 {
			Ts[t+1] = DeepIntFloat64MapCopy(Tr)
		} else {
			T = DeepIntFloat64MapCopy(Tr)
		}
		for k := range Tl {
			delete(Tl, k)
		}
		for k := range Tr {
			delete(Tr, k)
		}

	}
	return
}

func GenTauDiagonalDecomposeMatrices_Ver2(MatrixDimension int, TargetMaxDiagonalNo int) (Ts []map[int][]float64, err error) {
	maxNo := TargetMaxDiagonalNo
	n := MatrixDimension
	if maxNo%n != 0 {
		return nil, errors.New("targetMaxDiagonalNois not divisible by MatrixDimension")
	}
	dpNum1 := int(math.Ceil((float64(n-1) / float64(maxNo/n)))) - 1
	// gen original Sigma Matrix in diagonal structure.
	T, err := Gen_tao_diagonalVectors(n)
	if err != nil {
		return nil, err
	}
	Tl := make(map[int][]float64)
	Tr := make(map[int][]float64)
	Ts = make([]map[int][]float64, dpNum1+1)
	if dpNum1 == 0 {
		Ts[0] = DeepIntFloat64MapCopy(T)
		return Ts, err
	}
	for t := 0; t < dpNum1; t++ {
		for i := 0; i <= maxNo/n; i++ {
			key := i * n
			value := T[key]
			if value == nil {
				continue
			}
			if Tl[key] == nil {
				Tl[key] = make([]float64, n*n)
			}
			if Tr[0] == nil {
				Tr[0] = make([]float64, n*n)
			}
			for idx, entry := range value {
				if entry == 0 {
					continue
				}
				Tl[key][idx] = 1
				Tr[0][(idx+key+n*n)%(n*n)] = 1
			}
		}
		for i := maxNo/n + 1; i < n; i++ {
			key := i * n
			value := T[key]
			fixkey := maxNo
			if value == nil {
				continue
			}
			if Tl[key-fixkey] == nil {
				Tl[key-fixkey] = make([]float64, n*n)
			}
			if Tr[fixkey] == nil {
				Tr[fixkey] = make([]float64, n*n)
			}
			for idx, entry := range value {
				if entry == 0 {
					continue
				}
				Tl[key-fixkey][idx] = 1
				Tr[fixkey][(key+idx-(fixkey)+n*n)%(n*n)] = 1
			}
		}
		Ts[dpNum1+1-t-1] = DeepIntFloat64MapCopy(Tr)
		if t == dpNum1-1 {
			Ts[0] = DeepIntFloat64MapCopy(Tl)
		} else {
			T = DeepIntFloat64MapCopy(Tl)
		}
		for k := range Tl {
			delete(Tl, k)
		}
		for k := range Tr {
			delete(Tr, k)
		}

	}
	return
}

func Converge2DiagonalDecompose_Sigma(A [][]float64) (D []map[int][]float64, err error) {
	if len(A) != len(A[0]) {
		return nil, errors.New("input matrix is not a square one")
	}
	n := int(math.Sqrt(float64(len(A))))
	var U map[int][]float64
	U, err = Matrix2DiagonalVectors(A)
	if err != nil {
		return nil, err
	}

	diagVecs := make([]int, 0)
	diagVecs = append(diagVecs, 0)
	for i := 1; i <= int(math.Ceil(float64(n-1)/2.0)); i++ {
		diagVecs = append(diagVecs, i)
		diagVecs = append(diagVecs, -i+n*n)
	}

	D = make([]map[int][]float64, 2)
	D[0] = make(map[int][]float64)
	D[1] = make(map[int][]float64)
	for _, i := range diagVecs {
		D[0][i] = make([]float64, n*n)
		D[1][i] = make([]float64, n*n)
	}

	for i := -int(math.Ceil(float64(n-1) / 2.0)); i <= int(math.Ceil(float64(n-1)/2.0)); i++ {
		if i == 0 {
			continue
		}
		var j int
		if i > 0 {
			j = (i - n + n*n) % (n * n) // i>0, j<0 i-j = n
		} else {
			j = (i + n + n*n) % (n * n) // i<0, j>0 j-i = n
		}

		j1, j2, ok := findSums(diagVecs, j, n*n)
		if !ok {
			return nil, errors.New("cannot find sums")
		}
		for l, value := range U[j] {
			if value != 1 {
				continue
			}
			D[0][j1][(j+l-j1+n*n)%(n*n)] = 1
			D[1][j2][l] = 1
			z1 := (j + l - j1 + n*n) % (n * n)
			z2 := (j + l - j1 + n*n) % (n * n)
			for z1 != -1 {
				if U[(i+n*n)%(n*n)][z1] == 1 {
					D[0][0][(i+z1+n*n)%(n*n)] = 1
					D[1][(i+n*n)%(n*n)][z1] = 1
					z1 = (i + z1 + n*n) % (n * n)
				} else {
					z1 = -1
				}
			}
			for z2 != -1 {
				if U[(i+n*n)%(n*n)][(z2-i+n*n)%(n*n)] == 1 {
					D[0][(i+n*n)%(n*n)][(z2-i+n*n)%(n*n)] = 1
					D[1][0][(z2-i+n*n)%(n*n)] = 1
					z2 = (z2 - i + n*n) % (n * n)
				} else {
					z2 = -1
				}
			}

		}
		/*
			for l, value := range U[i] {
				if value != 1 {
					continue
				}
				occupied := false
				for _, A1diagVec := range D[0] {
					if A1diagVec[l] == 1 {
						occupied = true
						break
					}
				}
				if occupied {
					D[0][0][(i+l+n*n)%(n*n)] = 1
					D[1][i][l] = 1
				} else {
					D[0][i][l] = 1
					D[1][0][l] = 1
				}
			}
		*/
	}
	for l, value := range U[0] {
		if value != 1 {
			continue
		}
		D[0][0][l] = 1
		D[1][0][l] = 1
	}
	return
}

func Converge2DiagonalDecompose_Tao(T [][]float64) (D []map[int][]float64, err error) {
	if len(T) != len(T[0]) {
		return nil, errors.New("input matrix is not a square one")
	}
	n := int(math.Sqrt(float64(len(T))))
	var U map[int][]float64
	U, err = Matrix2DiagonalVectors(T)
	if err != nil {
		return nil, err
	}

	diagVecs := make([]int, 0)
	diagVecs = append(diagVecs, 0)
	for i := 1; i <= int(math.Ceil(float64(n-1)/2.0)); i++ {
		diagVecs = append(diagVecs, i*n)
	}

	D = make([]map[int][]float64, 2)
	D[0] = make(map[int][]float64)
	D[1] = make(map[int][]float64)
	for _, i := range diagVecs {
		D[0][i] = make([]float64, n*n)
		D[1][i] = make([]float64, n*n)
	}
	for i := 0; i <= int(math.Ceil(float64(n-1)/2.0)); i++ {
		for l, value := range U[i*n] {
			if value != 1 {
				continue
			}
			D[0][0][(i*n+l+n*n)%(n*n)] = 1
			D[1][i*n][l] = 1
		}
	}
	for i := int(math.Ceil(float64(n-1)/2.0)) + 1; i < n; i++ {
		i1n, i2n, ok := findSums(diagVecs, i*n, n*n)
		if !ok {
			return nil, errors.New("cannot find sums")
		}
		for l, value := range U[i*n] {
			if value != 1 {
				continue
			}
			D[0][i1n][(i*n+l-i1n+n*n)%(n*n)] = 1
			D[1][i2n][l] = 1
		}
	}

	return

}

func Converge2DiagonalDecompose(M [][]float64) (D []map[int][]float64, err error) {
	if len(M) != len(M[0]) {
		return nil, errors.New("input matrix is not a square one")
	}
	/*
		if diagonalvecNum >= len(M) {
			return nil, errors.New("not converge can be done since the required Number of diagonalvectors is too big")
		}
		if diagonalvecNum%2 == 0 {
			diagonalvecNum += 1
		}
	*/

	d := len(M)
	var U map[int][]float64
	U, err = Matrix2DiagonalVectors(M)
	if err != nil {
		return nil, err
	}
	keys := make([]int, 0, len(U))
	for k := range U {
		keys = append(keys, k)
	}
	subkeys := FindMinCompositeDiagVecSet(keys, d)
	// restkeys := difference(keys, subkeys)

	D = make([]map[int][]float64, 2)
	D[0] = make(map[int][]float64)
	D[1] = make(map[int][]float64)
	for _, i := range subkeys {
		D[0][i] = make([]float64, d)
		D[1][i] = make([]float64, d)
	}

	// Version 3
	subkeysMapNoZero := make(map[int]bool)
	for _, sk := range subkeys {
		if sk != 0 {
			subkeysMapNoZero[sk] = true
		}
	}
	subkeysNoZero := getKeys(subkeysMapNoZero)
	Opt1keysMap := make(map[int]bool)
	Opt2keysMap := make(map[int]bool)
	for _, s := range keys {
		s1, s2, ok := findSums(subkeysNoZero, s, d)
		if s1 == s2 && ok {
			Opt1keysMap[s] = true
		} else {
			if s == 0 && !ok {
				Opt1keysMap[s] = true
			} else {
				Opt2keysMap[s] = true
			}
		}
	}
	opt1keys := getKeys(Opt1keysMap)
	opt2keys := getKeys(Opt2keysMap)
	for _, r := range opt1keys {
		s1, s2, ok := findSums(subkeys, r, d)
		if !ok {
			return nil, err // debugtest
		}
		for i := 0; i < d; i++ {
			if U[r][i] == 1 {
				c := (r + i + d) % d
				l1 := (c - s1 + d) % d
				D[0][s1][l1] = 1
				D[1][s2][i] = 1
			}
		}
	}
	for _, r := range opt2keys {
		s1, s2, ok := findSums(subkeys, r, d)
		if !ok {
			return nil, err // debugtest
		}
		for i := 0; i < d; i++ {
			if U[r][i] == 1 {
				occupied := false
				c := (r + i + d) % d
				l1 := (c - s1 + d) % d
				for _, sk := range subkeys {
					if D[0][sk][l1] == 1 {
						occupied = true
						break
					}
				}
				if !occupied {
					D[0][s1][l1] = 1
					D[1][s2][i] = 1
				} else {
					c := (r + i + d) % d
					l2 := (c - s2 + d) % d
					D[0][s2][l2] = 1
					D[1][s1][i] = 1
				}

			}
		}
	}

	// Version 2, have some problem.
	/*
		for _, r := range restkeys {
			s1, s2, ok := findSums(subkeys, r, d)
			if !ok {
				return nil, err // debugtest
			}
			for i := 0; i < d; i++ {
				if U[r][i] == 1 {
					c := (r + i + d) % d
					l1 := (c - s1 + d) % d
					D[0][s1][l1] = 1
					D[1][s2][i] = 1
				}
			}
		}
		for _, k := range subkeys {
			for i := 0; i < d; i++ {
				if U[k][i] == 1 {
					occupied := false
					for _, k1 := range subkeys {
						if D[0][k1][i] == 1 {
							occupied = true
							break
						}
					}
					if occupied { // FIXME: Should implement kick out mechanism.
						c := (k + i + d) % d
						D[0][0][c] = 1
						D[1][k][i] = 1
					} else {
						D[0][k][i] = 1
						D[1][0][i] = 1
					}
				}
			}
		}
	*/

	/* Version 1, deprecated.

	keys := make([]int, 0, len(D[0]))
	for k := range D[0] {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	// coordinate conversion rule :
	// NormalCoordinate(i,j) -> diagonalCoordinate((j-i)mod d, (j-((j-i)mod d))mod d )
	// diagonalCoordinate(k,l) -> NormalCoordinate(((k+l)mod d)-k mod d , (k+l)mod d)
	for i := 0; i < d; i++ {
		idx := WhereistheOne(M[i])
		if idx != -1 && idx != i {
			occupied := false
			//localk := (idx - i + d) % d
			for k := 0; k < len(keys); k++ {
				keys[k] = (idx-keys[k]+d)%d - i
			}
			sort.Slice(keys, func(i, j int) bool {
				return math.Abs(float64(keys[i])) < math.Abs(float64(keys[j]))
			})

			for k := 0; k < len(keys); k++ {
				keys[k] = (-(keys[k] + i) + idx + d) % d
			}
			for _, k := range keys {
				occupied = false
				mapidx := (idx - k + d) % d // mapidx = (idx-k) mod d, compute the l=mapidx in diagonalCoordinate (k,l) of D[0]
				if D[0][k][mapidx] == 0 {
					// row := mapidx
					for k2, _ := range D[0] {
						// check to see if D[0]'s row is occupied already
						if D[0][k2][mapidx] != 0 {
							occupied = true
						}
					}
					if occupied {
						continue
					}
					key := (mapidx - i + d) % d // compute the diagonalCoordinate (k=key,l=no) of D[1]
					_, ok := D[1][key]
					if ok == false {
						occupied = true
						continue
					}
					no := (mapidx - key + d) % d
					for k2, _ := range D[1] {
						if D[1][k2][no] != 0 {
							occupied = true
						}
					}
					if occupied {
						continue
					}
					if D[1][key][no] == 0 {
						D[0][k][mapidx] = 1
						D[1][key][no] = 1
						break
					}
				}
			}
			if occupied {
				return D, errors.New("can not decompose the matrix.")
			}
		}
	}

	*/

	return
}

func WhereistheOne(Vec []float64) (idx int) {
	idx = -1
	for i := 0; i < len(Vec); i++ {
		if Vec[i] == 1 {
			idx = i
			return
		}
	}
	return
}

func FindMinCompositeDiagVecSet(SA []int, n int) []int {

	// Greedy Version
	S := make(map[int]bool)
	S[0] = true
	R := make(map[int]bool)
	for _, a := range SA {
		if a != 0 {
			R[a] = true
		}
	}
	for len(R) > 0 {
		r := getAnyKey(R)
		delete(R, r)

		found := false
		for s1 := range S {
			target := (r - s1 + n) % n
			if S[target] {
				found = true
				break
			}
		}

		if !found {
			S[r] = true
		}
	}
	Rkeys := difference(SA, getKeys(S))
	delete(S, 0)
	Skeys := getKeys(S)
	for _, key := range Skeys {
		_, _, ok := findSums(getKeys(S), key, n)
		if ok {
			delete(S, key)
			for _, rkey := range Rkeys {
				_, _, okk := findSums(getKeys(S), rkey, n)
				if !okk {
					S[key] = true
					break
				}
			}
		}
	}
	S[0] = true

	/* Greedy Version
	S := make(map[int]bool)
	S[0] = true
	R := make(map[int]bool)
	for _, a := range SA {
		if a != 0 {
			R[a] = true
		}
	}
	for len(R) > 0 {
		r := getAnyKey(R)
		delete(R, r)

		found := false
		for s1 := range S {
			target := (r - s1 + n) % n
			if S[target] {
				found = true
				break
			}
		}

		if !found {
			S[r] = true
		}
	}
	*/

	return getKeys(S)
}

func getAnyKey(m map[int]bool) int {
	for k := range m {
		return k
	}
	return -1
}

func getKeys(m map[int]bool) []int {
	keys := make([]int, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	return keys
}

func intersection(a, b []int) []int {
	// 将 a 转换成 map，方便判断是否包含某个元素
	m := make(map[int]bool)
	for _, v := range a {
		m[v] = true
	}

	var result []int
	for _, v := range b {
		// 如果 b 中的元素也在 a 中，则加入交集
		if m[v] {
			result = append(result, v)
		}
	}

	return result
}

func difference(a, b []int) []int {
	// 将 b 转换成 map，方便判断是否包含某个元素
	m := make(map[int]bool)
	for _, v := range b {
		m[v] = true
	}

	var result []int
	for _, v := range a {
		// 如果 a 中的元素不在 b 中，则加入差集
		if !m[v] {
			result = append(result, v)
		}
	}

	return result
}

func findSums(nums []int, s int, n int) (int, int, bool) {
	m := make(map[int]bool, len(nums))
	for i := 0; i < len(nums); i++ {
		m[nums[i]] = true
	}
	for _, num := range nums {
		complement := (s - num + n) % n
		if m[complement] {
			return num, complement, true
		}
	}
	return 0, 0, false
}

func findAllSums(nums []int, s int, n int) (total map[int]int) {
	total = make(map[int]int)
	m := make(map[int]bool, len(nums))
	for i := 0; i < len(nums); i++ {
		m[nums[i]] = true
	}
	for _, num := range nums {
		complement := (s - num + n) % n
		if m[complement] {
			total[num] = complement
		}
	}
	return
}

func SlowPlainMatrixMult(A [][]float64, B [][]float64, d int) (C [][]float64) {
	var sum float64
	C = make([][]float64, d)
	for i := 0; i < d; i++ {
		C[i] = make([]float64, d)
	}
	for i := 0; i < d; i++ {
		for j := 0; j < d; j++ {
			sum = 0
			for k := 0; k < d; k++ {
				// data[i]
				sum += A[j][k] * B[k][i]
			}
			C[j][i] = sum
		}
	}
	return
}

func SlowPlainMatrixAdd(A [][]float64, B [][]float64, d int) (C [][]float64) {
	C = make([][]float64, d)
	for i := 0; i < d; i++ {
		C[i] = make([]float64, d)
		for j := 0; j < d; j++ {
			C[i][j] = A[i][j] + B[i][j]
		}
	}
	return
}

func NewnonZeroEntry(k int, l int, dim int, diagVecs []int) (nonZeroEntry *NonZeroEntry) {
	nonZeroEntry = new(NonZeroEntry)
	nonZeroEntry.Trytime = 1
	nonZeroEntry.DiagVecCoordinate.X = k
	nonZeroEntry.DiagVecCoordinate.Y = l
	nonZeroEntry.NormalCoordinate.X = l
	nonZeroEntry.NormalCoordinate.Y = (k + l + dim) % dim
	nonZeroEntry.Choice = make(map[Point]int)
	total := findAllSums(diagVecs, k, dim)
	for s1, s2 := range total {
		diagVecpair := Point{X: s1, Y: s2}
		nonZeroEntry.Choice[diagVecpair] = 1
	}
	nonZeroEntry.Rest = len(total)
	return
}

func (nonZeroEntry *NonZeroEntry) ReInit(dim int) {
	/*
		diagVecpair, _ := latestEffectEntry.LatestChoice(dim)
		occupyRow := (latestEffectEntry.NormalCoordinate.Y - diagVecpair.X + dim) % dim
		occupyCol := (diagVecpair.Y + latestEffectEntry.NormalCoordinate.X + dim) % dim
		nonZeroEntry.Recover(occupyRow, occupyCol, dim)
	*/
	for diagVecpair_local, value := range nonZeroEntry.Choice {
		if value != 0 && value != 1 {
			nonZeroEntry.Choice[diagVecpair_local] = 1
			nonZeroEntry.Rest++
		}
	}
	nonZeroEntry.Trytime = 1
}

func (nonZeroEntry *NonZeroEntry) Forbid(row int, col int, dim int) (effected bool, safe bool) {
	effected = false
	safe = false
	for diagVecpair := range nonZeroEntry.Choice {
		var effected_local = false
		if (diagVecpair.X+row+dim)%dim == nonZeroEntry.NormalCoordinate.Y {
			nonZeroEntry.Choice[diagVecpair] = 0
			effected_local = true
			effected = true
		}
		if (col-diagVecpair.Y+dim)%dim == nonZeroEntry.NormalCoordinate.X {
			nonZeroEntry.Choice[diagVecpair] = 0
			effected_local = true
			effected = true
		}
		if nonZeroEntry.Choice[diagVecpair] == 1 {
			safe = true
		}
		if effected_local {
			nonZeroEntry.Rest--
		}
	}
	return
}

func (nonZeroEntry *NonZeroEntry) Recover(row int, col int, dim int) {
	for diagVecpair := range nonZeroEntry.Choice {
		var recovered = false
		if (diagVecpair.X+row+dim)%dim == nonZeroEntry.NormalCoordinate.Y {
			if nonZeroEntry.Choice[diagVecpair] != 1 {
				nonZeroEntry.Choice[diagVecpair] = 1
				recovered = true
			}
		}
		if (col-diagVecpair.Y+dim)%dim == nonZeroEntry.NormalCoordinate.X {
			if nonZeroEntry.Choice[diagVecpair] != 1 {
				nonZeroEntry.Choice[diagVecpair] = 1
				recovered = true
			}
		}
		if recovered {
			nonZeroEntry.Rest++
		}
	}
}

func (nonZeroEntry *NonZeroEntry) Choose(dim int) (occupyRow int, occupyCol int, ok bool) {
	ok = false
	if nonZeroEntry.Rest < 1 {
		return -1, -1, false
	}
	for diagVecpair, status := range nonZeroEntry.Choice {
		if status == 1 {
			ok = true
			nonZeroEntry.Trytime++
			nonZeroEntry.Rest--
			nonZeroEntry.Choice[diagVecpair] = nonZeroEntry.Trytime
			occupyRow = (nonZeroEntry.NormalCoordinate.Y - diagVecpair.X + dim) % dim
			occupyCol = (diagVecpair.Y + nonZeroEntry.NormalCoordinate.X + dim) % dim // get Norm(X,Y) and s2...
			break
		}
	}
	return
}

func (nonZeroEntry *NonZeroEntry) LatestChoice(dim int) (latestDiagVecpair Point, componentpair Point) {
	var latestValue = 1
	latestDiagVecpair = Point{X: -1, Y: -1}
	for diagVecpair, value := range nonZeroEntry.Choice {
		if latestValue < value {
			latestValue = value
			latestDiagVecpair = diagVecpair
		}
	}
	if latestDiagVecpair.X != -1 && latestDiagVecpair.Y != -1 {
		componentpair.X = (nonZeroEntry.NormalCoordinate.Y - latestDiagVecpair.X + dim) % dim
		componentpair.Y = nonZeroEntry.NormalCoordinate.X
	} else {
		componentpair.X = -1
		componentpair.Y = -1
	}

	return
}

func SortByRestChoice(slices []*NonZeroEntry) {
	sort.Slice(slices, func(i, j int) bool {
		return slices[i].Rest < slices[j].Rest
	})
}

func Decompose(M [][]float64, TargetDiagVec []int) (D []map[int][]float64, err error) {
	d := len(M)
	var U map[int][]float64
	U, err = Matrix2DiagonalVectors(M)
	if err != nil {
		return nil, err
	}

	S := make([]*NonZeroEntry, d)
	i := 0
	for diagVec, components := range U {
		for idx, value := range components {
			if value == 1 {
				S[i] = NewnonZeroEntry(diagVec, idx, d, TargetDiagVec)
				i++
			}
		}
	}
	SortByRestChoice(S)

	LastOccupied := make(map[Point]Point)
	bruteforceTime := 1
	for _, s := range S {
		LastOccupied[Point{X: s.DiagVecCoordinate.X, Y: s.DiagVecCoordinate.Y}] = Point{X: -1, Y: -1}
		bruteforceTime *= s.Rest
	}
	fmt.Printf("BruteForce Times:%d", bruteforceTime)

	Stack := make([]*NonZeroEntry, 0)
	Stack = append(Stack, S[0])
	S = S[1:]
	for len(Stack) < d {
		s := Stack[len(Stack)-1]
		// clean the previous effect:
		LastOccupiedRow := LastOccupied[Point{X: s.DiagVecCoordinate.X, Y: s.DiagVecCoordinate.Y}].X
		LastOccupiedCol := LastOccupied[Point{X: s.DiagVecCoordinate.X, Y: s.DiagVecCoordinate.Y}].Y
		if LastOccupiedRow != -1 && LastOccupiedCol != -1 {
			// Not first time choosing, recover its effect.
			for _, a := range S {
				a.Recover(LastOccupiedRow, LastOccupiedCol, d)
			}
		}

		// Rechoose new path
		occupiedRow, occupiedCol, ok := s.Choose(d)

		if ok {
			// Update infos:
			LastOccupied[Point{X: s.DiagVecCoordinate.X, Y: s.DiagVecCoordinate.Y}] = Point{X: occupiedRow, Y: occupiedCol}

			var safe bool
			for _, a := range S {
				_, safe = a.Forbid(occupiedRow, occupiedCol, d)
				/*
					if !safe {
						a.Recover(occupiedRow, occupiedCol, d)
						break
					}
				*/
			}
			if safe {
				SortByRestChoice(S)
				Stack = append(Stack, S[0])
				S = S[1:]
			}
			/* version 1
			var aChosed *NonZeroEntry
			var aChosedidx int
			var firstChosed = false
			var effected, safe bool
			for idx, a := range S {
				effected, safe = a.Forbid(occupiedRow, occupiedCol, d)
				if !safe {
					a.Recover(occupiedRow, occupiedCol, d)
					break
				}
				if effected && !firstChosed {
					aChosed = a
					aChosedidx = idx
					firstChosed = true
				}
			}
			if !firstChosed && safe {
				Stack = append(Stack, S[0])
				S = S[1:]
			} else if safe {
				Stack = append(Stack, aChosed)
				S = append(S[:aChosedidx], S[aChosedidx+1:]...)
			}
			*/

		} else { // pop it out.
			LastOccupied[Point{X: s.DiagVecCoordinate.X, Y: s.DiagVecCoordinate.Y}] = Point{X: -1, Y: -1}

			if len(Stack) == 0 {
				break
			}
			idx := len(Stack) - 1
			Stack = Stack[:idx]
			s.ReInit(d)
			S = append(S, s)
			/*
				Stemp := make([]*NonZeroEntry, len(S))
				copy(Stemp, S) // FIXED memory sharing issue. Can not use ":=" to do the copy.
				Stemp[0] = s
				S = append(Stemp[:1], S...)
			*/
		}
		fmt.Printf("last index of Stack: %d, Point (%d,%d)\n", len(Stack), Stack[len(Stack)-1].NormalCoordinate.X, Stack[len(Stack)-1].NormalCoordinate.Y)
	}
	s := Stack[len(Stack)-1]
	s.Choose(d)

	D = make([]map[int][]float64, 2)
	D[0] = make(map[int][]float64)
	D[1] = make(map[int][]float64)
	for _, i := range TargetDiagVec {
		D[0][i] = make([]float64, d)
		D[1][i] = make([]float64, d)
	}
	for _, entry := range Stack {
		diagVecpair, componentpair := entry.LatestChoice(d)
		D[0][diagVecpair.X][componentpair.X] = 1
		D[1][diagVecpair.Y][componentpair.Y] = 1
	}
	return

}

// return true if two matrix matches.
// return false otherwise
func MatrixCompare(A [][]float64, B [][]float64) (equal bool) {
	equal = true
	if len(A) != len(B) {
		equal = false
		return
	} else if len(A[0]) != len(B[0]) {
		equal = false
		return
	}
	for i := 0; i < len(A); i++ {
		for j := 0; j < len(A[0]); j++ {
			if A[i][j] != B[i][j] {
				equal = false
				return
			}
		}
	}
	return
}

func DeepIntFloat64MapCopy(src map[int][]float64) (dst map[int][]float64) {
	dst = make(map[int][]float64)
	for key, value := range src {
		copiedSlice := make([]float64, len(value))
		copy(copiedSlice, value)
		dst[key] = copiedSlice
	}
	return
}
