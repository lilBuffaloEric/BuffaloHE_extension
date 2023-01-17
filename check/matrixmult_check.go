package check

import (
	"math"
	auxio "project1-fhe_extension/auxiliary_io"
	mtrxmult "project1-fhe_extension/matrix_mult"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func SquareMatrix_product_check() {
	d := 3
	MatrixA := make([][]float64, d)
	for i := 0; i < d; i++ {
		MatrixA[i] = make([]float64, d)
		for j := 0; j < d; j++ {
			MatrixA[i][j] = float64(i*d + j)
		}
	}

	MatrixB := MatrixA
	auxio.Print_matrix_f64_full_2d(MatrixA, d, d)
	MatrixC, err := mtrxmult.SquareMatrix_product_permute_version(MatrixA, MatrixB)
	if err != nil {
		panic(err)
	}
	auxio.Print_matrix_f64_full_2d(MatrixC, d, d)

}

func SquareMatrix_encode_check() {
	d := 3
	g := 2
	Multi_MatrixA := make([][][]float64, g)
	for k := 0; k < g; k++ {
		Multi_MatrixA[k] = make([][]float64, d)
		for i := 0; i < d; i++ {
			Multi_MatrixA[k][i] = make([]float64, d)
			for j := 0; j < d; j++ {
				Multi_MatrixA[k][i][j] = float64(i*d + j)
			}
		}
	}
	a, err := mtrxmult.Row_orderingInv(Multi_MatrixA[0])
	if err != nil {
		panic(err)
	}
	auxio.Print_vector_f64_full(a)
	a, err = mtrxmult.Row_orderingInv_multiple(Multi_MatrixA)
	if err != nil {
		panic(err)
	}
	auxio.Print_vector_f64_full(a)
}

func SquareMatrix_permute_check() {
	var err error
	// create a 4x4 matrix matrixA
	d := 4
	MatrixA := make([][]float64, d)
	for i := 0; i < d; i++ {
		MatrixA[i] = make([]float64, d)
		for j := 0; j < d; j++ {
			MatrixA[i][j] = float64(i*d + j)
		}
	}
	auxio.Print_matrix_f64_full_2d(MatrixA, d, d)
	// encode matrixA into a vector using row ordering
	var a []float64
	a, err = mtrxmult.Row_orderingInv(MatrixA)
	if err != nil {
		panic(err)
	}
	// initialize encryption scheme.
	LogN := 14
	Q := []uint64{
		1152921504598720513,
		1099490000897, 1099498258433, 1099499175937, 1099499569153, 1099500617729}

	P := []uint64{1152921504606584833}
	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:         LogN,
		Q:            Q,
		P:            P,
		LogSlots:     int(math.Log2(float64(d * d))), // we only need 16 slots to represent a 4x4 matrix in a ciphertext.
		DefaultScale: 1 << 40,
	}); err != nil {
		panic(err)
	}
	// generate Keys:
	kgen := ckks.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPair()
	// Relinearization key
	rlk := kgen.GenRelinearizationKey(sk, 1)
	// Galois keys
	// Since we are testing Sigma, Tao, RowShift and ColShift permutation,
	// rotation steps {-d+1 ~ d-1} and  {1d,2d,...(d-1)d} are needed.
	steps := make([]int, d+d-1+d-1)
	j := 0
	for i := 0; i < d; i++ {
		steps[j] = i
		j++
		if i != 0 {
			steps[j] = -i + d*d // -i is original but we will have to set -i mod d^2 here.
			j++
		}
	}
	for i := d; i < d*d; i += d {
		steps[j] = i
		j++
	}
	galk := kgen.GenRotationKeysForRotations(steps[:], false, sk)

	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)

	// encoder
	encoder := ckks.NewEncoder(params)

	// encryption
	ptA := encoder.EncodeNew(a, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	ctA := encryptor.EncryptNew(ptA)

	// Check Sigma permutation
	var ctSigmaA *rlwe.Ciphertext
	ctSigmaA, err = mtrxmult.Sigma_linearTransform(params, rlk, galk, ctA, d)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix_full(params, sk, ctSigmaA, d, d)
	checkA := mtrxmult.Sigma_permute(MatrixA, d)
	auxio.Print_matrix_f64_full_2d(checkA, d, d)

	// Check Tao permutation
	var ctTaoA *rlwe.Ciphertext
	ctTaoA, err = mtrxmult.Tao_linearTransform(params, rlk, galk, ctA, d)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix_full(params, sk, ctTaoA, d, d)
	checkA = mtrxmult.Tao_permute(MatrixA, d)
	auxio.Print_matrix_f64_full_2d(checkA, d, d)

	// Check ColShift permutation
	k := 0
	var ctColShiftA *rlwe.Ciphertext
	ctColShiftA, err = mtrxmult.ColShift_linearTransform(params, rlk, galk, ctA, d, k)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix_full(params, sk, ctColShiftA, d, d)
	checkA = mtrxmult.Phi_permute(MatrixA, d)
	auxio.Print_matrix_f64_full_2d(checkA, d, d)

	// Check RowShift permutation
	var ctRowShiftA *rlwe.Ciphertext
	ctRowShiftA, err = mtrxmult.RowShift_linearTransform(params, rlk, galk, ctTaoA, d, k)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix_full(params, sk, ctRowShiftA, d, d)
	checkA = mtrxmult.Psi_permute(MatrixA, d)
	auxio.Print_matrix_f64_full_2d(checkA, d, d)

	// Check Matrix Multiplication:
	var ctMatrixMult *rlwe.Ciphertext
	ctB := ctA
	ctMatrixMult, err = mtrxmult.SquareMatrix_Mult_dbg(params, sk, rlk, galk, ctA, ctB, d)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix_full(params, sk, ctMatrixMult, d, d)
	checkA, err = mtrxmult.SquareMatrix_product_permute_version(MatrixA, MatrixA)
	if err != nil {
		panic(err)
	}
	auxio.Print_matrix_f64_full_2d(checkA, d, d)

	// Check Matrix masking:
	// Regard Plaintext Slots as a 4x4 square matrix
	// set a rectangle matrix C with size (m=3) x (l=2) as a submatrix of Plaintext square matrix.
	f := 1
	m := 3
	l := 2
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	MatrixC := make([][]float64, d)
	for i := 0; i < d; i++ {
		MatrixC[i] = make([]float64, d)
	}
	for i := 0; i < m; i++ {
		for j := 0; j < l; j++ {
			MatrixC[i][j] = float64(f)
			f++
		}
	}
	var C []float64
	C, err = mtrxmult.Col_orderingInv(MatrixC)
	if err != nil {
		panic(err)
	}
	ptC := encoder.EncodeNew(C, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	ctC := encryptor.EncryptNew(ptC)
	auxio.Quick_check_matrix_full(params, sk, ctC, d, d)
	var ctMaskedC *rlwe.Ciphertext
	ctMaskedC, err = mtrxmult.Matrix_masking_withColOrder(params, rlk, ctC, d, d, l, 1)
	if err != nil {
		panic(err)
	}
	evaluator.Rescale(ctMaskedC, params.DefaultScale(), ctMaskedC)
	auxio.Quick_check_matrix_full(params, sk, ctMaskedC, d, d)

	// Check Matrix Replication:
	// For simplicity, We generate all the possible Rotation keys for this operation.
	steps_forReplicate := make([]int, 2*d*d)
	for i := 0; i < d*d; i++ {
		steps_forReplicate[i] = i
		steps_forReplicate[i+d*d] = -i
	}
	galk_forReplicate := kgen.GenRotationKeysForRotations(steps_forReplicate, false, sk)
	// Use ctMaskedC to display the effect.
	// Row Replication:
	var ctReplicatedMaskedC *rlwe.Ciphertext
	ctReplicatedMaskedC, err = mtrxmult.Matrix_rowTotalSum(params, galk_forReplicate, ctMaskedC, d, d)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix_full(params, sk, ctReplicatedMaskedC, d, d)

	// MultiRow Replication:
	// set another rectangle matrix D with size (l=2) x (n=4) as a submatrix of Plaintext square matrix.
	f = 1
	n := 4
	MatrixD := make([][]float64, d)
	for i := 0; i < d; i++ {
		MatrixD[i] = make([]float64, d)
	}
	for i := 0; i < l; i++ {
		for j := 0; j < n; j++ {
			MatrixD[i][j] = float64(f)
			f++
		}
	}
	var D []float64
	D, err = mtrxmult.Col_orderingInv(MatrixD)
	if err != nil {
		panic(err)
	}
	ptD := encoder.EncodeNew(D, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	ctD := encryptor.EncryptNew(ptD)
	auxio.Quick_check_matrix_full(params, sk, ctD, d, d)
	var ctReplicatedD *rlwe.Ciphertext
	ctReplicatedD, err = mtrxmult.Matrix_multirowReplicate_withColOrder_dbg(params, sk, galk_forReplicate, ctD, d, d, l)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix_full(params, sk, ctReplicatedD, d, d)

	// Check Complete Row rotation in Col ordering scheme.
	var ctRotatedD *rlwe.Ciphertext
	for i := 0; i < d; i++ {
		ctRotatedD, err = mtrxmult.Matrix_rowRotation_withColOrder_dbg(params, sk, galk_forReplicate, ctReplicatedD, d, d, i)
		if err != nil {
			panic(err)
		}
		auxio.Quick_check_matrix_full(params, sk, ctRotatedD, d, d)
	}

	// Check Rectangle Matrix Multiplication:
	var ctCD *rlwe.Ciphertext
	ctCD, err = mtrxmult.RectangleMatrix_Mult_withColOrder_dbg(params, sk, galk_forReplicate, rlk, ctC, ctD, d, d, m, l, n)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix_full(params, sk, ctCD, d, d)

}
