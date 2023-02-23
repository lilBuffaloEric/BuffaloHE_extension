package check

/*
 * Check for funcs from matrix_mult, and more advanced routines like compute Covariant Matrix (for PCA)
 */

import (
	"fmt"
	"math"
	auxio "project1-fhe_extension/auxiliary_io"
	mtrxmult "project1-fhe_extension/matrix_mult"
	"runtime"
	"sync"
	"time"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func SquareMatrix_product_check() {
	d := 5
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
	// ------------------------------------------------------ testing
	var M map[int][]float64
	M, err = mtrxmult.Gen_trans_C_tao_diagonalVectors(d)
	if err != nil {
		panic(err)
	}
	print(M)

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
	var now time.Time
	// create a 4x4 matrix matrixA
	d := 128
	BSGSRatio := 4
	MatrixA := make([][]float64, d)
	for i := 0; i < d; i++ {
		MatrixA[i] = make([]float64, d)
		for j := 0; j < d; j++ {
			MatrixA[i][j] = float64(i*d + j)
		}
	}
	auxio.Print_matrix_f64_2d(MatrixA, d, d)
	// encode matrixA into a vector using row ordering
	var a []float64
	a, err = mtrxmult.Row_orderingInv(MatrixA)
	if err != nil {
		panic(err)
	}
	// initialize encryption scheme.
	LogN := 15
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

	// Check Sigma permutation (with BSGS)
	fmt.Printf("Checking Sigma permutation (with BSGS)...\n")
	var ctSigmaA, ctSigmaA_BSGS *rlwe.Ciphertext
	now = time.Now()
	ctSigmaA, err = mtrxmult.Sigma_linearTransform(params, rlk, galk, ctA, d) // NonBSGS version
	fmt.Printf("NonBSGS version Done in (%s) \n", time.Since(now))
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctSigmaA, d, d)
	now = time.Now()
	ctSigmaA_BSGS, err = mtrxmult.Sigma_linearTransformBSGS(params, rlk, galk, ctA, d, float64(BSGSRatio)) //BSGS version
	if err != nil {
		panic(err)
	}
	fmt.Printf("BSGS version with Ration %d Done in (%s) \n", BSGSRatio, time.Since(now))
	auxio.Quick_check_matrix(params, sk, ctSigmaA_BSGS, d, d)
	fmt.Printf("Compare with exact Sigma Permutation:\n")
	checkA := mtrxmult.Sigma_permute(MatrixA, d)
	auxio.Print_matrix_f64_2d(checkA, d, d)

	// Check Tao permutation
	fmt.Printf("Checking Tao permutation (with BSGS) ...\n")
	var ctTaoA, ctTaoA_BSGS *rlwe.Ciphertext
	now = time.Now()
	ctTaoA, err = mtrxmult.Tao_linearTransform(params, rlk, galk, ctA, d) // nonBSGS version
	fmt.Printf("NonBSGS version Done in (%s) \n", time.Since(now))
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctTaoA, d, d)
	now = time.Now()
	ctTaoA_BSGS, err = mtrxmult.Tao_linearTransformBSGS(params, rlk, galk, ctA, d, float64(BSGSRatio)) // BSGS version
	fmt.Printf("BSGS version with Ration %d Done in (%s) \n", BSGSRatio, time.Since(now))
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctTaoA_BSGS, d, d)
	fmt.Printf("Compare with exact Tao Permutation:\n")
	checkA = mtrxmult.Tao_permute(MatrixA, d)
	auxio.Print_matrix_f64_2d(checkA, d, d)

	// Check ColShift permutation
	fmt.Printf("Checking ColShift permutation ...\n")
	k := 0
	var ctColShiftA *rlwe.Ciphertext
	ctColShiftA, err = mtrxmult.ColShift_linearTransform(params, rlk, galk, ctA, d, k)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctColShiftA, d, d)
	checkA = mtrxmult.Phi_permute(MatrixA, d)
	auxio.Print_matrix_f64_2d(checkA, d, d)

	// Check RowShift permutation
	fmt.Printf("Checking RowShift permutation ...\n")
	var ctRowShiftA *rlwe.Ciphertext
	ctRowShiftA, err = mtrxmult.RowShift_linearTransform(params, rlk, galk, ctA, d, k)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctRowShiftA, d, d)
	checkA = mtrxmult.Psi_permute(MatrixA, d)
	auxio.Print_matrix_f64_2d(checkA, d, d)

	// Check Square Matrix Multiplication:
	fmt.Printf("Checking Square Matrix Multiplication (with BSGS)...\n")
	var ctMatrixMult, ctMatrixMult_BSGS *rlwe.Ciphertext
	ctB := ctA
	now = time.Now()
	ctMatrixMult, err = mtrxmult.SquareMatrix_Mult(params, rlk, galk, ctA, ctB, d) // nonBSGS version
	if err != nil {
		panic(err)
	}
	fmt.Printf("nonBSGS version SquareMatrix Mult done in (%s)\n", time.Since(now))
	auxio.Quick_check_matrix(params, sk, ctMatrixMult, d, d)
	now = time.Now()
	ctMatrixMult_BSGS, err = mtrxmult.SquareMatrix_MultBSGS_dbg(params, sk, rlk, galk, ctA, ctB, d, float64(BSGSRatio), float64(BSGSRatio)) // BSGS version
	if err != nil {
		panic(err)
	}
	fmt.Printf("BSGS version with SigmaBSGSRatio %d , TaoBSGSRatio %d done in (%s)\n", BSGSRatio, BSGSRatio, time.Since(now))
	auxio.Quick_check_matrix(params, sk, ctMatrixMult_BSGS, d, d)
	fmt.Printf("Compare with exact result of SquareMatrix Mult:\n")
	checkA, err = mtrxmult.SquareMatrix_product_permute_version(MatrixA, MatrixA)
	if err != nil {
		panic(err)
	}
	auxio.Print_matrix_f64_2d(checkA, d, d)

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
	auxio.Quick_check_matrix(params, sk, ctC, d, d)
	var ctMaskedC *rlwe.Ciphertext
	ctMaskedC, err = mtrxmult.Matrix_masking_withColOrder(params, rlk, ctC, d, d, l, 1)
	if err != nil {
		panic(err)
	}
	evaluator.Rescale(ctMaskedC, params.DefaultScale(), ctMaskedC)
	auxio.Quick_check_matrix(params, sk, ctMaskedC, d, d)

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
	ctReplicatedMaskedC, err = mtrxmult.Matrix_rowTotalSum_withColOrder(params, galk_forReplicate, ctMaskedC, d, d)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctReplicatedMaskedC, d, d)

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
	auxio.Quick_check_matrix(params, sk, ctD, d, d)
	var ctReplicatedD *rlwe.Ciphertext
	ctReplicatedD, err = mtrxmult.Matrix_multirowReplicate_withColOrder_dbg(params, sk, galk_forReplicate, ctD, d, d, l)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctReplicatedD, d, d)

	// Check Complete Row rotation in Col ordering scheme.
	var ctRotatedD *rlwe.Ciphertext
	for i := 0; i < d; i++ {
		ctRotatedD, err = mtrxmult.Matrix_rowRotation_withColOrder_dbg(params, sk, galk_forReplicate, ctReplicatedD, d, d, i)
		if err != nil {
			panic(err)
		}
		auxio.Quick_check_matrix(params, sk, ctRotatedD, d, d)
	}

	// Check Rectangle Matrix Multiplication:
	var ctCD *rlwe.Ciphertext
	ctCD, err = mtrxmult.RectangleMatrix_Mult_withColOrder_dbg(params, sk, galk_forReplicate, rlk, ctC, ctD, d, d, m, l, n)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctCD, d, d)

	// Check (Square) Matrix Transposition:
	var ctTransA *rlwe.Ciphertext
	ctTransA, err = mtrxmult.Transpose_linearTransform(params, rlk, galk_forReplicate, ctA, d)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctTransA, d, d)

	// Check (Square) Matrix Transposition + Tao permutation:
	var ctTransCtaoA *rlwe.Ciphertext
	ctTransCtaoA, err = mtrxmult.TransposeCtao_linearTransform(params, rlk, galk_forReplicate, ctA, d)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctTransCtaoA, d, d)
	// compare with the result of doing transpostion and tao permutation seperately:
	var ctSeperateA *rlwe.Ciphertext
	ctSeperateA, err = mtrxmult.Tao_linearTransform(params, rlk, galk_forReplicate, ctTransA, d)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctSeperateA, d, d)

}

func OnlyMatrixMult_check() {
	var err error
	var now time.Time
	// create a 4x4 matrix matrixA
	d := 4
	SigmaBSGSRatio := 4
	TaoBSGSRatio := 4
	MatrixA := make([][]float64, d)
	for i := 0; i < d; i++ {
		MatrixA[i] = make([]float64, d)
		for j := 0; j < d; j++ {
			MatrixA[i][j] = float64(i*d + j)
		}
	}
	auxio.Print_matrix_f64_2d(MatrixA, d, d)
	// encode matrixA into a vector using row ordering
	var a []float64
	a, err = mtrxmult.Row_orderingInv(MatrixA)
	if err != nil {
		panic(err)
	}
	var checkA [][]float64
	// initialize encryption scheme.
	LogN := 15
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

	// Check Square Matrix Multiplication:
	fmt.Printf("Checking Square Matrix Multiplication (with BSGS)...\n")
	var ctMatrixMult, ctMatrixMult_BSGS *rlwe.Ciphertext
	ctB := ctA
	now = time.Now()
	ctMatrixMult, err = mtrxmult.SquareMatrix_Mult(params, rlk, galk, ctA, ctB, d) // nonBSGS version
	if err != nil {
		panic(err)
	}
	fmt.Printf("nonBSGS version SquareMatrix Mult done in (%s)\n", time.Since(now))
	auxio.Quick_check_matrix(params, sk, ctMatrixMult, d, d)
	now = time.Now()
	ctMatrixMult_BSGS, err = mtrxmult.SquareMatrix_MultBSGS_dbg(params, sk, rlk, galk, ctA, ctB, d, float64(SigmaBSGSRatio), float64(TaoBSGSRatio)) // BSGS version
	if err != nil {
		panic(err)
	}
	fmt.Printf("BSGS version with SigmaBSGSRatio %d , TaoBSGSRatio %d done in (%s)\n", SigmaBSGSRatio, TaoBSGSRatio, time.Since(now))
	auxio.Quick_check_matrix(params, sk, ctMatrixMult_BSGS, d, d)
	fmt.Printf("Compare with exact result of SquareMatrix Mult:\n")
	checkA, err = mtrxmult.SquareMatrix_product_permute_version(MatrixA, MatrixA)
	if err != nil {
		panic(err)
	}
	auxio.Print_matrix_f64_2d(checkA, d, d)
}

func PCA_check() {
	var err error
	matrixCols := 112
	matrixRows := 128
	rows := 60000
	cols := 785
	d := 1 << 7
	maxSubRoutes := 4
	runtime.GOMAXPROCS(3)

	// initialize encryption scheme.
	LogN := 15
	Q := []uint64{
		1152921504598720513,
		1099490000897, 1099498258433, 1099499175937, 1099499569153, 1099500617729}

	P := []uint64{1152921504606584833}
	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:         LogN,
		Q:            Q,
		P:            P,
		LogSlots:     int(math.Log2(float64(d * d))),
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
	steps := make([]int, d+d-1+d-1+2*d-2)
	j := 0
	for i := 0; i < d; i++ {
		steps[j] = i
		j++
		if i != 0 {
			steps[j] = -i + d*d // -i is original but we will have to set -i mod d^2 here.
			j++
			steps[j] = (d - 1) * i
			j++
			steps[j] = (d-1)*(-i) + d*d
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

	// evaluator
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})

	// Open the file
	// var csvfile *os.File
	csvfilepath := "mnist_train.csv"

	// "D:\\cardistry\\discovery@GY\\软件学院 2th\\机器学习基础\\机器学习实验三\\mnist\\mnist_train.csv"

	/*
		csvfile, err = os.Open(csvfilepath)
		if err != nil {
			log.Fatalln("Couldn't open the csv file", err)
		}
	*/
	// defer csvfile.Close()

	// Construct cols-1/matrixCols x cols-1/matrixCols ciphertexts:
	Cov := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols)
	for i := 0; i < (cols-1)/matrixCols; i++ {
		Cov[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	}
	// Construct [rows/d] matrices with size d x d
	ORSB1 := make([][][]float64, int(math.Ceil(float64(rows)/float64(d))))
	// ORSB2 := make([][][]float64, int(math.Ceil(float64(rows)/float64(d))))
	for i := 0; i < int(math.Ceil(float64(rows)/float64(d))); i++ {
		ORSB1[i] = make([][]float64, d)
		// ORSB2[i] = make([][]float64, d)
		for j := 0; j < d; j++ {
			ORSB1[i][j] = make([]float64, d)
			// ORSB2[i][j] = make([]float64, d)
		}
	}

	// Construct ((cols-1)/matrixCols)^2 Ciphertexts representing a cols x cols Covariance Matrix
	ctCov := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctMean := make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctMMT := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols)
	for i := 0; i < (cols-1)/matrixCols; i++ {
		ctCov[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
		ctMMT[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	}
	// We need some more temp Ciphertext to help...
	var ctTemp *rlwe.Ciphertext
	ctSum := encryptor.EncryptNew(auxio.Encode_single_int(params, 0, params.MaxLevel()))
	ctORSB := make([]*rlwe.Ciphertext, int(math.Ceil(float64(rows)/float64(matrixRows))))
	ctSubRbuff := make([]*rlwe.Ciphertext, maxSubRoutes)

	// Need to Store some of the ciphertexts into disk.
	var Covtracker = auxio.Filetracker{Fp: "ctCov", Off: 0, Hierarchy: make([]int, 2)}
	var Meantracker = auxio.Filetracker{Fp: "ctMean", Off: 0, Hierarchy: make([]int, 1)}
	var Datatracker = auxio.Filetracker{Fp: "ctData", Off: 0, Hierarchy: make([]int, 2)}
	Covtracker.Hierarchy[1] = (cols - 1) / matrixCols
	Datatracker.Hierarchy[1] = int(math.Ceil(float64(rows) / float64(matrixRows)))
	var tempbytes = 0
	// Start Computing X*X^T and Mean
	var wg sync.WaitGroup
	for k := 1; k < cols; k += matrixCols {
		for j := 1; j < rows+1; j += matrixRows {
			// extract a matrixRows x matrixCols matrix

			/* SingleThread Mode
			var sb [][]string
			var sbf64 [][]float64
			var MatrixSB []float64
			var ptMtrxSB *rlwe.Plaintext

			if j+matrixRows > rows {
				sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, k, rows-j+1, matrixCols)
			} else {
				sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, k, matrixRows, matrixCols)
			}
			sbf64, err = auxio.Switch2d_str2f64(sb)
			if err != nil {
				panic(err)
			}
			MatrixSB, err = mtrxmult.Row_orderingInvZeroPad(sbf64, d)
			ptMtrxSB = encoder.EncodeNew(MatrixSB, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
			ctORSB[j/matrixRows] = encryptor.EncryptNew(ptMtrxSB)
			copy(ORSB1[j/matrixRows], sbf64)
			*/

			// MultiThread Mode
			wg.Add(1)
			go func(d int, ORSB1 [][][]float64, j int) {
				defer wg.Done()
				var sb [][]string
				var sbf64 [][]float64
				var MatrixSB []float64
				var ptMtrxSB *rlwe.Plaintext
				encryptor_subroute := ckks.NewEncryptor(params, pk)
				if j+matrixRows > rows+1 {
					sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, k, rows-j+1, matrixCols)
				} else {
					sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, k, matrixRows, matrixCols)
				}
				sbf64, err = auxio.Switch2d_str2f64(sb)
				if err != nil {
					panic(err)
				}
				MatrixSB, err = mtrxmult.Row_orderingInvZeroPad(sbf64, d)
				ptMtrxSB = encoder.EncodeNew(MatrixSB, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
				ctORSB[j/matrixRows] = encryptor_subroute.EncryptNew(ptMtrxSB)
				copy(ORSB1[j/matrixRows], sbf64)
				fmt.Printf("the %d th routine done\n", j)
				// return
			}(d, ORSB1, j)
		}
		wg.Wait()

		for i := 1; i < cols; i += matrixCols {
			// somehow we need to constraint the number of subroutes...
			subroutesNum := maxSubRoutes
			for j := 1; j < rows+1; j += matrixRows {
				// extract a matrixRows x matrixCols matrix

				/* Single Thread Mode
				var sb [][]string
				var sbf64 [][]float64
				var MatrixSB []float64
				var ptMtrxSB *rlwe.Plaintext
				var ctMtrxSB *rlwe.Ciphertext

				if j+matrixRows > rows {
					sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, i, rows-j+1, matrixCols)
				} else {
					sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, i, matrixRows, matrixCols)
				}
				sbf64, err = auxio.Switch2d_str2f64(sb)
				if err != nil {
					panic(err)
				}
				MatrixSB, err = mtrxmult.Row_orderingInvZeroPad(sbf64, d)
				ptMtrxSB = encoder.EncodeNew(MatrixSB, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
				ctMtrxSB = encryptor.EncryptNew(ptMtrxSB)
				ctMtrxSB, err = mtrxmult.Transpose_linearTransform(params, rlk, galk, ctMtrxSB, d)
				ctTemp, err = mtrxmult.SquareMatrix_Mult(params, rlk, galk, ctORSB[j/matrixRows], ctMtrxSB, d)
				err = evaluator.Rescale(ctTemp, ctSum.Scale, ctTemp)
				ctSum = evaluator.AddNew(ctTemp, ctSum)
				*/

				// MultiThread Mode
				wg.Add(1)
				subroutesNum--
				go func(d int, j int, sbrNum int) {
					defer wg.Done()
					var sb [][]string
					var sbf64 [][]float64
					var MatrixSB []float64
					var ptMtrxSB *rlwe.Plaintext
					encryptor_subroute := ckks.NewEncryptor(params, pk)
					// evaluator_subroute := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
					if j+matrixRows > rows+1 {
						sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, k, rows-j+1, matrixCols)
					} else {
						sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, k, matrixRows, matrixCols)
					}
					sbf64, err = auxio.Switch2d_str2f64(sb)
					if err != nil {
						panic(err)
					}
					MatrixSB, err = mtrxmult.Row_orderingInvZeroPad(sbf64, d)
					ptMtrxSB = encoder.EncodeNew(MatrixSB, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
					ctSubRbuff[sbrNum] = encryptor_subroute.EncryptNew(ptMtrxSB)
					ctSubRbuff[sbrNum], err = mtrxmult.Transpose_linearTransform(params, rlk, galk, ctSubRbuff[sbrNum], d)
					ctTemp, err = mtrxmult.SquareMatrix_Mult(params, rlk, galk, ctORSB[j/matrixRows], ctSubRbuff[sbrNum], d)
					// err = evaluator_subroute.Rescale(ctTemp, ctSum.Scale, ctTemp)
					ctSubRbuff[sbrNum] = ctTemp
					fmt.Printf("the %d th iteration done, this is subroutine No.%d\n", j, sbrNum)
					// return
				}(d, j, subroutesNum)
				if subroutesNum <= 0 {
					wg.Wait()
					subroutesNum = maxSubRoutes
					for j := 0; j < maxSubRoutes; j++ {
						err = evaluator.Rescale(ctSubRbuff[j], ctSum.Scale, ctSubRbuff[j])
						ctSum = evaluator.AddNew(ctSubRbuff[j], ctSum)
					}
				}
			}
			wg.Wait()
			for j := subroutesNum; j < maxSubRoutes; j++ {
				err = evaluator.Rescale(ctSubRbuff[j], ctSum.Scale, ctSubRbuff[j])
				ctSum = evaluator.AddNew(ctSubRbuff[j], ctSum)
			}
			evaluator.MultByConst(ctSum, 1/rows, ctSum)
			/*
				for j := 1; j < rows; j += matrixRows {
					err = evaluator.Rescale(ctTemp, ctSum.Scale, ctTemp)
					ctSum = evaluator.AddNew(ctSubRbuff[j/matrixRows], ctSum)
				}
			*/
			fmt.Printf("one subMatrix Inner product Done \n")
			// ctCov[(i-1)/matrixCols][(k-1)/matrixCols] = ctSum // might delete when inner store is on
			// Inner store:
			tempbytes, err = auxio.WriteCt2file(Covtracker.Fp, ctSum, Covtracker.Off)
			Covtracker.Off += tempbytes
			Covtracker.Hierarchy[0] = ctSum.MarshalBinarySize()
			ctSum = encryptor.EncryptZeroNew(params.MaxLevel())
		}
		fmt.Printf("X*X^T %d/7 process Done\n", (k-1)/matrixCols+1)
		// At the same time we compute 1/7 of the Mean vector.
		ctMean[(k-1)/matrixCols] = encryptor.EncryptNew(auxio.Encode_single_int(params, 0, params.MaxLevel()))
		for j := 1; j < rows+1; j += matrixRows {
			evaluator.Add(ctMean[(k-1)/matrixCols], ctORSB[(j-1)/matrixRows], ctMean[(k-1)/matrixCols])
		}
		ctMean[(k-1)/matrixCols], err = mtrxmult.Matrix_rowTotalSum_withColOrder(params, galk, ctMean[(k-1)/matrixCols], d, d)
		if err != nil {
			panic(err)
		}
		evaluator.MultByConst(ctMean[(k-1)/matrixCols], 1/rows, ctMean[(k-1)/matrixCols])
		fmt.Printf("Mean %d/7 process Done\n", (k-1)/matrixCols+1)
		// We want to store the ctORSB,ctCov,ctMean into disk, to reduce the memory occupation...
		// Store the mean:
		tempbytes, err = auxio.WriteCt2file(Meantracker.Fp, ctMean[(k-1)/matrixCols], Meantracker.Off)
		Meantracker.Off += tempbytes
		Meantracker.Hierarchy[0] = ctMean[(k-1)/matrixCols].MarshalBinarySize()
		ctMean[(k-1)/matrixCols] = nil // release the memory.

		/* Store the Covariance, this may be done in the inner loop instead
		tempbytes,err = auxio.WriteCtVec2file(Covtracker.Fp, ctCov[(k-1)/matrixCols],Covtracker.Off)
		Covtracker.Off += tempbytes
		Covtracker.Hierarchy[0] = tempbytes
		*/

		// Store the Data:
		tempbytes, err = auxio.WriteCtVec2file(Datatracker.Fp, ctORSB, Datatracker.Off)
		Datatracker.Off += tempbytes
		Datatracker.Hierarchy[0] = ctORSB[0].MarshalBinarySize()

	}
	// Release the ctORSB's memory, since all of the subblock ciphertexts are stored in disk.
	for i := 0; i < len(ctORSB); i++ {
		ctORSB[i] = nil
	}

	// Compute Mean*Mean^T, first retrieving the Mean from disk.
	_, err = auxio.ReadCtVec4file(Meantracker.Fp, Meantracker.Hierarchy[0], ctMean, 0)
	for i := 0; i < len(ctMean); i++ {
		ctTemp, err = mtrxmult.Transpose_linearTransform(params, rlk, galk, ctMean[i], d)
		for j := 0; j < len(ctMean); j++ {
			ctMMT[j][i] = evaluator.MulRelinNew(ctTemp, ctMean[j])
		}
	}

	// Retrieve the X*X^T from disk, and do Cov = X*X^T-Mean*Mean^T
	for i := 0; i < len(ctCov); i++ {
		_, err = auxio.ReadCtVec4file(Covtracker.Fp, Covtracker.Hierarchy[0], ctCov[i], i*Covtracker.Hierarchy[0]*Covtracker.Hierarchy[1])
		for j := 0; j < len(ctCov); j++ {
			evaluator.Sub(ctCov[i][j], ctMMT[i][j], ctCov[i][j])
			auxio.Quick_check_matrix(params, sk, ctCov[i][j], d, d)
		}
	}

}

/*
func PCA_check() {
	var err error

	var btp *bootstrapping.Bootstrapper
	var kgen rlwe.KeyGenerator
	var encoder ckks.Encoder
	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var encryptor rlwe.Encryptor
	var decryptor rlwe.Decryptor

	// Bootstrapping parameters
	// Two sets of four parameters each, DefaultParametersSparse and DefaultParametersDense,
	// (each index 0 to 3) ensuring 128 bit of security are available in
	// github.com/tuneinsight/lattigo/v4/ckks/bootstrapping/default_params.go
	//
	// LogSlots is hardcoded to 15 in the parameters, but can be changed from 1 to 15.
	// When changing LogSlots make sure that the number of levels allocated to CtS and StC is
	// smaller or equal to LogSlots.

	paramSet := bootstrapping.DefaultParametersSparse[1] // bootstrapping.DefaultParametersDense[1]
	ckksParams := paramSet.SchemeParams
	// ckksParams.LogN = 15
	ckksParams.LogSlots = 14
	btpParams := paramSet.BootstrappingParams

	params, err := ckks.NewParametersFromLiteral(ckksParams)
	if err != nil {
		panic(err)
	}

	fmt.Println()
	fmt.Printf("CKKS parameters: logN = %d, logSlots = %d, H(%d; %d), logQP = %d, levels = %d, scale= 2^%f, sigma = %f \n", params.LogN(), params.LogSlots(), params.HammingWeight(), btpParams.EphemeralSecretWeight, params.LogQP(), params.QCount(), math.Log2(params.DefaultScale().Float64()), params.Sigma())

	// Scheme context and keys
	kgen = ckks.NewKeyGenerator(params)

	sk, pk = kgen.GenKeyPair()
	// Relinearization key
	rlk := kgen.GenRelinearizationKey(sk, 1)
	// Galois keys
	// Since we are testing Sigma, Tao, RowShift and ColShift permutation,
	// rotation steps {-d+1 ~ d-1} , {1d,2d,...(d-1)d} and {-d(d-1),...,(d-1)(d-1)} are needed.
	d := 1 << (params.LogSlots() / 2)
	matrixRows := 112
	matrixCols := 128
	steps := make([]int, d+d-1+d-1+d-1+d-3)
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
	for i := -d; i < -1; i++ {
		steps[j] = i*(d-1) + d*d
		j++
	}
	for i := 2; i < (d - 1); i++ {
		steps[j] = i * (d - 1)
		j++
	}
	galk := kgen.GenRotationKeysForRotations(steps[:], false, sk)

	encoder = ckks.NewEncoder(params)
	decryptor = ckks.NewDecryptor(params, sk)
	encryptor = ckks.NewEncryptor(params, pk)

	fmt.Println()
	fmt.Println("Generating bootstrapping keys...")
	evk := bootstrapping.GenEvaluationKeys(btpParams, params, sk)
	fmt.Println("Done")

	if btp, err = bootstrapping.NewBootstrapper(params, btpParams, evk); err != nil {
		panic(err)
	}

	// Open the file
	var csvfile *os.File
	csvfilepath := ""
	csvfile, err = os.Open(csvfilepath)
	if err != nil {
		log.Fatalln("Couldn't open the csv file", err)
	}
	defer csvfile.Close()
	rows := 60000
	cols := 785

	// Parse the file
	reader := csv.NewReader(csvfile)
	// Construct cols-1/matrixCols x cols-1/matrixCols ciphertexts:
	Cov := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols)
	for i := 0; i < (cols-1)/matrixCols; i++ {
		Cov[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	}
	// Construct [rows/d] matrices with size d x d
	ORSB1 := make([][][]string, int(math.Ceil(float64(rows)/float64(d))))
	for i := d; i < int(math.Ceil(float64(rows)/float64(d))); i++ {
		ORSB1[i] = make([][]string, d)
		for j := d; j < d; j++ {
			ORSB1[i][j] = make([]string, d)
		}
	}
	// Skip the first row:
	_, err = reader.Read()
	var line []string
	// Iterate through the records to get
	for i := 1; i < cols; i += matrixCols {
		for j := 0; j < rows; j++ {
			line, err = reader.Read()
			for k := i; k < i+matrixCols; k++ {
				ORSB1[j/matrixRows][j%matrixRows][k%matrixCols] = line[k]
			}
		}
		fmt.Printf("One Row of Sub blocks completed")
	}

}
*/
