package check

/*
 * Check for funcs from matrix_mult, and more advanced routines like compute Covariant Matrix (for PCA)
 */

import (
	"fmt"
	"math"
	"runtime/debug"
	"sync"
	"time"

	auxio "project1-fhe_extension_v1.0/auxiliary_io"
	mtrxmult "project1-fhe_extension_v1.0/matrix_mult"
	nonpolyfunc "project1-fhe_extension_v1.0/nonpoly_func"

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
	Sigma_BSGSRatio := 8
	DSigma_BSGSRatio := 1
	Tau_BSGSRatio := 8
	DTau_BSGSRatio := 1
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

	// create the Map of Sigma LinearTransformation
	var SigmaDiagonalMap map[int][]float64
	SigmaDiagonalMap, err = mtrxmult.Gen_sigma_diagonalVecotrs(d)
	if err != nil {
		panic(err)
	}
	// create the Map of Tau LinearTransformation

	var TauDiagonalMap map[int][]float64
	TauDiagonalMap, err = mtrxmult.Gen_tao_diagonalVectors(d)
	if err != nil {
		panic(err)
	}

	// create the Decomposed Maps of Sigma LinearTransformation
	SigmaMatrix := mtrxmult.DiagonalVectors2Matrix(SigmaDiagonalMap, d*d)
	DSigmaDiagonalMaps, err := mtrxmult.Converge2DiagonalDecompose_Sigma(SigmaMatrix)
	if err != nil {
		panic(err)
	}
	// create the Decomposed Maps of Tau LinearTransformation

	TauMatrix := mtrxmult.DiagonalVectors2Matrix(TauDiagonalMap, d*d)
	DTauDiagonalMaps, err := mtrxmult.Converge2DiagonalDecompose_Tao(TauMatrix)
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

	// Check Sigma permutation (with BSGS and Decomposed BSGS version)
	fmt.Printf("Checking Sigma permutation (with BSGS)...\n")
	var ctSigmaA, ctSigmaA_BSGS, ctSigmaA_DBSGS *rlwe.Ciphertext

	now = time.Now()
	ctSigmaA, err = mtrxmult.Sigma_linearTransform_dbg(params, rlk, galk, ctA, d) // NonBSGS version
	fmt.Printf("NonBSGS version Done in (%s) \n", time.Since(now))
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctSigmaA, d, d)

	now = time.Now()
	ctSigmaA_BSGS, err = mtrxmult.General_linearTransform_dbg(params, rlk, galk, ctA, SigmaDiagonalMap, float64(Sigma_BSGSRatio)) //BSGS version
	if err != nil {
		panic(err)
	}
	fmt.Printf("BSGS version with Ratio %d Done in (%s) \n", Sigma_BSGSRatio, time.Since(now))
	auxio.Quick_check_matrix(params, sk, ctSigmaA_BSGS, d, d)

	now = time.Now()
	ctSigmaA_DBSGS, err = mtrxmult.General_linearTransform_dbg(params, rlk, galk, ctA, DSigmaDiagonalMaps[0], float64(DSigma_BSGSRatio))
	ctSigmaA_DBSGS, err = mtrxmult.General_linearTransform_dbg(params, rlk, galk, ctSigmaA_DBSGS, DSigmaDiagonalMaps[1], float64(DSigma_BSGSRatio))
	fmt.Printf("Decomposed BSGS version with Ratio %d Done in (%s) \n ", Sigma_BSGSRatio, time.Since(now))
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctSigmaA_DBSGS, d, d)

	fmt.Printf("Compare with exact Sigma Permutation:\n")
	checkA := mtrxmult.Sigma_permute(MatrixA, d)
	auxio.Print_matrix_f64_2d(checkA, d, d)

	// Check Tao permutation
	fmt.Printf("Checking Tao permutation (with BSGS) ...\n")
	var ctTaoA, ctTaoA_BSGS, ctTauA_DBSGS *rlwe.Ciphertext

	now = time.Now()
	ctTaoA, err = mtrxmult.Tao_linearTransform(params, rlk, galk, ctA, d) // nonBSGS version
	fmt.Printf("NonBSGS version Done in (%s) \n", time.Since(now))
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctTaoA, d, d)

	now = time.Now()
	ctTaoA_BSGS, err = mtrxmult.General_linearTransform_dbg(params, rlk, galk, ctA, TauDiagonalMap, float64(Tau_BSGSRatio)) // BSGS version
	fmt.Printf("BSGS version with Ration %d Done in (%s) \n", Tau_BSGSRatio, time.Since(now))
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctTaoA_BSGS, d, d)

	now = time.Now()
	ctTauA_DBSGS, err = mtrxmult.General_linearTransform_dbg(params, rlk, galk, ctA, DTauDiagonalMaps[0], float64(DTau_BSGSRatio))
	ctTauA_DBSGS, err = mtrxmult.General_linearTransform_dbg(params, rlk, galk, ctTauA_DBSGS, DTauDiagonalMaps[1], float64(DTau_BSGSRatio))
	fmt.Printf("Decomposed BSGS version with Ratio %d Done in (%s) \n", Tau_BSGSRatio, time.Since(now))
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix(params, sk, ctTauA_DBSGS, d, d)

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
	var ctMatrixMult, ctMatrixMult_BSGS, ctMatrixMult_DBSGS *rlwe.Ciphertext
	ctB := ctA
	now = time.Now()
	ctMatrixMult, err = mtrxmult.SquareMatrix_Mult(params, rlk, galk, ctA, ctB, d) // nonBSGS version
	if err != nil {
		panic(err)
	}
	fmt.Printf("nonBSGS version SquareMatrix Mult done in (%s)\n", time.Since(now))
	auxio.Quick_check_matrix(params, sk, ctMatrixMult, d, d)

	now = time.Now()
	ctMatrixMult_BSGS, err = mtrxmult.SquareMatrix_MultBSGS_dbg(params, sk, rlk, galk, ctA, ctB, d, float64(Sigma_BSGSRatio), float64(Tau_BSGSRatio)) // BSGS version
	if err != nil {
		panic(err)
	}
	fmt.Printf("BSGS version with SigmaBSGSRatio %d , TaoBSGSRatio %d done in (%s)\n", Sigma_BSGSRatio, Tau_BSGSRatio, time.Since(now))
	auxio.Quick_check_matrix(params, sk, ctMatrixMult_BSGS, d, d)

	now = time.Now()
	ctMatrixMult_DBSGS, err = mtrxmult.SquareMatrix_MultBSGS_DecomposeVer_dbg(params, sk, rlk, galk, ctA, ctB, DSigmaDiagonalMaps, DTauDiagonalMaps, float64(Sigma_BSGSRatio), float64(Tau_BSGSRatio))
	if err != nil {
		panic(err)
	}
	fmt.Printf("Decomposed BSGS version with SigmaBSGSRatio %d, TauBSGSRatio %d done in (%s)\n", Sigma_BSGSRatio, Tau_BSGSRatio, time.Since(now))
	auxio.Quick_check_matrix(params, sk, ctMatrixMult_DBSGS, d, d)

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

func PCA_check(loadkeysFromDisk bool) {
	var err error
	matrixCols := 128
	matrixRows := 128
	rows := 60000
	cols := 257

	// BSGSRatio:
	/*
		// BGSGRatio4Transpose := 2
		BSGSRatio4Sigma := 2
		BSGSRatio4Tau := 2
	*/
	// Maximum Routines:
	maxSubRoutes := 7

	// N1
	DSigma_BSGSN1 := 16
	DTau_BSGSN1 := 16
	N1Trans := 16

	// runtime.GOMAXPROCS(6)
	d := 1 << 7

	// MaxDiagVec
	SigmaMaxDiagVec := d / 4
	TauMaxDiagVec := (d / 4) * d

	// necessary file paths:
	var KeysPath = "Keys"
	var Rtk4Transpose_path = "Rtk4Transpose"
	var Rtk4MtrxMult_path = "Rtk4MtrxMult"
	var Rtk4RowtotalSum_path = "Rtk4RowtotalSum"
	var ctCov_path = "ctCov"
	var ctMean_path = "ctMean"
	var csvfilepath = "mnist.csv"

	// initialize encryption scheme params.
	var params ckks.Parameters
	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var rlk *rlwe.RelinearizationKey
	var kgen rlwe.KeyGenerator

	// generate Keys, we can load keys from disk or generate it.
	if loadkeysFromDisk {
		params, sk, pk, rlk, err = LoadKeys_check()
		if err != nil {
			panic(err)
		}
	} else {
		params, err = ckks.NewParametersFromLiteral(PN15QP720S14)
		if err != nil {
			panic(err)
		}
		kgen = ckks.NewKeyGenerator(params)
		sk, pk = kgen.GenKeyPair()
		// Relinearization key
		rlk = kgen.GenRelinearizationKey(sk, 1)

		// Print Keys' size:
		fmt.Printf("Encryption Scheme Params occupy %d bytes\n", params.MarshalBinarySize())
		fmt.Printf("Secret Key:%d bytes, Public Key:%d bytes\n", sk.MarshalBinarySize(), pk.MarshalBinarySize())
		fmt.Printf("Relinearization Key: %d bytes\n", rlk.MarshalBinarySize())

		// Store this set of keys, for convenience:
		keysTracker := auxio.NewTracker(KeysPath, -1)
		n := 0
		n, err = keysTracker.StoreUpdateOne(&params)
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store params %d bytes\n", n)
		}
		n, err = keysTracker.StoreUpdateOne(sk)
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store sk %d bytes\n", n)
		}
		n, err = keysTracker.StoreUpdateOne(pk)
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store pk %d bytes\n", n)
		}
		n, err = keysTracker.StoreUpdateOne(rlk)
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store rlk %d bytes\n", n)
		}
		n, err = keysTracker.StoreFinish()
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store filetracker %d bytes\n", n)
		}
	}

	inputLevel := 9

	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)

	// Decryptor
	// decryptor := ckks.NewDecryptor(params, sk)

	// encoder
	encoder := ckks.NewEncoder(params)

	// evaluator
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})

	// create the Map and LinearTransform object of Transpose LinearTransformation
	fmt.Printf("create the Map and LinearTransform object of Transpose LinearTransformation\n")
	var TransposeDiagonalMap map[int][]float64
	TransposeDiagonalMap, err = mtrxmult.Gen_transpose_diagonalVectors(d)
	if err != nil {
		panic(err)
	}
	TransposeLT := ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, TransposeDiagonalMap, params.MaxLevel(), params.DefaultScale(), N1Trans, (d - 1), params.LogSlots())

	// create the Map of Sigma LinearTransformation
	fmt.Printf("create the Map of Sigma LinearTransformation\n")
	DSigmaDiagonalMaps1, DSigmaDiagonalMaps2, err := mtrxmult.GenSigmaDiagnalDecomposeMatrices(d, SigmaMaxDiagVec)
	if err != nil {
		panic(err)
	}

	// create the Map of Tau LinearTransformation
	fmt.Printf("create the Map of Tau LinearTransformation\n")
	DTauDiagonalMaps, err := mtrxmult.GenTauDiagonalDecomposeMatrices(d, TauMaxDiagVec)
	if err != nil {
		panic(err)
	}

	// create the Decomposed Maps and LinearTransform object of Sigma LinearTransformation
	fmt.Printf("create the Decomposed Maps and LinearTransform object of Sigma LinearTransformation\n")
	DSigmaLTs1 := make([]ckks.LinearTransform, len(DSigmaDiagonalMaps1))
	DSigmaLTs2 := make([]ckks.LinearTransform, len(DSigmaDiagonalMaps2))
	for i := range DSigmaLTs1 {
		if i == len(DSigmaLTs1)-1 {
			DSigmaLTs1[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DSigmaDiagonalMaps1[i], params.MaxLevel(), params.DefaultScale(), DSigma_BSGSN1, 1, params.LogSlots())
		} else {
			DSigmaLTs1[i] = ckks.GenLinearTransform(encoder, DSigmaDiagonalMaps1[i], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		}
	}
	for i := range DSigmaLTs2 {
		if i == len(DSigmaLTs2)-1 {
			DSigmaLTs2[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DSigmaDiagonalMaps2[i], params.MaxLevel(), params.DefaultScale(), DSigma_BSGSN1, 1, params.LogSlots())
		} else {
			DSigmaLTs2[i] = ckks.GenLinearTransform(encoder, DSigmaDiagonalMaps2[i], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		}
	}

	// create the Decomposed Maps of Tau LinearTransformation
	fmt.Printf("create the Decomposed Maps of Tau LinearTransformation\n")
	DTauLTs := make([]ckks.LinearTransform, len(DTauDiagonalMaps))
	for i := range DTauLTs {
		if i == len(DTauLTs)-1 {
			DTauLTs[i] = ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, DTauDiagonalMaps[i], params.MaxLevel(), params.DefaultScale(), DTau_BSGSN1, d, params.LogSlots())
		} else {
			DTauLTs[i] = ckks.GenLinearTransform(encoder, DTauDiagonalMaps[i], params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		}
	}

	// create the Map of Plaintext for ColShift Diagonals
	fmt.Printf("create the Map of Plaintext for ColShift Diagonals\n")
	ColShiftLTs := make([]ckks.LinearTransform, 0)
	for k := 1; k < d; k++ {
		U, err := mtrxmult.Gen_colShift_diagonalVectors(d, k)
		if err != nil {
			panic(err)
		}
		LT := ckks.GenLinearTransform(encoder, U, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
		ColShiftLTs = append(ColShiftLTs, LT)
	}

	// Create Rotation keys correspond to Transposition, MatrixMult and RowTotalSum
	var galk4Transpose *rlwe.RotationKeySet
	var galk4MtrxMult *rlwe.RotationKeySet
	var galk4RowtotalSum *rlwe.RotationKeySet

	// Create File Tracker to Store or load these keys
	fmt.Printf("create File Tracker to Store or Load the Rotation keys\n")
	var RtkTracker4Transpose *auxio.Filetracker
	var RtkTracker4MtrxMult *auxio.Filetracker
	var RtkTracker4RowtotalSum *auxio.Filetracker

	if loadkeysFromDisk {
		// we are not going to imediately get these keys from disk when they already exists in disk.
		fmt.Printf("we only Init the Tracker for rotationKeySets, but are not going to imediately get these keys from disk.\n")
		RtkTracker4Transpose = auxio.NewTracker4File(Rtk4Transpose_path)
		RtkTracker4MtrxMult = auxio.NewTracker4File(Rtk4MtrxMult_path)
		RtkTracker4RowtotalSum = auxio.NewTracker4File(Rtk4RowtotalSum_path)
	} else {
		// Create File Tracker to Store or load these keys
		fmt.Printf("Init File Tracker to Store or Load the Rotation keys\n")
		RtkTracker4Transpose = auxio.NewTracker(Rtk4Transpose_path, -1)
		RtkTracker4MtrxMult = auxio.NewTracker(Rtk4MtrxMult_path, -1)
		RtkTracker4RowtotalSum = auxio.NewTracker(Rtk4RowtotalSum_path, -1)

		// kgen = ckks.NewKeyGenerator(params)

		// Create Rotation keys correspond to Transposition, MatrixMult and RowTotalSum
		fmt.Printf("create Rotation keys correspond to Transposition, MatrixMult and RowTotalSum\n")
		fmt.Printf("creating Rotation keys for Transposition\n")
		Rotation4Transpose := TransposeLT.Rotations4ArithmeticSeq()
		galk4Transpose = kgen.GenRotationKeysForRotations(Rotation4Transpose, false, sk)

		fmt.Printf("creating Rotation keys for MatrixMult\n")
		var MtrxMultRotmap = make(map[int]bool)
		for k := range DSigmaLTs1 {
			for _, i := range DSigmaLTs1[k].Rotations4ArithmeticSeq() {
				if !MtrxMultRotmap[i] {
					MtrxMultRotmap[i] = true
				}
			}
		}
		for k := range DSigmaLTs2 {
			for _, i := range DSigmaLTs2[k].Rotations4ArithmeticSeq() {
				if !MtrxMultRotmap[i] {
					MtrxMultRotmap[i] = true
				}
			}
		}
		for k := range DTauLTs {
			for _, i := range DTauLTs[k].Rotations4ArithmeticSeq() {
				if !MtrxMultRotmap[i] {
					MtrxMultRotmap[i] = true
				}
			}
		}
		var Rotation4MtrxMult = make([]int, 0)
		for rot := range MtrxMultRotmap {
			Rotation4MtrxMult = append(Rotation4MtrxMult, rot)
		}
		Rotation4MtrxMult = append(Rotation4MtrxMult, -d+d*d)
		galk4MtrxMult = kgen.GenRotationKeysForRotations(Rotation4MtrxMult, false, sk)

		fmt.Printf("creating Rotation keys for RowtotalSum\n")
		Rotation4RowTotalSum := make([]int, d)
		for i := d; i < d*d; i = (i << 1) {
			Rotation4RowTotalSum = append(Rotation4RowTotalSum, i)
		}
		// Rotation4RowTotalSum = append(Rotation4RowTotalSum, -d+d*d)
		galk4RowtotalSum = kgen.GenRotationKeysForRotations(Rotation4RowTotalSum, false, sk)

		// Store the rotation keys if they are freshly generated.
		fmt.Printf("Store the rotation keys since they are freshly generated.\n")
		_, err = RtkTracker4Transpose.StoreUpdateOne(galk4Transpose)
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store galk4Transpose %d bytes\n", galk4Transpose.MarshalBinarySize())
		}
		_, err = RtkTracker4Transpose.StoreFinish()
		if err != nil {
			panic(err)
		}

		_, err = RtkTracker4MtrxMult.StoreUpdateOne(galk4MtrxMult)
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store galk4MtrxMult %d bytes\n", galk4MtrxMult.MarshalBinarySize())
		}
		_, err = RtkTracker4MtrxMult.StoreFinish()
		if err != nil {
			panic(err)
		}

		_, err = RtkTracker4RowtotalSum.StoreUpdateOne(galk4RowtotalSum)
		if err != nil {
			panic(err)
		} else {
			fmt.Printf("Store galk4RowtotalSum %d bytes\n", galk4RowtotalSum.MarshalBinarySize())
		}
		_, err = RtkTracker4RowtotalSum.StoreFinish()
		if err != nil {
			panic(err)
		}
	}

	// Free the memory:
	galk4Transpose = nil
	galk4MtrxMult = nil
	galk4RowtotalSum = nil
	debug.FreeOSMemory()

	// Construct cols-1/matrixCols x cols-1/matrixCols ciphertexts:
	Cov := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols)
	for i := 0; i < (cols-1)/matrixCols; i++ {
		Cov[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	}
	// Construct [rows/d] matrices with size d x d
	ORSB1 := make([][][]float64, int(math.Ceil(float64(rows)/float64(d))))
	for i := 0; i < int(math.Ceil(float64(rows)/float64(d))); i++ {
		ORSB1[i] = make([][]float64, d)
		for j := 0; j < d; j++ {
			ORSB1[i][j] = make([]float64, d)
		}
	}

	// Construct ((cols-1)/matrixCols)^2 Ciphertexts representing a cols x cols Covariance Matrix
	ctCov := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctMean := make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctMMT := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols) // ciphertext for MeanVec · MeanVec^T
	for i := 0; i < (cols-1)/matrixCols; i++ {
		ctCov[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
		ctMMT[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	}
	// We need some more temp Ciphertext to help...
	// var ctTemp *rlwe.Ciphertext
	// var ctSum *rlwe.Ciphertext
	// ctSum := encryptor.EncryptNew(auxio.Encode_single_float64(params, float64(0), params.MaxLevel(), rlwe.NewScale(1)))
	// ctSum := encryptor.EncryptZeroNew(params.MaxLevel())                                   // ciphertext for inner product of two Rows of Subblock
	ctORSBT := make([]*rlwe.Ciphertext, int(math.Ceil(float64(rows)/float64(matrixRows))))                  // ciphetexts for One Row of Subblock Transposed
	ctORSB := make([]*rlwe.Ciphertext, int(math.Ceil(float64(rows)/float64(matrixRows))))                   // ciphetexts for One Row of Subblock.
	var ctORSBT_SigmaDecomp = make([][]*rlwe.Ciphertext, int(math.Ceil(float64(rows)/float64(matrixRows)))) // Sigma Decomposed ciphetexts for One Row of Subblock Transposed.
	// ctSubRbuff := make([]*rlwe.Ciphertext, maxSubRoutes)                                   // subroutine buffer.
	// ctOCSB := make([]*rlwe.Ciphertext, int(math.Ceil(float64(cols)/float64(matrixCols))))
	var ctOCSBSum []*rlwe.Ciphertext
	// ctOCSBSum := make([]*rlwe.Ciphertext, int(math.Ceil(float64(cols)/float64(matrixCols))))
	IdentityVec := make([]float64, 1<<params.LogSlots())
	for i := 0; i < len(IdentityVec); i++ {
		IdentityVec[i] = 0.5
	}
	// ptIdentity := encoder.EncodeNew(IdentityVec, params.MaxLevel(), params.DefaultScale(), params.LogSlots())

	// Need to Store some of the ciphertexts into disk.
	var Covtracker = auxio.NewTracker(ctCov_path, -1)
	var Meantracker = auxio.NewTracker(ctMean_path, -1)
	// var Datatracker = auxio.NewTracker("ctData", -1)

	// Start Computing X*X^T and Mean
	fmt.Printf("Start Computing X*X^T and Mean\n")
	var elapsed time.Duration
	var now time.Time
	var wg sync.WaitGroup
	for k := 1; k < cols; k += matrixCols {
		// Load the Rotation keys for Transposition
		RtkTracker4Transpose = auxio.NewTracker4File(Rtk4Transpose_path)
		galk4Transpose = new(rlwe.RotationKeySet)
		_, err = RtkTracker4Transpose.ReadUpdateOne(galk4Transpose)
		if err != nil {
			panic(err)
		}
		fmt.Printf("Retrieve the Rotation keys with %d MB for Transposition\n", galk4Transpose.MarshalBinarySize()/(1024*1024))

		// Generate evaluation key for transposition and other operation:
		elk := rlwe.EvaluationKey{Rlk: rlk, Rtks: galk4Transpose}
		// Extract the Transposition of OneRowOfSubblock.
		now = time.Now()
		fmt.Printf("Extract the Transposition of %d th OneRowOfSubblock, and by the way compute their Sigma transformation result in all column shift\n", (k-1)/matrixCols)
		subrouteNum := maxSubRoutes
		for j := 1; j < rows+1; j += matrixRows {
			// extract a matrixRows x matrixCols matrix
			// MultiThread Mode
			wg.Add(1)
			subrouteNum--
			go func(d int, ORSB1 [][][]float64, j int) {
				now := time.Now()
				defer wg.Done()
				var sb [][]string            // One subblock in string type
				var sbf64 [][]float64        // One subblick in float64 type
				var MatrixSB []float64       // Row odering Vector for one subblock
				var ptMtrxSB *rlwe.Plaintext // Plaintext for one subblock
				var ctMtrxSB *rlwe.Ciphertext
				encoder_subroute := ckks.NewEncoder(params)
				encryptor_subroute := ckks.NewEncryptor(params, pk)
				evaluator_subroute := ckks.NewEvaluator(params, elk)
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
				ptMtrxSB = encoder_subroute.EncodeNew(MatrixSB, inputLevel, params.DefaultScale(), params.LogSlots())
				ctMtrxSB = encryptor_subroute.EncryptNew(ptMtrxSB)
				// auxio.Quick_decode_matrix_full(params, ptMtrxSB, d, d)
				ctORSB[j/matrixRows] = ctMtrxSB                                                                       // dbg testing
				ctORSBT[j/matrixRows] = evaluator_subroute.LinearTransform4ArithmeticSeqNew(ctMtrxSB, TransposeLT)[0] // encrypt the Transposed Plaintext for one subblock into the OneRowOfSubblock ciphertext list.
				if err != nil {
					panic(err)
				}
				err = evaluator_subroute.Rescale(ctORSBT[j/matrixRows], params.DefaultScale(), ctORSBT[j/matrixRows])
				if err != nil {
					panic(err)
				}
				// copy(ORSB1[j/matrixRows], sbf64)
				fmt.Printf("the %d th routine done in %s \n", j, time.Since(now))
				fmt.Printf("ctORSBT at level %d, Scale %f\n", ctORSBT[j/matrixRows].Level(), math.Log2(ctORSBT[j/matrixRows].Scale.Float64()))
				// auxio.Quick_check_matrix_full(params, sk, ctORSBT[j/matrixRows], d, d)
				// return
			}(d, ORSB1, j)
			if subrouteNum <= 0 {
				wg.Wait()
				subrouteNum = maxSubRoutes
			}
		}
		wg.Wait()
		fmt.Printf("Extractiong & Transposition complete in %s with one ciphertext: %d bytes, take up space %d MB in total\n", time.Since(now), ctORSBT[0].MarshalBinarySize(), ctORSBT[0].MarshalBinarySize()*len(ctORSBT)*2/(1024*1024)) // dbg testing
		elapsed += time.Since(now)
		// Free the Rotation keys for Transposition
		galk4Transpose = nil

		// We then compute the Sum of One Row of Subblocks to prepare for the Aggregation of Mean vector computation
		now = time.Now()
		fmt.Printf("We then compute the Sum of One Row of Subblocks to prepare for the Aggregation of Mean vector computation\n")
		ctMean[(k-1)/matrixCols] = encryptor.EncryptZeroNew(params.MaxLevel()) // dbg testing
		for j := 1; j < rows+1; j += matrixRows {
			evaluator.Add(ctMean[(k-1)/matrixCols], ctORSB[(j-1)/matrixRows], ctMean[(k-1)/matrixCols])
			ctORSB[(j-1)/matrixRows] = nil // free the space.
		}
		fmt.Printf("Sum of Subblocks complete in %s\n", time.Since(now))
		elapsed += time.Since(now)
		debug.FreeOSMemory()

		// Retrieve the Rotation keys for MatrixMult
		RtkTracker4MtrxMult = auxio.NewTracker4File(Rtk4MtrxMult_path)
		galk4MtrxMult = new(rlwe.RotationKeySet)
		_, err = RtkTracker4MtrxMult.ReadUpdateOne(galk4MtrxMult)
		if err != nil {
			panic(err)
		}
		fmt.Printf("Retrieve the Rotation keys with %d MB for MatrixMult", galk4MtrxMult.MarshalBinarySize()/(1024*1024))

		// Compute the Sigma Decomposition of the ciphertexts in ctORSBT:
		//now = time.Now()
		//fmt.Printf("Compute the Sigma Decomposition of the ciphertexts in ctORSBT\n")
		// ctORSBT_SigmaDecomp, err = mtrxmult.Sigma_linearTransform_MultiThread_dbg(params, sk, rlk, galk4MtrxMult, ctORSBT, DSigmaLTs1, DSigmaLTs2, ColShiftLTs, d, len(ctORSBT))
		//fmt.Printf("Sigma Decomposition of the ciphertexts in ctORSBT Done in %s\n", time.Since(now))
		//elapsed += time.Since(now)

		/* dbg testing
		for i := 0; i < len(ctORSBT); i++ {
			auxio.Quick_check_matrix(params, sk, ctORSBT[i], d, d)
		}
		*/

		// Start Computing Matrix Mult for each element in ctORSBT, each of them will be multiplied by elements in ctOCSB in paralle.
		now = time.Now()
		fmt.Printf("Start Computing Matrix Mult for each element in ctORSBT, each of them will be multiplied by elements in ctOCSB in paralle.\n")
		ctOCSBSum = make([]*rlwe.Ciphertext, (k-1)/matrixCols+1)

		subrouteNum_outer := maxSubRoutes / ((k-1)/matrixCols + 1)
		cond := sync.NewCond(&sync.Mutex{})
		j := 1
		for idx := range ctORSBT_SigmaDecomp {
			if ctORSBT_SigmaDecomp[idx] == nil {
				wg.Wait()
				nowtmp := time.Now()
				var ctORSBT_SigmaDecompTmp [][]*rlwe.Ciphertext
				if idx+maxSubRoutes > len(ctORSBT) {
					fmt.Printf("Computing SigmaLinTrans for ctORSBT[%d:%d]\n", idx, len(ctORSBT))
					ctORSBT_SigmaDecompTmp, err = mtrxmult.Sigma_linearTransform_MultiThread_dbg(params, sk, rlk, galk4MtrxMult, ctORSBT[idx:], DSigmaLTs1, DSigmaLTs2, ColShiftLTs, d, maxSubRoutes)
				} else {
					fmt.Printf("Computing SigmaLinTrans for ctORSBT[%d:%d]\n", idx, idx+maxSubRoutes)
					ctORSBT_SigmaDecompTmp, err = mtrxmult.Sigma_linearTransform_MultiThread_dbg(params, sk, rlk, galk4MtrxMult, ctORSBT[idx:idx+maxSubRoutes], DSigmaLTs1, DSigmaLTs2, ColShiftLTs, d, maxSubRoutes)
				}
				if err != nil {
					panic(err)
				}
				fmt.Printf("ctORSBT_SigmaDecompTmp occupying %d MB\n", ctORSBT_SigmaDecompTmp[0][0].MarshalBinarySize()*len(ctORSBT_SigmaDecompTmp)*len(ctORSBT_SigmaDecompTmp[0])/(1024*1024))
				copy(ctORSBT_SigmaDecomp[idx:idx+len(ctORSBT_SigmaDecompTmp)], ctORSBT_SigmaDecompTmp)
				fmt.Printf("SigmaLinTrans complete in %s\n", time.Since(nowtmp))
				for i := 0; i < idx; i++ {
					ctORSBT_SigmaDecomp[i] = nil
				}
				nowtmp = time.Now()
				debug.FreeOSMemory()
				elapsed -= time.Since(nowtmp)
			}
			wg.Add(1)
			subrouteNum_outer--
			go func(idx int, j int) {
				defer wg.Done()
				ctOCSB := make([]*rlwe.Ciphertext, (k-1)/matrixCols+1)
				var wg_subroutine sync.WaitGroup
				// Retrieve OCSB of one row
				subrouteNum_inner := ((k-1)/matrixCols + 1)
				for i := 1; i <= k; i += matrixCols {
					fmt.Printf("Retrieve dataset[%d][%d], and load it into ctOCSB\n", idx, (i-1)/matrixCols)
					wg_subroutine.Add(1)
					subrouteNum_inner--
					go func(d int, i int, j int) {
						defer wg_subroutine.Done()
						var sb [][]string
						var sbf64 [][]float64
						var MatrixSB []float64
						var ptMtrxSB *rlwe.Plaintext
						encryptor_subroute := ckks.NewEncryptor(params, pk)
						encoder_subroute := ckks.NewEncoder(params)
						// evaluator_subroute := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
						if j+matrixRows > rows+1 {
							sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, i, rows-j+1, matrixCols)
						} else {
							sb, err = auxio.GetOneCSVSubblock(csvfilepath, j, i, matrixRows, matrixCols)
						}
						sbf64, err = auxio.Switch2d_str2f64(sb)
						if err != nil {
							panic(err)
						}
						MatrixSB, err = mtrxmult.Row_orderingInvZeroPad(sbf64, d)
						ptMtrxSB = encoder_subroute.EncodeNew(MatrixSB, inputLevel, params.DefaultScale(), params.LogSlots())
						ctOCSB[(i-1)/matrixCols] = encryptor_subroute.EncryptNew(ptMtrxSB)
					}(d, i, j)
					if subrouteNum_inner <= 0 {
						wg_subroutine.Wait()
						subrouteNum_inner = ((k-1)/matrixCols + 1)
					}
				}
				wg_subroutine.Wait()
				// Compute MtrxMult between element in ctORSBT and elements in ctOCSB:
				fmt.Printf("Compute MtrxMult between ctORSBT[%d] and elements in %d-th ctOCSB, and add them to the previous result in ctOCSBSum:\n", idx, idx)
				ctOCSBTemp, err := mtrxmult.SquareMatrix_MultBSGS_DecomposeVer3_NeedLTs_dbg(params, sk, rlk, galk4MtrxMult, ctORSBT_SigmaDecomp[idx], ctOCSB, DTauLTs, d, len(ctOCSB))
				// ctOCSBTemp, err := mtrxmult.SquareMatrix_MultBSGS_DecomposeVer2_NeedLTs_dbg(params, sk, rlk, galk4MtrxMult, ct, ctOCSB, DSigmaLTs1, DSigmaLTs2, DTauLTs, ColShiftLTs, d, len(ctOCSB))
				if err != nil {
					panic(err)
				}
				cond.L.Lock()
				if ctOCSBSum[0] == nil {
					fmt.Printf("Let ctOCSBTemp[%d] as the first value of ctOCBSum\n", idx)
					for i := range ctOCSBTemp {
						ctOCSBSum[i] = ctOCSBTemp[i].CopyNew()
						//fmt.Printf("First ctOCSBSum[%d] value:\n", i)
						//auxio.Quick_check_matrix(params, sk, ctOCSBSum[i], d, d)
					}
				} else {
					fmt.Printf("Adding ctOCSBTemp[%d] into ctOCBSum\n", idx)
					for i := range ctOCSBTemp {
						//fmt.Printf("New value that will be added to ctOCSBSum[%d]:\n", i)
						//auxio.Quick_check_matrix(params, sk, ctOCSBTemp[i], d, d)
						evaluator.Add(ctOCSBTemp[i], ctOCSBSum[i], ctOCSBSum[i])
						//fmt.Printf("After adding a new one, ctOCSBSum[%d] value:\n", i)
						//auxio.Quick_check_matrix(params, sk, ctOCSBSum[i], d, d)
					}
				}
				cond.L.Unlock()
			}(idx, j)
			if subrouteNum_outer <= 0 {
				wg.Wait()
				subrouteNum_outer = maxSubRoutes / ((k-1)/matrixCols + 1)
			}
			j += matrixRows
		}
		wg.Wait()
		fmt.Printf("Scale the result in ctOCSBSum, and load them into ctCov:\n")
		for i := range ctOCSBSum {
			err = evaluator.Rescale(ctOCSBSum[i], params.DefaultScale(), ctOCSBSum[i])
			if err != nil {
				panic(err)
			}
			//fmt.Printf("ctCov[%d][%d] before scaling has scale %f, level %d\n", (k-1)/matrixCols, i, math.Log2(ctOCSBSum[i].Scale.Float64()), ctOCSBSum[i].Level())
			// auxio.Quick_check_matrix(params, sk, ctOCSBSum[i], d, d)
			// ctOCSBSum[i] = encryptor.EncryptNew(encoder.EncodeNew(auxio.DecryptDecode(params, sk, ctOCSBSum[i]), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
			evaluator.MultByConst(ctOCSBSum[i], 1.0/float64(rows), ctOCSBSum[i])
			err = evaluator.Rescale(ctOCSBSum[i], params.DefaultScale(), ctOCSBSum[i])
			if err != nil {
				panic(err)
			}
			fmt.Printf("ctCov[%d][%d] now has scale %f, level %d\n", (k-1)/matrixCols, i, math.Log2(ctOCSBSum[i].Scale.Float64()), ctOCSBSum[i].Level())
			ctCov[(k-1)/matrixCols][i] = ctOCSBSum[i].CopyNew()
			auxio.Quick_check_matrix(params, sk, ctCov[(k-1)/matrixCols][i], d, d)
		}

		galk4MtrxMult = nil
		fmt.Printf("X*X^T %d/7 process Done in %s\n", (k-1)/matrixCols+1, time.Since(now))
		elapsed += time.Since(now)
		debug.FreeOSMemory()

		// At the same time we computie 1/7 of the Mean vector by Aggregating the OneRowOfSubblock ciphertexts and RowTotalSum.
		now = time.Now()
		fmt.Printf("At the same time we computie 1/7 of the Mean vector by Aggregating the OneRowOfSubblock ciphertexts and RowTotalSum.\n")
		// First Retrieve the RotationKeys for RowtotalSum :
		RtkTracker4RowtotalSum = auxio.NewTracker4File(Rtk4RowtotalSum_path)
		galk4RowtotalSum = new(rlwe.RotationKeySet)
		_, err = RtkTracker4RowtotalSum.ReadUpdateOne(galk4RowtotalSum)
		if err != nil {
			panic(err)
		}
		ctMean[(k-1)/matrixCols], err = mtrxmult.Matrix_Aggregate_withRowOrder(params, sk, galk4RowtotalSum, ctMean[(k-1)/matrixCols], d, d, 0)
		if err != nil {
			panic(err)
		}
		// auxio.Quick_check_matrix_full(params, sk, ctMean[(k-1)/matrixCols], d, d) // dbg testing
		evaluator.MultByConst(ctMean[(k-1)/matrixCols], 1.0/float64(rows), ctMean[(k-1)/matrixCols])
		fmt.Printf("Check the MeanVec\n")
		auxio.Quick_check_matrix(params, sk, ctMean[(k-1)/matrixCols], d, d) // dbg testing
		fmt.Printf("Mean %d/7 process Done in %s\n", (k-1)/matrixCols+1, time.Since(now))
		elapsed += time.Since(now)
		// We want to store the ctMean into disk, to reduce the memory occupation...
		// Store the mean:
		Meantracker.StoreUpdateOne(ctMean[(k-1)/matrixCols])
		ctMean[(k-1)/matrixCols] = nil                  // release the memory.
		galk4RowtotalSum = nil                          // free the Rotation key for RowtotalSum
		for i := 0; i < len(ctORSBT_SigmaDecomp); i++ { // free the decomposition space
			ctORSBT_SigmaDecomp[i] = nil
		}
		debug.FreeOSMemory()

	}

	// Load the Rotation keys for Transposition
	fmt.Printf("Retrieve the Rotation keys for Transposition\n")
	RtkTracker4Transpose = auxio.NewTracker4File(Rtk4Transpose_path)
	galk4Transpose = new(rlwe.RotationKeySet)
	_, err = RtkTracker4Transpose.ReadUpdateOne(galk4Transpose)
	if err != nil {
		panic(err)
	}
	// Generate evaluation key for transposition and other operation:
	elk := rlwe.EvaluationKey{Rlk: rlk, Rtks: galk4Transpose}
	evaluator = ckks.NewEvaluator(params, elk)

	// Using the symmetric property of Covariance to compute the rest of its part.
	now = time.Now()
	fmt.Printf("Using the symmetric property of Covariance to compute the rest of its part.\n")
	for i := 0; i < len(ctCov); i++ {
		for j := 0; j < len(ctCov); j++ {
			if i < j { // in the previous loop we've got ctCov[k][i] where i<k
				err = evaluator.Rescale(ctCov[j][i], params.DefaultScale(), ctCov[j][i])
				fmt.Printf("The ctCov[%d][%d], it will be used to compute ctCov[%d][%d]\n", j, i, i, j)
				auxio.Quick_check_matrix(params, sk, ctCov[j][i], d, d)
				if err != nil {
					panic(err)
				}
				ctCov[i][j] = evaluator.LinearTransform4ArithmeticSeqNew(ctCov[j][i], TransposeLT)[0]
				err = evaluator.Rescale(ctCov[i][j], params.DefaultScale(), ctCov[i][j])
				if err != nil {
					panic(err)
				}
				fmt.Printf("Syncronize the level and scale of ctCov[%d][%d] and ctCov[%d][%d]\n", i, j, j, i)
				evaluator.SetScale(ctCov[j][i], ctCov[i][j].Scale) // Syncronize the level and scale of ctCov[i][j] and ctCov[j][i]
				fmt.Printf("ctCov[%d][%d] now has scale %f, level %d, same as ctCov[%d][%d] scale %f, level %d\n", j, i, math.Log2(ctCov[j][i].Scale.Float64()), ctCov[j][i].Level(), i, j, math.Log2(ctCov[i][j].Scale.Float64()), ctCov[i][j].Level())
			}
		}
	}
	for i := 0; i < len(ctCov); i++ {
		err = evaluator.Rescale(ctCov[i][i], params.DefaultScale(), ctCov[i][i])
		if err != nil {
			panic(err)
		}
		if len(ctCov) > 1 {
			evaluator.SetScale(ctCov[i][i], ctCov[0][1].Scale) // Syncronize the level and scale of ctCov[i][i] and arbitary ctCov[i][j] where i < j.
			fmt.Printf("ctCov[%d][%d] is synchronized to scale %f, level %d\n", i, i, math.Log2(ctCov[i][i].Scale.Float64()), ctCov[i][i].Level())
		}
		for j := 0; j < len(ctCov); j++ {
			fmt.Printf("ctCov[%d][%d] stored into disk:\n", i, j)
			auxio.Quick_check_matrix(params, sk, ctCov[i][j], d, d)
			_, err = Covtracker.StoreUpdateOne(ctCov[i][j])
			if err != nil {
				panic(err)
			}
		}
	}
	fmt.Printf("Rest subblocks of Covariance computed and synchronized in %s\n", time.Since(now))
	elapsed += time.Since(now)

	fmt.Printf("Covariance Process Done in %s\n", elapsed)

	// Finish all the storing process:
	_, err = Covtracker.StoreFinish()
	if err != nil {
		panic(err)
	}
	_, err = Meantracker.StoreFinish()
	if err != nil {
		panic(err)
	}

	// Release the ctORSBT's memory, since all of the subblock ciphertexts are stored in disk.
	for i := 0; i < len(ctORSBT); i++ {
		ctORSBT[i] = nil
	}

}

func CovMatrix_check() {
	// initialize encryption scheme params.
	var err error
	matrixCols := 112
	// matrixRows := 128
	// rows := 256
	cols := 785

	// Size of the Square Matrix.
	d := 1 << 7

	// N1
	N1Trans := 16

	var params ckks.Parameters
	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var rlk *rlwe.RelinearizationKey
	// var kgen rlwe.KeyGenerator

	params, sk, pk, rlk, err = LoadKeys_check()
	if err != nil {
		panic(err)
	}

	// encoder
	encoder := ckks.NewEncoder(params)

	// create the Map and LinearTransform object of Transpose LinearTransformation
	fmt.Printf("create the Map and LinearTransform object of Transpose LinearTransformation\n")
	var TransposeDiagonalMap map[int][]float64
	TransposeDiagonalMap, err = mtrxmult.Gen_transpose_diagonalVectors(d)
	if err != nil {
		panic(err)
	}
	TransposeLT := ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, TransposeDiagonalMap, params.MaxLevel(), params.DefaultScale(), N1Trans, (d - 1), params.LogSlots())

	// create FileTrackers
	var Meantracker *auxio.Filetracker
	var Covtracker *auxio.Filetracker
	var RtkTracker4Transpose *auxio.Filetracker
	var RtkTracker4RowtotalSum *auxio.Filetracker
	var galk4Transpose *rlwe.RotationKeySet
	var galk4RowtotalSum *rlwe.RotationKeySet
	var galk4ColtotalSum *rlwe.RotationKeySet
	var RtkTracker4ColtotalSum *auxio.Filetracker

	// Construct ((cols-1)/matrixCols)^2 Ciphertexts representing a cols x cols Covariance Matrix
	ctCov := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctMean := make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctMMT := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols) // ciphertext for MeanVec · MeanVec^T
	ctEigVec := make([][]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctEigVal := make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	for i := 0; i < (cols-1)/matrixCols; i++ {
		ctCov[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
		ctMMT[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
		ctEigVec[i] = make([]*rlwe.Ciphertext, (cols-1)/matrixCols)

	}
	var ctTemp *rlwe.Ciphertext
	IdentityVec := make([]float64, 1<<params.LogSlots())
	for i := 0; i < len(IdentityVec); i++ {
		IdentityVec[i] = 0.5
	}
	// ptIdentity := encoder.EncodeNew(IdentityVec, params.MaxLevel(), params.DefaultScale(), params.LogSlots())

	// Compute Mean*Mean^T, first retrieving the Mean from disk.
	//_, err = auxio.ReadCtVec4file(Meantracker.Fp, Meantracker.Hierarchy[0], ctMean, 0)
	Meantracker = auxio.NewTracker4File("ctMean")
	RtkTracker4Transpose = auxio.NewTracker4File("Rtk4Transpose")
	RtkTracker4RowtotalSum = auxio.NewTracker4File("Rtk4RowtotalSum")
	//RtkTracker4ColtotalSum = auxio.NewTracker4File("Rtk4ColtotalSum")

	galk4Transpose = new(rlwe.RotationKeySet)
	RtkTracker4Transpose.ReadUpdateOne(galk4Transpose)
	galk4RowtotalSum = new(rlwe.RotationKeySet)
	galk4ColtotalSum = new(rlwe.RotationKeySet)
	RtkTracker4RowtotalSum.ReadUpdateOne(galk4RowtotalSum)
	//RtkTracker4ColtotalSum.ReadUpdateOne(galk4ColtotalSum)

	// Create RotationKeys for ColumnTotalSum

	kgen := ckks.NewKeyGenerator(params)
	RtkTracker4ColtotalSum = auxio.NewTracker("Rtk4ColtotalSum", -1)
	Rots4ColtotalSum := make([]int, 0)
	for i := 1; i < d; i = (i << 1) {
		Rots4ColtotalSum = append(Rots4ColtotalSum, i) // wairing to check
		Rots4ColtotalSum = append(Rots4ColtotalSum, -i)
	}
	galk4ColtotalSum = kgen.GenRotationKeysForRotations(Rots4ColtotalSum, false, sk)
	_, err = RtkTracker4ColtotalSum.StoreUpdateOne(galk4ColtotalSum)
	if err != nil {
		panic(err)
	}
	_, err = RtkTracker4ColtotalSum.StoreFinish()
	if err != nil {
		panic(err)
	}

	// evaluator
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk4Transpose})
	encryptor := ckks.NewEncryptor(params, pk)
	decryptor := ckks.NewDecryptor(params, sk)

	for i := 0; i < len(ctMean); i++ {
		ctMean[i] = new(rlwe.Ciphertext)
		Meantracker.ReadUpdateOne(ctMean[i])
		evaluator.Rescale(ctMean[i], params.DefaultScale(), ctMean[i])
		fmt.Printf("The %d-th ctMean with scale %f, level %d: \n", i, math.Log2(ctMean[i].Scale.Float64()), ctMean[i].Level())
	}

	for i := 0; i < len(ctMean); i++ {

		ctTemp = evaluator.LinearTransform4ArithmeticSeqNew(ctMean[i], TransposeLT)[0]
		fmt.Printf("The %d-th ctMean^T with scale %f, level %d: \n", i, math.Log2(ctTemp.Scale.Float64()), ctTemp.Level())
		auxio.Quick_check_matrix_full(params, sk, ctTemp, d, d)
		if err != nil {
			panic(err)
		}
		for j := 0; j < len(ctMean); j++ {
			// fmt.Printf("The i:%d, j:%d ctMMT: \n", i, j)
			ctMMT[i][j] = evaluator.MulRelinNew(ctTemp, ctMean[j])
			evaluator.Rescale(ctMMT[i][j], params.DefaultScale(), ctMMT[i][j])
			fmt.Printf("The ctMMT[%d][%d] with scale %f, level %d: \n", i, j, math.Log2(ctMMT[i][j].Scale.Float64()), ctMMT[i][j].Level())
		}
	}
	// free the Transpose Rotation keys:
	galk4Transpose = nil

	// Retrieve the X*X^T from disk, and do Cov = X*X^T-Mean*Mean^T
	Covtracker = auxio.NewTracker4File("ctCov")
	for i := 0; i < len(ctCov); i++ {
		//_, err = auxio.ReadCtVec4file(Covtracker.Fp, Covtracker.Hierarchy[0], ctCov[i], i*Covtracker.Hierarchy[0]*Covtracker.Hierarchy[1])
		for j := 0; j < len(ctCov); j++ {
			fmt.Printf("The i:%d, j:%d ctMM^T with scale %f, level %d: \n", i, j, math.Log2(ctMMT[i][j].Scale.Float64()), ctMMT[i][j].Level())
			auxio.Quick_check_matrix_full(params, sk, ctMMT[i][j], d, d)
			ctCov[i][j] = new(rlwe.Ciphertext)
			Covtracker.ReadUpdateOne(ctCov[i][j])
			fmt.Printf("The Original i:%d, j:%d ctX^TX has scale %f, level %d: \n", i, j, math.Log2(ctCov[i][j].Scale.Float64()), ctCov[i][j].Level())
			// Synchronizing scale Method 1:
			/*
				for k := 0; k <= ctCov[i][j].Level()-ctMMT[i][j].Level(); k++ {
					evaluator.Mul(ctCov[i][j], ptIdentity, ctCov[i][j])
					evaluator.MultByConst(ctCov[i][j], 2, ctCov[i][j])
					err = evaluator.Rescale(ctCov[i][j], params.DefaultScale(), ctCov[i][j])
					if err != nil {
						panic(err)
					}
					fmt.Printf("Drop ctX^TX to Scale %f, level %d\n", math.Log2(ctCov[i][j].Scale.Float64()), ctCov[i][j].Level())
				}
			*/
			// Synchronizing scale Method 2:
			/*
				err = evaluator.Rescale(ctCov[i][j], params.DefaultScale(), ctCov[i][j])
				if err != nil {
					panic(err)
				}
			*/
			evaluator.SetScale(ctMMT[i][j], ctCov[i][j].Scale)

			// evaluator.Rescale(ctCov[i][j], params.DefaultScale(), ctCov[i][j])
			fmt.Printf("The i:%d, j:%d ctX^TX with scale %f, level %d: \n", i, j, math.Log2(ctCov[i][j].Scale.Float64()), ctCov[i][j].Level())
			auxio.Quick_check_matrix_full(params, sk, ctCov[i][j], d, d)
			evaluator.Sub(ctCov[i][j], ctMMT[i][j], ctCov[i][j])
			fmt.Printf("The i:%d, j:%d ctCov: with scale %f, level %d: \n", i, j, math.Log2(ctCov[i][j].Scale.Float64()), ctCov[i][j].Level())
			//auxio.Quick_check_matrix_full(params, sk, ctCov[i][j], d, d)
			// do one recryption,this is for the following PowerMethod.
			// ctCov[i][j] = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctCov[i][j]), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
			auxio.Quick_check_matrix_full(params, sk, ctCov[i][j], d, d)
		}
	}

	// Compute the k most dominant eigenVectors and the corresponding eigenValues, We will skip the Boostrapping procedure, and use recryption instead.
	// create a Vi

	ptVi := auxio.Encode_single_float64(params, 1.0, params.MaxLevel(), params.DefaultScale())
	ctVi := make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctVip1 := make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	ctViT := make([]*rlwe.Ciphertext, (cols-1)/matrixCols)
	/*
		for i := 0; i < len(ctVi); i++ {
			ctVi[i] = encryptor.EncryptNew(ptVi)
		}
	*/

	// PowerMethod Iteration, iteratively compute: Vi <- Vi^T * Cov , notice that Cov is a semmetric matrix.
	var ctSum *rlwe.Ciphertext
	// var ctSign *rlwe.Ciphertext
	TargetEigVecNum := 2
	IterNum := 6
	NewtonIterNum := 15
	VecMod := 0 // We will have two Vector Mod: ColVector Mod (1) & RowVector Mod (0)

	for k := 0; k < TargetEigVecNum; k++ {
		for i := 0; i < len(ctVi); i++ {
			ctVi[i] = encryptor.EncryptNew(ptVi)
		}
		VecMod = 0
		for t := 0; t < IterNum; t++ {
			fmt.Printf("------------------ this is the %dth iteration ---------------\n", t)

			if t == IterNum-1 {
				NewtonIterNum = 20
			}

			if VecMod == 0 {
				/*
					if ctVi[0].Level() < 6 {
						for i := range ctVi {
							ctVi[i] = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctVi[i]), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
						}
					}
				*/
				for i := 0; i < len(ctCov); i++ {
					// We will do ctCov[j][i] * ctVi[j] , but will actually do ColAggregate(\sum(ctCov[i][j] * ctVi[j])) to complete this task, where ctCov[i][j]^T = ctCov[j][i]
					fmt.Printf("ctVi[%d]:\n", 0)
					auxio.Quick_check_infos(ctVi[0], "ctVi")
					ctSum = evaluator.MulNew(ctCov[i][0], ctVi[0]) // Consume one level
					fmt.Printf("ctCov[%d][%d]:\n", i, 0)
					auxio.Quick_check_infos(ctCov[i][0], "ctCov")
					// auxio.Quick_check_matrix(params, sk, ctCov[i][0], d, d)
					// auxio.Quick_check_matrix(params, sk, ctVi[0], d, d)
					fmt.Printf("MulAndAddResult:\n")
					auxio.Quick_check_infos(ctSum, "ctSum")
					// auxio.Quick_check_matrix(params, sk, ctSum, d, d)
					for j := 1; j < len(ctCov); j++ {
						fmt.Printf("ctCov[%d][%d]:\n", i, j)
						//auxio.Quick_check_matrix(params, sk, ctCov[i][j], d, d)
						fmt.Printf("ctVi[%d]:\n", j)
						// auxio.Quick_check_matrix(params, sk, ctVi[j], d, d)
						evaluator.MulAndAdd(ctCov[i][j], ctVi[j], ctSum)
						fmt.Printf("MulAndAddResult:\n")
						//auxio.Quick_check_matrix(params, sk, ctSum, d, d)
					}
					evaluator.Relinearize(ctSum, ctSum)
					//fmt.Printf("Before Aggregation:\n")
					//auxio.Quick_check_matrix(params, sk, ctSum, d, d)

					// Do one recryption
					if ctSum.Level() < params.MaxLevel() {
						ctSum = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctSum), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
					}

					err = evaluator.Rescale(ctSum, params.DefaultScale(), ctSum)
					if err != nil {
						panic(err)
					}
					ctSum, err = mtrxmult.Matrix_Aggregate_withRowOrder(params, sk, galk4ColtotalSum, ctSum, d, d, 1) // Consume one level
					if err != nil {
						panic(err)
					}
					fmt.Printf("The %dth ctSum Aggregated with scale %f level %d: \n", i, math.Log2(ctSum.Scale.Float64()), ctSum.Level())
					//auxio.Quick_check_matrix_full(params, sk, ctSum, d, d)

					ctVip1[i] = ctSum.CopyNew()
					fmt.Printf("The %dth Vip1 with scale %f level %d: \n", i, math.Log2(ctVip1[i].Scale.Float64()), ctVip1[i].Level())
					//auxio.Quick_check_matrix_full(params, sk, ctVip1[i], d, d)
					ctSum = nil

				}
				// now the ctVi has become the ColVector mod, we switch the sign:
				VecMod = 1

			} else if VecMod == 1 {
				/*
					if ctVi[0].Level() < 6 {
						fmt.Printf("Do one Recryption\n")
						for i := range ctVi {
							ctVi[i] = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctVi[i]), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
						}
					}
				*/
				for i := 0; i < len(ctCov); i++ {
					// we will do ctCov[i][j] * ctVi[j], but will actually do RowAggregate(\sum(ctCov[j][i] * ctVi[j])) to complete this task, where ctCov[i][j]^T = ctCov[j][i]
					ctSum = evaluator.MulNew(ctCov[0][i], ctVi[0]) // Consume one level
					fmt.Printf("ctCov[%d][%d]:\n", 0, i)
					//auxio.Quick_check_matrix(params, sk, ctCov[0][i], d, d)
					fmt.Printf("ctVi[%d]:\n", 0)
					//auxio.Quick_check_matrix(params, sk, ctVi[0], d, d)
					fmt.Printf("MulAndAddResult:\n")
					//auxio.Quick_check_matrix(params, sk, ctSum, d, d)
					for j := 1; j < len(ctCov); j++ {
						fmt.Printf("ctCov[%d][%d]:\n", j, i)
						//auxio.Quick_check_matrix(params, sk, ctCov[j][i], d, d)
						fmt.Printf("ctVi[%d]:\n", j)
						//auxio.Quick_check_matrix(params, sk, ctVi[j], d, d)
						evaluator.MulAndAdd(ctCov[j][i], ctVi[j], ctSum)
						fmt.Printf("MulAndAddResult:\n")
						//auxio.Quick_check_matrix(params, sk, ctSum, d, d)
					}
					evaluator.Relinearize(ctSum, ctSum)
					//fmt.Printf("Before Aggregation:\n")
					//auxio.Quick_check_matrix(params, sk, ctSum, d, d)

					// Do one recryption
					if ctSum.Level() < params.MaxLevel() {
						ctSum = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctSum), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
					}

					err = evaluator.Rescale(ctSum, params.DefaultScale(), ctSum)
					if err != nil {
						panic(err)
					}
					ctSum, err = mtrxmult.Matrix_Aggregate_withRowOrder(params, sk, galk4RowtotalSum, ctSum, d, d, 0)
					if err != nil {
						panic(err)
					}
					fmt.Printf("The %dth ctSum Aggregated with scale %f level %d: \n", i, math.Log2(ctSum.Scale.Float64()), ctSum.Level())
					auxio.Quick_check_matrix(params, sk, ctSum, d, d)

					ctVip1[i] = ctSum.CopyNew()
					fmt.Printf("The %dth Vip1 with scale %f level %d: \n", i, math.Log2(ctVip1[i].Scale.Float64()), ctVip1[i].Level())
					auxio.Quick_check_matrix(params, sk, ctVip1[i], d, d)
					ctSum = nil
					// ctSum = encryptor.EncryptZeroNew(params.MaxLevel())
				}
				// now the ctVi has become the RowVecotr mod, we switch the sign:
				VecMod = 0
			}

			// if this is not the last turn, we will do the normalisation. The last turn is only for computing eigen value and do not update eigen vector.
			if t < IterNum-1 {
				// Computing Inner product.
				ctSum = evaluator.MulNew(ctVip1[0], ctVip1[0])
				for i := 1; i < len(ctVip1); i++ { // Sum ctVip1
					evaluator.MulAndAdd(ctVip1[i], ctVip1[i], ctSum)
				}
				ctSum = evaluator.RelinearizeNew(ctSum)
				err = evaluator.Rescale(ctSum, params.DefaultScale(), ctSum) // Consume one level
				if err != nil {
					panic(err)
				}

				// Aggregate ctSum
				fmt.Printf("Before Inner Product:\n")
				auxio.Quick_check_matrix(params, sk, ctSum, d, d)
				// VecMod Now represents the ctVip1's form.
				if VecMod == 1 {
					ctSum, err = mtrxmult.Matrix_Aggregate_withRowOrder(params, sk, galk4RowtotalSum, ctSum, d, d, 0)
				} else if VecMod == 0 {
					ctSum, err = mtrxmult.Matrix_Aggregate_withRowOrder(params, sk, galk4ColtotalSum, ctSum, d, d, 1) // Consume one level
				}
				if err != nil {
					panic(err)
				}
				fmt.Printf("Inner Product with Scale %f, Level %d:\n", math.Log2(ctSum.Scale.Float64()), ctSum.Level())
				auxio.Quick_check_matrix(params, sk, ctSum, d, d)

				// Taylor Init Guess:

				ctGuess, err := nonpolyfunc.TaylorInitNew(params, rlk, ctSum, nonpolyfunc.Inv_sqrt_taylor2_0to2pow10[:]) // Consume one level *2
				if err != nil {
					panic(err)
				}
				fmt.Printf("Taylor Guess with Scale %f, Level %d:\n", math.Log2(ctGuess.Scale.Float64()), ctGuess.Level())
				auxio.Quick_check_matrix(params, sk, ctGuess, d, d)

				// Do one Recryption:
				ctGuess = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctGuess), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
				ctSum = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctSum), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))

				// InvSqrt by Newton:
				ctInvsqrt, err := nonpolyfunc.InvSqrtByNewton_dbg(params, rlk, pk, sk, ctSum, ctGuess, NewtonIterNum) // Consume one level * 3 * NewtonIterNum
				if err != nil {
					panic(err)
				}
				fmt.Printf("InvSqrt with Scale %f, Level %d:\n", math.Log2(ctInvsqrt.Scale.Float64()), ctInvsqrt.Level())
				auxio.Quick_check_matrix(params, sk, ctInvsqrt, d, d)

				// Do one Recryption:
				ctInvsqrt = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctInvsqrt), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
				fmt.Printf("InvSqrt after recryption with Scale %f, Level %d:\n", math.Log2(ctInvsqrt.Scale.Float64()), ctInvsqrt.Level())

				// Normalisation
				for i := 0; i < len(ctVip1); i++ {
					evaluator.MulRelin(ctVip1[i], ctInvsqrt, ctVip1[i])
					err = evaluator.Rescale(ctVip1[i], params.DefaultScale(), ctVip1[i]) // Consume one level
					if err != nil {
						panic(err)
					}
				}
				fmt.Printf("ctVip1 after normalisation has Scale %f, Level %d, this will be used to update ctVi\n", math.Log2(ctVip1[0].Scale.Float64()), ctVip1[0].Level())

				// Update ctVi with ctVip1
				for i := 0; i < len(ctVip1); i++ {
					ctVi[i] = ctVip1[i].CopyNew()
					fmt.Printf("The %dth iteration's Vip1 %d result:\n", t, i)
					// auxio.Quick_check_matrix_full(params, sk, ctVi[i], d, d)
				}

			} else { // if this is the last turn, compute the eigen value corresponding to the eigen vector.
				// Recall that at this point, ctVip1 should contain Cov @ ctVi, the eigen value we want can be computed by :
				// eig_val = <(Cov @ ctVi), ctVi> / <ctVi, ctVi> = <ctVip1, ctVi> / <ctVi,ctVi>, if we have confidence in the
				// normalisation procedure in previous turn, than this equation can be reduced to eig_val = <ctVip1,ctVi> / 1
				// After this, we will have to compute the Shifted Covariance Matrix: ShiftedCov = Cov - eigval * eigvec^T @ eigvec

				fmt.Printf("This is the last iteration of %dth EigenVector computation\n", k)
				// So we will Compute the Inner product <ctVip1,ctVi>
				// VecMod Now represents the ctVip1's form.
				if ((VecMod + 1) & (2 - 1)) == 1 { // if ctVi is now in column mod, then we create a replication of its row mod version
					// this will cost one level.
					for i := 0; i < len(ctVip1); i++ {
						ctViT[i], err = mtrxmult.ReplicateVec_Switch_Axis(params, sk, galk4RowtotalSum, ctVi[i], d, (VecMod+1)&(2-1)) // Consume one level
						fmt.Printf("The ctViT with Scale %f, level %d:\n", math.Log2(ctViT[i].Scale.Float64()), ctViT[i].Level())
						auxio.Quick_check_matrix(params, sk, ctViT[i], d, d)
						fmt.Printf("The ctVi  with Scale %f, level %d:\n", math.Log2(ctVi[i].Scale.Float64()), ctVi[i].Level())
						auxio.Quick_check_matrix(params, sk, ctVi[i], d, d)
						if err != nil {
							panic(err)
						}
					}
				} else if ((VecMod + 1) & (2 - 1)) == 0 { // if ctVi is now in row mod, then we create a replication of its column mod version
					// this will cost one level.
					for i := 0; i < len(ctVip1); i++ {
						ctViT[i], err = mtrxmult.ReplicateVec_Switch_Axis(params, sk, galk4ColtotalSum, ctVi[i], d, (VecMod+1)&(2-1)) // Consume one level *2
						if err != nil {
							panic(err)
						}
					}
				}

				ctSum = evaluator.MulNew(ctVip1[0], ctViT[0])
				for i := 1; i < len(ctVip1); i++ { // Sum ctVip1
					evaluator.MulAndAdd(ctVip1[i], ctViT[i], ctSum)
				}
				ctSum = evaluator.RelinearizeNew(ctSum)
				err = evaluator.Rescale(ctSum, params.DefaultScale(), ctSum) // Consume one level
				if err != nil {
					panic(err)
				}
				ctSum = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctSum), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
				fmt.Printf("Inner product before aggregation has scale %f, level %d after recryption\n", math.Log2(ctSum.Scale.Float64()), ctSum.Level())

				// Aggregate ctSum then we get the eig_val.
				// fmt.Printf("Before Inner Product:\n")
				// auxio.Quick_check_matrix(params, sk, ctSum, d, d)
				// VecMod Now represents the ctVip1's form.
				if VecMod == 1 { // since ctVip1 is in column mod, then ctSum is the multiple of ctVip1 and ctViT that is in Column mod
					ctSum, err = mtrxmult.Matrix_Aggregate_withRowOrder(params, sk, galk4RowtotalSum, ctSum, d, d, 0)
				} else if VecMod == 0 { // since ctVip1 is in Row mod, then ctSum is the multiple of ctVip1 and ctViT that is in row mod
					ctSum, err = mtrxmult.Matrix_Aggregate_withRowOrder(params, sk, galk4ColtotalSum, ctSum, d, d, 1) // Consume one level
				}
				if err != nil {
					panic(err)
				}
				fmt.Printf("Aggregate ctSum then we get the eig_val with Scale %f, level %d:\n", math.Log2(ctSum.Scale.Float64()), ctSum.Level())
				// DO one Recryption.
				auxio.Quick_check_matrix(params, sk, ctSum, d, d)

				// Store the eig_val:
				ctEigVal[k] = ctSum.CopyNew()
				// Store the eig_Vec:
				for i := 0; i < len(ctVi); i++ {
					ctEigVec[k][i] = ctVi[i].CopyNew()
				}

				// Update Cov as the Shifted Cov, We will use container ctMMT to temporarily store the eigvec^T @ eigvec
				for i := 0; i < len(ctEigVec[k]); i++ {
					fmt.Printf("Check the ctEigVec[%d][%d] to see wether it is really of form %d:\n", k, i, (VecMod+1)&(2-1))
					// auxio.Quick_check_matrix_full(params, sk, ctEigVec[k][i], d, d)
					ctTemp = ctViT[i]
					fmt.Printf("The ctEigVec^T[%d][%d] is with scale %f, level %d: \n", k, i, math.Log2(ctTemp.Scale.Float64()), ctTemp.Level())
					auxio.Quick_check_matrix(params, sk, ctTemp, d, d)
					if err != nil {
						panic(err)
					}
					for j := 0; j < len(ctEigVec[k]); j++ {
						// fmt.Printf("The i:%d, j:%d ctMMT: \n", i, j)
						if ((VecMod + 1) & (2 - 1)) == 1 {
							fmt.Printf("Internal turn: The ctEigVec^T[%d][%d], it will be multiplied with ctEigVec[%d][%d] to become ctMMT[%d][%d]\n:", k, i, k, j, j, i)
							ctMMT[j][i] = evaluator.MulRelinNew(ctTemp, ctEigVec[k][j])
							evaluator.Rescale(ctMMT[j][i], params.DefaultScale(), ctMMT[j][i]) // Consume one level
							ctMMT[j][i] = evaluator.MulRelinNew(ctMMT[j][i], ctEigVal[k])
							evaluator.Rescale(ctMMT[j][i], params.DefaultScale(), ctMMT[j][i]) // Consume one level
							fmt.Printf("The eigvec^T @ eigvec * eigval[%d][%d] with scale %f, level %d: \n", i, j, math.Log2(ctMMT[j][i].Scale.Float64()), ctMMT[j][i].Level())
							if ctMMT[j][i].Level() >= ctCov[j][i].Level() {
								ctMMT[j][i].SetScale(ctCov[j][i].Scale)
							} else {
								ctCov[j][i].SetScale(ctMMT[j][i].Scale)
							}
							auxio.Quick_check_matrix(params, sk, ctMMT[j][i], d, d)

							fmt.Printf("The Original ctCov[%d][%d] with scale %f, level %d: \n", i, j, math.Log2(ctCov[j][i].Scale.Float64()), ctCov[j][i].Level())
							auxio.Quick_check_matrix(params, sk, ctCov[j][i], d, d)
							ctCov[j][i] = evaluator.SubNew(ctCov[j][i], ctMMT[j][i])
							fmt.Printf("The Shifted ctCov[%d][%d] with scale %f, level %d: \n", i, j, math.Log2(ctCov[j][i].Scale.Float64()), ctCov[j][i].Level())
							auxio.Quick_check_matrix(params, sk, ctCov[j][i], d, d)
						} else if ((VecMod + 1) & (2 - 1)) == 0 { // idx should be flipped here .
							ctMMT[i][j] = evaluator.MulRelinNew(ctTemp, ctEigVec[k][j])
							evaluator.Rescale(ctMMT[i][j], params.DefaultScale(), ctMMT[i][j]) // Consume one level
							ctMMT[i][j] = evaluator.MulRelinNew(ctMMT[i][j], ctEigVal[k])
							evaluator.Rescale(ctMMT[i][j], params.DefaultScale(), ctMMT[i][j]) // Consume one level
							fmt.Printf("The eigvec^T @ eigvec * eigval[%d][%d] with scale %f, level %d: \n", i, j, math.Log2(ctMMT[i][j].Scale.Float64()), ctMMT[i][j].Level())
							if ctMMT[i][j].Level() >= ctCov[i][j].Level() {
								ctMMT[i][j].SetScale(ctCov[i][j].Scale)
							} else {
								ctCov[i][j].SetScale(ctMMT[i][j].Scale)
							}
							auxio.Quick_check_matrix(params, sk, ctMMT[i][j], d, d)

							fmt.Printf("The Original ctCov[%d][%d] with scale %f, level %d: \n", i, j, math.Log2(ctCov[j][i].Scale.Float64()), ctCov[i][j].Level())
							auxio.Quick_check_matrix(params, sk, ctCov[i][j], d, d)
							ctCov[i][j] = evaluator.SubNew(ctCov[i][j], ctMMT[i][j])
							fmt.Printf("The Shifted ctCov[%d][%d] with scale %f, level %d: \n", i, j, math.Log2(ctCov[i][j].Scale.Float64()), ctCov[i][j].Level())
							auxio.Quick_check_matrix(params, sk, ctCov[i][j], d, d)
						}

					}
				}
			}

		}
	}

}

func LoadKeys_check() (params ckks.Parameters, sk *rlwe.SecretKey, pk *rlwe.PublicKey, rlk *rlwe.RelinearizationKey, err error) {
	// var err error
	// var params ckks.Parameters
	// var sk *rlwe.SecretKey
	// var pk *rlwe.PublicKey
	// var rlk *rlwe.RelinearizationKey
	// var galk4Transpose *rlwe.RotationKeySet
	// var galk4MtrxMult *rlwe.RotationKeySet
	// var galk4RowtotalSum *rlwe.RotationKeySet

	// Allocate Memory for keys
	sk = new(rlwe.SecretKey)
	pk = new(rlwe.PublicKey)
	rlk = new(rlwe.RelinearizationKey)

	// create Keys' File Tracker.
	var keysTracker = auxio.NewTracker4File("Keys")
	// var RtkTracker4Transpose = auxio.NewTracker4File("")
	// var RtkTracker4MtrxMult = auxio.NewTracker4File("")
	// var RtkTracker4RowtotalSum = auxio.NewTracker4File("")
	var n int // record IO bytes.

	// Load sk,pk and rlk.
	n, err = keysTracker.ReadUpdateOne(&params)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read params %d bytes\n", n)
	}
	n, err = keysTracker.ReadUpdateOne(sk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read sk %d bytes\n", n)
	}
	n, err = keysTracker.ReadUpdateOne(pk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read pk %d bytes\n", n)
	}
	n, err = keysTracker.ReadUpdateOne(rlk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read rlk %d bytes\n", n)
	}
	//
	return

}

func Transepose_check() {
	var err error
	var params ckks.Parameters
	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var rlk *rlwe.RelinearizationKey
	var galk4Transpose *rlwe.RotationKeySet
	// var now time.Time
	// create a 4x4 matrix matrixA
	d := 128
	// Sigma_BSGSRatio := 2
	// Tau_BSGSRatio := 2
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

	params, sk, pk, rlk, err = LoadKeys_check()
	if err != nil {
		panic(err)
	}
	encoder := ckks.NewEncoder(params)
	encryptor := ckks.NewEncryptor(params, pk)

	// Load the Rotation keys for Transposition
	fmt.Printf("Retrieve the Rotation keys for Transposition\n")
	RtkTracker4Transpose := auxio.NewTracker4File("Rtk4Transpose")
	galk4Transpose = new(rlwe.RotationKeySet)
	_, err = RtkTracker4Transpose.ReadUpdateOne(galk4Transpose)
	if err != nil {
		panic(err)
	}
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk4Transpose})

	// create the Map and LinearTransform object of Transpose LinearTransformation
	fmt.Printf("create the Map and LinearTransform object of Transpose LinearTransformation\n")
	var TransposeDiagonalMap map[int][]float64
	TransposeDiagonalMap, err = mtrxmult.Gen_transpose_diagonalVectors(d)
	if err != nil {
		panic(err)
	}
	TransposeLT := ckks.GenLinearTransformBSGS4ArithmeticSeq(encoder, TransposeDiagonalMap, params.MaxLevel(), params.DefaultScale(), 16, (d - 1), params.LogSlots())
	ptA := encoder.EncodeNew(a, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	ctA := encryptor.EncryptNew(ptA)
	auxio.Quick_check_matrix_full(params, sk, ctA, d, d)
	ctAT := evaluator.LinearTransform4ArithmeticSeqNew(ctA, TransposeLT)[0]
	auxio.Quick_check_matrix_full(params, sk, ctAT, d, d)

}
