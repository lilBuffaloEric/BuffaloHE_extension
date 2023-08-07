package check

/*
 * Check for funcs in auxiliary_io
 */

import (
	"fmt"

	auxio "project1-fhe_extension_v1.0/auxiliary_io"
	mtrxmult "project1-fhe_extension_v1.0/matrix_mult"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

var (
	// PN15QP320S14 is a default parameter set for logN=15, logQP=320 and LogSlots=14
	PN15QP320S14 = ckks.ParametersLiteral{
		LogN: 15,
		Q: []uint64{1152921504598720513,
			1099490000897, 1099498258433,
			1099499175937, 1099499569153,
			1099500617729}, // 60 + 40*5
		P:            []uint64{1152921504606584833}, // 60
		LogSlots:     14,
		DefaultScale: 1 << 40,
	}
	// PN15QP320S14 is a default parameter set for logN=15, logQP=440 and LogSlots=14
	PN15QP440S14 = ckks.ParametersLiteral{
		LogN: 15,
		Q: []uint64{
			0x4000000120001, 0x10000140001, 0xffffe80001, // 50 + 6 x 40
			0x10000290001, 0xffffc40001, 0x100003e0001,
			0x10000470001,
		},
		P:            []uint64{0x40000001b0001, 0x3ffffffdf0001, 0x4000000270001}, // 50, 50, 50
		LogSlots:     14,
		DefaultScale: 1 << 40,
	}
	// PN15QP880 is a default parameter set for logN=15 and logQP=720 and LogSlots=14
	PN15QP720S14 = ckks.ParametersLiteral{
		LogN: 15,
		Q: []uint64{0x4000000120001, 0x10000140001, 0xffffe80001, // 50 + 13 x 40
			0x10000290001, 0xffffc40001, 0x100003e0001,
			0x10000470001, 0x100004b0001, 0xffffb20001,
			0x10000500001, 0x10000650001, 0xffff940001,
			0xffff8a0001, 0xffff820001},
		P:            []uint64{0x40000001b0001, 0x3ffffffdf0001, 0x4000000270001}, // 50, 50, 50
		LogSlots:     14,
		DefaultScale: 1 << 40,
	}
)

func Csvio_check() {

	// check CSV Readin
	csvfilepath := "D:\\cardistry\\discovery@GY\\软件学院 2th\\机器学习基础\\机器学习实验三\\mnist\\mnist_train.csv"
	cols := 785
	rows := 60000
	matrixCols := 112
	matrixRows := 128
	for k := 1; k < cols; k += matrixCols {
		for j := 1; j < rows; j += matrixRows {
			if j+matrixRows > rows {
				auxio.GetOneCSVSubblock(csvfilepath, j, k, rows-j+1, matrixCols)
			} else {
				auxio.GetOneCSVSubblock(csvfilepath, j, k, matrixRows, matrixCols)
			}
		}

	}

}

func Matrixvisualize_check() {
	d := 8
	// New version Sigma Decomposition check:
	Z1s, Z2s, err := mtrxmult.GenSigmaDiagnalDecomposeMatrices(d, d/4)
	if err != nil {
		panic(err)
	}
	Z1 := mtrxmult.DiagonalVectors2Matrix(Z1s[len(Z1s)-1], d*d)
	for i := len(Z1s) - 2; i >= 0; i-- {
		Z1 = mtrxmult.SlowPlainMatrixMult(mtrxmult.DiagonalVectors2Matrix(Z1s[i], d*d), Z1, d*d)
	}
	Z2 := mtrxmult.DiagonalVectors2Matrix(Z2s[len(Z2s)-1], d*d)
	for i := len(Z2s) - 2; i >= 0; i-- {
		Z2 = mtrxmult.SlowPlainMatrixMult(mtrxmult.DiagonalVectors2Matrix(Z2s[i], d*d), Z2, d*d)
	}
	Z := mtrxmult.SlowPlainMatrixAdd(Z1, Z2, d*d)
	auxio.Print_matrix_f64_full_2d(Z, d*d, d*d, 0)
	Ori_DS, err := mtrxmult.Gen_sigma_diagonalVecotrs(d)
	if err != nil {
		panic(err)
	}
	Ori_Z := mtrxmult.DiagonalVectors2Matrix(Ori_DS, d*d)
	print(mtrxmult.MatrixCompare(Z, Ori_Z))

	// New version Tau Decompostion check:
	Ts, err := mtrxmult.GenTauDiagonalDecomposeMatrices(d, d/4*d)
	if err != nil {
		panic(err)
	}
	T := mtrxmult.DiagonalVectors2Matrix(Ts[len(Ts)-1], d*d)
	for i := len(Ts) - 2; i >= 0; i-- {
		T = mtrxmult.SlowPlainMatrixMult(mtrxmult.DiagonalVectors2Matrix(Ts[i], d*d), T, d*d)
	}
	Ori_DT, err := mtrxmult.Gen_tao_diagonalVectors(d)
	if err != nil {
		panic(err)
	}
	Ori_T := mtrxmult.DiagonalVectors2Matrix(Ori_DT, d*d)
	print(mtrxmult.MatrixCompare(T, Ori_T))

	U, err := mtrxmult.Gen_tao_diagonalVectors(d)
	if err != nil {
		panic(err)
	}
	A := mtrxmult.DiagonalVectors2Matrix(U, d*d)
	//auxio.Print_matrix_f64_full_2d(A, d*d, d*d)
	/*
		U2, _ := mtrxmult.Gen_tao_diagonalVectors(d)
		B := mtrxmult.DiagonalVectors2Matrix(U2, d*d)
		auxio.Print_matrix_f64_full_2d(B, d*d, d*d)
		UT, _ := mtrxmult.Gen_transpose_diagonalVectors(d)
		if err != nil {
			panic(err)
		}
		AT := mtrxmult.DiagonalVectors2Matrix(UT, d*d)
		auxio.Print_matrix_f64_full_2d(AT, d*d, d*d)
	*/
	/*
		for k := 1; k < d; k++ {
			UC, _ := mtrxmult.Gen_colShift_diagonalVectors(d, k)
			C := mtrxmult.DiagonalVectors2Matrix(UC, d*d)
			CA := mtrxmult.SlowPlainMatrixMult(C, A, d*d)
			auxio.Print_matrix_f64_full_2d(CA, d*d, d*d)
			UR, _ := mtrxmult.Gen_rowShift_diagonalVectors(d, k)
			R := mtrxmult.DiagonalVectors2Matrix(UR, d*d)
			RB := mtrxmult.SlowPlainMatrixMult(R, B, d*d)
			auxio.Print_matrix_f64_full_2d(RB, d*d, d*d)
		}

		TargetDiagVec := make([]int, int(math.Ceil(float64((len(U)-1)/2)/2))*2+1)
		for i, j := 0, 0; i < len(TargetDiagVec); i++ {
			if i == 0 {
				TargetDiagVec[i] = 0
				j++
			} else {
				TargetDiagVec[i] = j
				j++
				i++
				TargetDiagVec[i] = -j + d*d

			}
		} // (16-1)
	*/

	Uconverged, _ := mtrxmult.Converge2DiagonalDecompose_Tao(A)
	D0 := mtrxmult.DiagonalVectors2Matrix(Uconverged[0], d*d)
	D1 := mtrxmult.DiagonalVectors2Matrix(Uconverged[1], d*d)
	//auxio.Print_matrix_f64_full_2d(D0, d*d, d*d)
	//auxio.Print_matrix_f64_full_2d(D1, d*d, d*d)
	CoveredA := mtrxmult.SlowPlainMatrixMult(D1, D0, d*d)
	//auxio.Print_matrix_f64_full_2d(CoveredA, d*d, d*d)

	print(mtrxmult.MatrixCompare(CoveredA, A)) // not necessary.
}

/*
func CiphertextIO_check() {
	// initialize encryption scheme.
	var err error
	d := 4
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
	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)
	// encoder
	// encoder := ckks.NewEncoder(params)
	// Evaluator
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})

	// encryption
	ctArr := make([]*rlwe.Ciphertext, 3)
	for i := 0; i < 3; i++ {
		ctArr[i] = evaluator.AddConstNew(encryptor.EncryptZeroNew(params.MaxLevel()), i)
		auxio.Quick_check_vector(params, sk, ctArr[i])
	}
	// Create Tracker for ctArr:
	var Tracker = auxio.Filetracker{Fp: "CiphertextIO_test", Off: 0, Hierarchy: make([]int, 1)}
	Tracker.Hierarchy[0] = ctArr[0].MarshalBinarySize()

	// Write ctArr into file.
	Tracker.Off, err = auxio.WriteCtVec2file(Tracker.Fp, ctArr, 0)
	if err != nil {
		panic(err)
	}
	fmt.Printf("Write in bytes %d\n", Tracker.Off)

	// Zero ctArr to see difference.
	fmt.Printf("Now zero the ciphertext vector to see the difference.\n")
	for i := 0; i < 3; i++ {
		ctArr[i] = encryptor.EncryptZeroNew(params.MaxLevel())
		auxio.Quick_check_vector(params, sk, ctArr[i])
	}

	// Read ctArr[1:3] into file, so we set the offset to the end of first ciphertext.
	Tracker.Off, err = auxio.ReadCtVec4file(Tracker.Fp, Tracker.Hierarchy[0], ctArr[1:3], Tracker.Hierarchy[0])
	if err != nil {
		panic(err)
	}
	fmt.Printf("Read in bytes %d\n", Tracker.Off)

	// Check the correctness
	for i := 0; i < 3; i++ {
		auxio.Quick_check_vector(params, sk, ctArr[i])
	}
}
*/

func DiskIO_check() {

	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var rlk *rlwe.RelinearizationKey
	var galk *rlwe.RotationKeySet
	var keysTracker = auxio.NewTracker("Keys", -1)
	// var CiphertextsTracker = auxio.NewTracker("Cihpertexts",-1)
	params, err := ckks.NewParametersFromLiteral(ckks.PN14QP438)
	if err != nil {
		panic(err)
	}
	kgen := ckks.NewKeyGenerator(params)
	sk, pk = kgen.GenKeyPair()
	rlk = kgen.GenRelinearizationKey(sk, 1)
	Rotations := make([]uint64, 2)
	Rotations[0] = 1
	Rotations[1] = 2
	galk = kgen.GenRotationKeys(Rotations, sk)

	// store these keys into file "Keys"
	n := 0
	n, err = keysTracker.StoreUpdateOne(sk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Store sk %d bytes", n)
	}
	n, err = keysTracker.StoreUpdateOne(pk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Store pk %d bytes", n)
	}
	n, err = keysTracker.StoreUpdateOne(rlk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Store rlk %d bytes", n)
	}
	n, err = keysTracker.StoreUpdateOne(galk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Store galk %d bytes", n)
	}
	n, err = keysTracker.StoreFinish()
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Store filetracker %d bytes", n)
	}
	// Get these keys from file "Keys"
	var Retrieved_sk *rlwe.SecretKey
	var Retrieved_pk *rlwe.PublicKey
	var Retrieved_rlk *rlwe.RelinearizationKey
	var Rectrieved_galk *rlwe.RotationKeySet
	Retrieved_sk = new(rlwe.SecretKey)
	Retrieved_pk = new(rlwe.PublicKey)
	Retrieved_rlk = new(rlwe.RelinearizationKey)
	Rectrieved_galk = new(rlwe.RotationKeySet)

	keysTracker = nil
	keysTracker = auxio.NewTracker4File("Keys")
	n, err = keysTracker.ReadUpdateOne(Retrieved_sk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read sk %d bytes", n)
	}
	n, err = keysTracker.ReadUpdateOne(Retrieved_pk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read pk %d bytes", n)
	}
	n, err = keysTracker.ReadUpdateOne(Retrieved_rlk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read rlk %d bytes", n)
	}
	n, err = keysTracker.ReadUpdateOne(Rectrieved_galk)
	if err != nil {
		panic(err)
	} else {
		fmt.Printf("Read galk %d bytes", n)
	}
	// test the original keys and retrieved keys' functionality
	encoder := ckks.NewEncoder(params)
	encryptor := ckks.NewEncryptor(params, pk)
	vec := make([]float64, params.MaxSlots())
	for i := 0; i < len(vec); i++ {
		vec[i] = float64(i)
	}
	pt := encoder.EncodeNew(vec, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	ct := encryptor.EncryptNew(pt)
	auxio.Quick_check_vector(params, Retrieved_sk, ct)
	// test Store many and Read many
	/*
		evaulator := ckks.NewEvaluator(params,rlwe.EvaluationKey{Rlk: Retrieved_rlk})
		ct_List := make([]*rlwe.Ciphertext,3)
		for i:=0;i<len(ct_List);i++{
			ct_List[i] = evaulator.AddConstNew(ct,1)
		}
		CiphertextsTracker.StoreUpdateMany(ct_List)
	*/

}

func Correctness_check() {
	// initialize encryption scheme params.
	var params ckks.Parameters
	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var rlk *rlwe.RelinearizationKey
	var err error
	params, sk, pk, rlk, err = LoadKeys_check()
	if err != nil {
		panic(err)
	}
	var d = 1 << (params.LogSlots() / 2)
	var Covtracker = auxio.NewTracker4File("ctCov")
	var ctCovSubblock = new(rlwe.Ciphertext)
	_, err = Covtracker.ReadUpdateOne(ctCovSubblock)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_matrix_full(params, sk, ctCovSubblock, d, d)
	print(pk)  // not necessary
	print(rlk) // not necessary
}
