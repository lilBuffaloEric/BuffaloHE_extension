package check

/*
 * Check for funcs in auxiliary_io
 */

import (
	"math"
	"os"
	auxio "project1-fhe_extension/auxiliary_io"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
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
	// rlk := kgen.GenRelinearizationKey(sk, 1)
	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)
	// encoder
	// encoder := ckks.NewEncoder(params)
	// Evaluator
	// evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})

	// encryption
	ptA := auxio.Encode_single_int(params, 1, params.MaxLevel())
	ctA := encryptor.EncryptNew(ptA)
	ctArr := make([]*rlwe.Ciphertext, 3)
	for i := 0; i < 3; i++ {
		ctArr[i] = ctA
	}
	file, _ := os.Create("CiphertextIO_test.txt")

	var MptA []byte
	var MctA []byte
	var MctArr = make([][]byte, 3)
	MptA, err = ptA.MarshalBinary()
	MctA, err = ctA.MarshalBinary()

	file.Write(MptA)
	file.Write(MctA)

	for i := 0; i < 3; i++ {
		MctArr[i], err = ctArr[i].MarshalBinary()
		file.Write(MctArr[i])
	}
	file.Close()

	file, err = os.Open("CiphertextIO_test.txt")
	CMptA := make([]byte, len(MptA))
	CMctA := make([]byte, len(MctA))
	file.ReadAt(CMctA, int64(len(MptA)))
	file.ReadAt(CMptA, 0)
	ctA.UnmarshalBinary(CMctA)
	ptA.UnmarshalBinary(CMptA)
	auxio.Quick_check_vector(params, sk, ctA)
	auxio.Quick_decode_vector(params, ptA)

}
