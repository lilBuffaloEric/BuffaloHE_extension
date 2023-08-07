package main

import "project1-fhe_extension_v1.0/check"

func main() {
	// ------------------------------ Internal Components Check (Not all available) ------------------------------
	// check.CovMatrix_check()
	// check.MatrixMultX_check(true)
	// check.LinTransX_check()
	// check.LinTransX2_check()
	// check.LinTransX3_check()
	// check.Transepose_check()
	// check.OnlyMatrixMult_check()
	// check.Invsqrt_check() // checked
	// check.SquareMatrix_product_check() // checked
	// check.SquareMatrix_encode_check() // checked
	// check.SquareMatrix_permute_check() // checked
	// check.Piecewise_check() // checked
	// check.Softmax_check() // checked
	// check.Csvio_check()
	// check.PCA_check(false) // 202356: update N1 as 32 for all Lintrans.
	// check.PCA_largeSpace_check(true)
	// check.CiphertextIO_check()
	// check.Matrixvisualize_check()
	// check.DiskIO_check()
	// check.Correctness_check()
	// check.ErrorCtrlX_check(true)
	// check.StoreAndLoadKeys_btpVersion_check(true)

	// ------------------------------ Available for testing ------------------------------
	// check.LinTransX3_check() // Routine Testing the Diagonal Convergence Decomposition of the Lintrans in homomorphic matrix multiplication.
	// check.PCA_btpVersion_check(false, false) // Routine Computing Cov' = 1/N*X^T*X and U = u*u^T
	check.PCA_largeSpace_Ver2_check(false) // Routine Computing Cov' = 1/N*X^T*X and U = u*u^T (theoratically faster, but need more space to cache intermediate homomorphic matrix-mult results in memory)
	// check.CovMatrix_btpVersion_check(false, false) // Routine Computing Cov = Cov' - U and the approximate eigenvectors of Cov. Use bootstrpping if reencryptionMode is set to "Ture", but correctness cannot be guaranteed.
	// check.PandaPCA_check(true) // Routine simulating the Panda's scheme (2021). This can only be reguarded as a speed testing, and no correctness can be guaranteed.
}

/*
import (
	"errors"
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

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

func Sigma_linearTransform(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, d int) (ctOut *rlwe.Ciphertext, err error) {
	if d <= 0 || float64(d*d) != math.Pow(2, float64(params.LogSlots())) {
		return nil, errors.New("d<0 or d^2 != 2^logSlots")
	}
	var U_sigma map[int][]float64
	U_sigma, err = Gen_sigma_diagonalVecotrs(d)
	if err != nil {
		return nil, err
	}
	Scale := ctIn.Scale
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	sigmaLT := ckks.GenLinearTransform(encoder, U_sigma, ctIn.Level(), Scale, params.LogSlots())
	ctOut_list := evaluator.LinearTransformNew(ctIn, sigmaLT)
	ctOut = ctOut_list[0]
	err = evaluator.Rescale(ctOut, ctIn.Scale, ctOut)
	if err != nil {
		return nil, err
	}
	fmt.Printf("%d levels consumed for LinearTransform\n", ctIn.Level()-ctOut.Level())
	return
}

func main() {
	var err error
	d := 4
	MatrixA := make([]float64, d*d)
	for i := 0; i < d; i++ {
		for j := 0; j < d; j++ {
			MatrixA[i*d+j] = float64(i*d + j)
		}
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
	// rotation steps {-d+1 ~ d-1}
	steps := make([]int, d+d-1+d-1)
	j := 0
	for i := 0; i < d; i++ {
		steps[j] = i
		j++
		if i != 0 {
			steps[j] = -i // FIXME: -i is original but will cause panic, so we will have to set -i mod d^2 here.
			j++
		}
	}
	galk := kgen.GenRotationKeysForRotations(steps[:], false, sk)

	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)

	// encoder
	encoder := ckks.NewEncoder(params)
	// encryption
	ptA := encoder.EncodeNew(MatrixA, params.MaxLevel(), params.DefaultScale(), params.LogSlots())
	ctA := encryptor.EncryptNew(ptA)
	// Check Sigma permutation
	var ctSigmaA *rlwe.Ciphertext
	ctSigmaA, err = Sigma_linearTransform(params, rlk, galk, ctA, d)
	if err != nil {
		panic(err)
	}
	fmt.Print(ctSigmaA) // not necessary...
}
*/
