package check

import (
	auxio "project1-fhe_extension_v1.0/auxiliary_io"
	piecewisefunc "project1-fhe_extension_v1.0/piecewise_func"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func Piecewise_check() {
	params, err := ckks.NewParametersFromLiteral(ckks.PN15QP880)
	if err != nil {
		panic(err)
	}

	// encoder := ckks.NewEncoder(params)

	// Keys
	kgen := ckks.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPair()

	// Relinearization key
	rlk := kgen.GenRelinearizationKey(sk, 1)

	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)

	// Decryptor
	// decryptor := ckks.NewDecryptor(params, sk)

	// Evaluator
	// evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})

	ptA := auxio.Encode_single_float64(params, -0.001, params.MaxLevel(), params.DefaultScale())
	//ptB := auxio.Encode_single_float64(params, 0.3, params.MaxLevel(), params.DefaultScale())
	ctA := encryptor.EncryptNew(ptA)
	//ctB := encryptor.EncryptNew(ptB)
	var ctSignA *rlwe.Ciphertext
	ctSignA, err = piecewisefunc.Sign_evaluate(params, rlk, ctA, piecewisefunc.Sign10)
	//ctSignA, err = piecewisefunc.ReLU_evaluate(params, rlk, ctA, piecewisefunc.Sign7)
	//ctSignA, err = piecewisefunc.Compare_evaluate(params, rlk, ctA, ctB, piecewisefunc.Sign7)

	if err != nil {
		panic(err)
	}
	auxio.Quick_check_vector(params, sk, ctSignA)

}
