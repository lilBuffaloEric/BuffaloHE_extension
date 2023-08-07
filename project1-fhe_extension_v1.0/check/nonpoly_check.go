package check

/*
 * Check for funcs in nonpoly_func
 */

import (
	"fmt"

	auxio "project1-fhe_extension_v1.0/auxiliary_io"
	nonpolyfunc "project1-fhe_extension_v1.0/nonpoly_func"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func Invsqrt_check() {
	params, err := ckks.NewParametersFromLiteral(ckks.PN15QP880)
	if err != nil {
		panic(err)
	}

	encoder := ckks.NewEncoder(params)

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

	// Values to encrypt
	values := make([]float64, params.Slots())
	for i := range values {
		values[i] = float64(1<<21)/float64(params.Slots())*float64(i) + 0.00000001
	}

	fmt.Printf("CKKS parameters: logN = %d, logQ = %d, levels = %d, scale= %f, sigma = %f \n",
		params.LogN(), params.LogQP(), params.MaxLevel()+1, params.DefaultScale().Float64(), params.Sigma())
	fmt.Println()
	fmt.Print("Target Values    : ")
	auxio.Print_vector_f64(values, len(values))
	// Plaintext creation and encoding process
	plaintext := encoder.EncodeNew(values, params.MaxLevel(), params.DefaultScale(), params.LogSlots())

	// Encryption process
	ciphertext := encryptor.EncryptNew(plaintext)

	a, b := 0.00000001, float64(1<<21)
	taylor_deg := 2

	fmt.Printf("Evaluation of the function Inverse Square Root for every slot in range (%0.2f, %0.2f] (degree for taylor initialization: %d)\n", a, b, taylor_deg)

	// Evaluation process
	fmt.Printf("We first use taylor expansion (deg: %d) at point (%0.2f+%0.2f)/2+1 to have a relative good guess:\n", taylor_deg, a, b)
	var ctGuess *rlwe.Ciphertext
	ctGuess, err = nonpolyfunc.TaylorInitNew(params, rlk, ciphertext, nonpolyfunc.Inv_sqrt_taylor1_0to2pow21[:])
	if err != nil {
		panic(err)
	}

	fmt.Printf("Guessing Values:  ")
	auxio.Quick_check_vector(params, sk, ctGuess)
	fmt.Printf("%d levels has been consumed.\n", ciphertext.Level()-ctGuess.Level())

	IterNum := 5
	var ctApprox *rlwe.Ciphertext
	fmt.Printf("Then we will use Newton method to approximate the exact Inverse Square Root result, with %d iteration", IterNum)
	ctApprox, err = nonpolyfunc.InvSqrtByNewton_dbg(params, rlk, pk, sk, ciphertext, ctGuess, IterNum)
	if err != nil {
		panic(err)
	}

	fmt.Printf("Finally Approximate Result:  ")
	auxio.Quick_check_vector(params, sk, ctApprox)
	fmt.Printf("%d levels has been consumed. \n", ctGuess.Level()-ctApprox.Level())

	auxio.Print_vector_f64(values, len(values))

}

func Softmax_check() {
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

	ptA := auxio.Encode_single_float64(params, 2, params.MaxLevel(), params.DefaultScale())
	//ptB := auxio.Encode_single_float64(params, 0.3, params.MaxLevel(), params.DefaultScale())
	ctA := encryptor.EncryptNew(ptA)
	//ctB := encryptor.EncryptNew(ptB)
	var ctInvA *rlwe.Ciphertext
	ctInvA, err = nonpolyfunc.InvByGoldSchmidt(params, sk, pk, rlk, ctA, 3, 4)
	if err != nil {
		panic(err)
	}
	auxio.Quick_check_vector(params, sk, ctInvA)

}
