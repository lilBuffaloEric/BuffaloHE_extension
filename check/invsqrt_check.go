package check

import (
	"fmt"
	auxio "project1-fhe_extension/auxiliary_io"
	invsqrt "project1-fhe_extension/inverse_sqrt"

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
		values[i] = 1024.0 / float64(params.Slots()) * (float64(i) + 0.000001)
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

	a, b := 0.0, 1024.0
	taylor_deg := 2

	fmt.Printf("Evaluation of the function Inverse Square Root for every slot in range (%0.2f, %0.2f] (degree for taylor initialization: %d)\n", a, b, taylor_deg)

	// Evaluation process
	fmt.Printf("We first use taylor expansion (deg: %d) at point (%0.2f+%0.2f)/2+1 to have a relative good guess:\n", taylor_deg, a, b)
	var ctGuess *rlwe.Ciphertext
	ctGuess, err = invsqrt.TaylorInitNew(params, rlk, ciphertext, invsqrt.Inv_sqrt_taylor2[:])
	if err != nil {
		panic(err)
	}

	fmt.Printf("Guessing Values:  ")
	auxio.Quick_check_vector(params, sk, ctGuess)
	fmt.Printf("%d levels has been consumed.\n", ciphertext.Level()-ctGuess.Level())

	IterNum := 2
	var ctApprox *rlwe.Ciphertext
	fmt.Printf("Then we will use Newton method to approximate the exact Inverse Square Root result, with %d iteration", IterNum)
	ctApprox, err = invsqrt.InvSqrtByNewton_dbg(params, rlk, sk, ciphertext, ctGuess, IterNum)
	if err != nil {
		panic(err)
	}

	fmt.Printf("Finally Approximate Result:  ")
	auxio.Quick_check_vector(params, sk, ctApprox)
	fmt.Printf("%d levels has been consumed. \n", ctGuess.Level()-ctApprox.Level())

	auxio.Print_vector_f64(values, len(values))

}
