package nonpolyfunc

import (
	"fmt"
	auxio "project1-fhe_extension/auxiliary_io"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

var (
	Inv_sqrt_taylor2 = [3]complex128{
		complex(0.082783231973869, 0), complex(-1.075804435198659e-04, 0), complex(6.291247875189161e-08, 0),
	}
	Inv_sqrt_taylor3 = [4]complex128{ // FIXME: Incorrect coefficient.
		complex(-1.021969094115431e-10, 0), complex(2.201936756316206e-07, 0), complex(-1.882657761597653e-04, 0), complex(0.096580437302848, 0),
	}
	Inv_sqrt_taylor7 = [8]complex128{ // FIXME: Incorrect coefficient.
		complex(1, 0), complex(1, 0), complex(1, 0), complex(1, 0), complex(1, 0), complex(1, 0), complex(1, 0), complex(1, 0),
	}
)

func TaylorInitNew(params ckks.Parameters, rlk *rlwe.RelinearizationKey, ctIn *rlwe.Ciphertext, apprxpoly []complex128) (ctOut *rlwe.Ciphertext, err error) {
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	poly := ckks.NewPoly(apprxpoly)
	ctOut, err = evaluator.EvaluatePoly(ctIn, poly, ctIn.Scale)
	return
}

// FIXME: Updtaed code, but hasn't been run yet.
func InvSqrtByNewton(params ckks.Parameters, rlk *rlwe.RelinearizationKey, ctIn *rlwe.Ciphertext, ctGuess *rlwe.Ciphertext, IterNum int) (ctOut *rlwe.Ciphertext, err error) {
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	fmt.Printf("Step 0.0: pre-compute pt_Int_3 = [3,3,...,3] and pt_Frac_1_2 = [1/2, 1/2,..., 1/2] \n")
	// ---------------- Acutal Operation ----------------
	pt_Int_3 := auxio.Encode_single_int(params, 3, params.MaxLevel())
	pt_Frac_1_2 := auxio.Encode_single_float64(params, 0.5, params.MaxLevel(), params.DefaultScale())
	auxio.Quick_decode_vector(params, pt_Int_3)    // debug: check the plaintext value
	auxio.Quick_decode_vector(params, pt_Frac_1_2) // debug: check the plaintext value
	// --------------------------------------------------
	ctGuess_duplicate := ctGuess.CopyNew()
	for i := 0; i < IterNum; i++ {
		fmt.Printf("Round %d approximating InvSqrt by Newton method...\n", i)
		fmt.Printf("Step %d.1: compute ctGuess_pow = ctGuess^2 \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess_pow2 := evaluator.MulRelinNew(ctGuess_duplicate, ctGuess_duplicate)
		err = evaluator.Rescale(ctGuess_pow2, ctGuess_duplicate.Scale, ctGuess_pow2)
		print(ctGuess_pow2.Level())
		if err != nil {
			ctOut = nil
			return
		}
		// auxio.Quick_check_vector(params, sk, ctGuess_pow2) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		fmt.Printf("Step %d.2: compute ctGuess_t_In = ctGuess * ctIn \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess_t_In := evaluator.MulRelinNew(ctGuess_duplicate, ctIn) // internally balances the level.
		print(ctGuess_t_In.Level())
		err = evaluator.Rescale(ctGuess_t_In, ctGuess_duplicate.Scale, ctGuess_t_In)
		print(ctGuess_t_In.Level())
		if err != nil {
			ctOut = nil
			return
		}
		// auxio.Quick_check_vector(params, sk, ctGuess_t_In) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		fmt.Printf("Step %d.3: compute ctGuessPow3_t_In = ctGuess_pow2 * ctGuess_t_In \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuessPow3_t_In := evaluator.MulRelinNew(ctGuess_pow2, ctGuess_t_In)
		err = evaluator.Rescale(ctGuessPow3_t_In, ctGuess_pow2.Scale, ctGuessPow3_t_In)
		print(ctGuessPow3_t_In.Level())
		if err != nil {
			ctOut = nil
			return
		}
		// auxio.Quick_check_vector(params, sk, ctGuessPow3_t_In) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		fmt.Printf("Step %d.4: compute ctGuess_t_3 = ctGuess * 3 \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess_t_3 := evaluator.MulNew(ctGuess_duplicate, pt_Int_3)
		// auxio.Quick_check_vector(params, sk, ctGuess_t_3) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		fmt.Printf("Step %d.5: compute ctGuess3_m_GuessPow3tIn = ctGuess_t_3 - ctGuessPow3_t_In \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess3_m_GuessPow3tIn := evaluator.SubNew(ctGuess_t_3, ctGuessPow3_t_In)
		// auxio.Quick_check_vector(params, sk, ctGuess3_m_GuessPow3tIn) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		fmt.Printf("Step %d.6: compute ctGuess3mGuessPow3tIn_t_1frac2 = 1/2 * ctGuess3_m_GuessPow3tIn \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess3mGuessPow3tIn_t_1frac2 := evaluator.MulNew(ctGuess3_m_GuessPow3tIn, pt_Frac_1_2)
		err = evaluator.Rescale(ctGuess3mGuessPow3tIn_t_1frac2, ctGuess3_m_GuessPow3tIn.Scale, ctGuess3mGuessPow3tIn_t_1frac2)
		if err != nil {
			ctOut = nil
			return
		}
		// auxio.Quick_check_vector(params, sk, ctGuess3mGuessPow3tIn_t_1frac2) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		fmt.Printf("Step %d.7: update ctGuess = ctGuess3mGuessPow3tIn_t_1frac2 \n", i)
		print(ctGuess3mGuessPow3tIn_t_1frac2.Level())
		ctGuess_duplicate = ctGuess3mGuessPow3tIn_t_1frac2.CopyNew()
		print(ctGuess_duplicate.Level())
		// auxio.Quick_check_vector(params, sk, ctGuess_duplicate) // debug: check the ciphertext's plaintext value
	}
	ctOut = ctGuess_duplicate.CopyNew()
	return
}

// Debug Version of InvSqrtByNewton(*)
func InvSqrtByNewton_dbg(params ckks.Parameters, rlk *rlwe.RelinearizationKey, sk *rlwe.SecretKey, ctIn *rlwe.Ciphertext, ctGuess *rlwe.Ciphertext, IterNum int) (ctOut *rlwe.Ciphertext, err error) {
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	fmt.Printf("Step 0.0: pre-compute pt_Int_3 = [3,3,...,3] and pt_Frac_1_2 = [1/2, 1/2,..., 1/2] \n")
	// ---------------- Acutal Operation ----------------
	pt_Int_3 := auxio.Encode_single_int(params, 3, params.MaxLevel())
	pt_Frac_1_2 := auxio.Encode_single_float64(params, 0.5, params.MaxLevel(), params.DefaultScale())
	auxio.Quick_decode_vector(params, pt_Int_3)    // debug: check the plaintext value
	auxio.Quick_decode_vector(params, pt_Frac_1_2) // debug: check the plaintext value
	// --------------------------------------------------
	ctGuess_duplicate := ctGuess.CopyNew()
	for i := 0; i < IterNum; i++ {
		fmt.Printf("Round %d approximating InvSqrt by Newton method...\n", i)
		fmt.Printf("Step %d.1: compute ctGuess_pow = ctGuess^2 \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess_pow2 := evaluator.MulRelinNew(ctGuess_duplicate, ctGuess_duplicate)
		err = evaluator.Rescale(ctGuess_pow2, ctGuess_duplicate.Scale, ctGuess_pow2)
		print(ctGuess_pow2.Level())
		if err != nil {
			ctOut = nil
			return
		}
		auxio.Quick_check_vector(params, sk, ctGuess_pow2) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		fmt.Printf("Step %d.2: compute ctGuess_t_In = ctGuess * ctIn \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess_t_In := evaluator.MulRelinNew(ctGuess_duplicate, ctIn) // internally balances the level.
		print(ctGuess_t_In.Level())
		err = evaluator.Rescale(ctGuess_t_In, ctGuess_duplicate.Scale, ctGuess_t_In)
		print(ctGuess_t_In.Level())
		if err != nil {
			ctOut = nil
			return
		}
		auxio.Quick_check_vector(params, sk, ctGuess_t_In) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		fmt.Printf("Step %d.3: compute ctGuessPow3_t_In = ctGuess_pow2 * ctGuess_t_In \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuessPow3_t_In := evaluator.MulRelinNew(ctGuess_pow2, ctGuess_t_In)
		err = evaluator.Rescale(ctGuessPow3_t_In, ctGuess_pow2.Scale, ctGuessPow3_t_In)
		print(ctGuessPow3_t_In.Level())
		if err != nil {
			ctOut = nil
			return
		}
		auxio.Quick_check_vector(params, sk, ctGuessPow3_t_In) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		fmt.Printf("Step %d.4: compute ctGuess_t_3 = ctGuess * 3 \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess_t_3 := evaluator.MulNew(ctGuess_duplicate, pt_Int_3)
		auxio.Quick_check_vector(params, sk, ctGuess_t_3) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		fmt.Printf("Step %d.5: compute ctGuess3_m_GuessPow3tIn = ctGuess_t_3 - ctGuessPow3_t_In \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess3_m_GuessPow3tIn := evaluator.SubNew(ctGuess_t_3, ctGuessPow3_t_In)
		auxio.Quick_check_vector(params, sk, ctGuess3_m_GuessPow3tIn) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		fmt.Printf("Step %d.6: compute ctGuess3mGuessPow3tIn_t_1frac2 = 1/2 * ctGuess3_m_GuessPow3tIn \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess3mGuessPow3tIn_t_1frac2 := evaluator.MulNew(ctGuess3_m_GuessPow3tIn, pt_Frac_1_2)
		err = evaluator.Rescale(ctGuess3mGuessPow3tIn_t_1frac2, ctGuess3_m_GuessPow3tIn.Scale, ctGuess3mGuessPow3tIn_t_1frac2)
		if err != nil {
			ctOut = nil
			return
		}
		auxio.Quick_check_vector(params, sk, ctGuess3mGuessPow3tIn_t_1frac2) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		fmt.Printf("Step %d.7: update ctGuess = ctGuess3mGuessPow3tIn_t_1frac2 \n", i)
		print(ctGuess3mGuessPow3tIn_t_1frac2.Level())
		ctGuess_duplicate = ctGuess3mGuessPow3tIn_t_1frac2.CopyNew()
		print(ctGuess_duplicate.Level())
		auxio.Quick_check_vector(params, sk, ctGuess_duplicate) // debug: check the ciphertext's plaintext value
	}
	ctOut = ctGuess_duplicate.CopyNew()
	return
}
