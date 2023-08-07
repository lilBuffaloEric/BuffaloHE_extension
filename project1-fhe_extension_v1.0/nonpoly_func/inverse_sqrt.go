package nonpolyfunc

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v4/ckks/bootstrapping"
	auxio "project1-fhe_extension_v1.0/auxiliary_io"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

var (
	Inv_sqrt_taylor2_0to2pow10 = [3]complex128{ // Interval [0.001,2^10] with order truncated to 2
		complex(0.082783231973869, 0), complex(-1.075804435198659e-04, 0), complex(6.291247875189161e-08, 0),
	}
	Inv_sqrt_taylor2_0to2pow21 = [3]complex128{ // Interval [0.1^7, 2^21] with order truncated to 2
		complex(0.001831053814386, 0), complex(-1.164151552936788e-09, 0), complex(3.330661132954803e-16, 0),
	}
	Inv_sqrt_taylor1_0to2pow21 = [3]complex128{
		complex(0.001464843051509, 0), complex(-4.656606211747153e-10, 0), complex(0, 0),
	}
	Inv_sqrt_taypor1_0to2pow16 = [3]complex128{ // Interval [0.1^7, 2^16] with order truncated to 1, a 15 round iteration of Lazy Newton method is suggested.
		complex(0.008286281154378, 0), complex(-8.428983850973325e-08, 0), complex(0, 0),
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
	// fmt.Printf("Step 0.0: pre-compute pt_Int_3 = [3,3,...,3] and pt_Frac_1_2 = [1/2, 1/2,..., 1/2] \n")
	// ---------------- Acutal Operation ----------------
	// pt_Int_3 := auxio.Encode_single_int(params, 3, params.MaxLevel()) // option 1
	// pt_Frac_1_2 := auxio.Encode_single_float64(params, 0.5, params.MaxLevel(), params.DefaultScale()) // option 1
	// auxio.Quick_decode_vector(params, pt_Int_3)    // debug: check the plaintext value // option 1
	// auxio.Quick_decode_vector(params, pt_Frac_1_2) // debug: check the plaintext value // option 1
	// --------------------------------------------------
	ctGuess_duplicate := ctGuess.CopyNew()
	for i := 0; i < IterNum; i++ {
		// fmt.Printf("Round %d approximating InvSqrt by Newton method...\n", i)
		// fmt.Printf("Step %d.1: compute ctGuess_pow = ctGuess^2 \n", i)
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
		// fmt.Printf("Step %d.2: compute ctGuess_t_In = ctGuess * ctIn \n", i)
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
		// fmt.Printf("Step %d.3: compute ctGuessPow3_t_In = ctGuess_pow2 * ctGuess_t_In \n", i)
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
		// fmt.Printf("Step %d.4: compute ctGuess_t_3 = ctGuess * 3 \n", i)
		// ---------------- Acutal Operation ----------------
		scale := rlwe.NewScale(ctGuess_duplicate.Scale.Float64() * ctGuess_t_In.Scale.Float64()) // have to synchronize the scale.
		pt_Int_3 := auxio.Encode_single_float64(params, 3.0, ctGuess.Level(), scale)
		ctGuess_t_3 := evaluator.MulNew(ctGuess_duplicate, pt_Int_3)
		err = evaluator.Rescale(ctGuess_t_3, ctGuessPow3_t_In.GetScale(), ctGuess_t_3)
		if err != nil {
			panic(err)
		}
		// auxio.Quick_check_vector(params, sk, ctGuess_t_3) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		// fmt.Printf("Step %d.5: compute ctGuess3_m_GuessPow3tIn = ctGuess_t_3 - ctGuessPow3_t_In \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess3_m_GuessPow3tIn := evaluator.SubNew(ctGuess_t_3, ctGuessPow3_t_In)
		// auxio.Quick_check_vector(params, sk, ctGuess3_m_GuessPow3tIn) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		// fmt.Printf("Step %d.6: compute ctGuess3mGuessPow3tIn_t_1frac2 = 1/2 * ctGuess3_m_GuessPow3tIn \n", i)
		// ---------------- Acutal Operation ----------------
		// ctGuess3mGuessPow3tIn_t_1frac2 := evaluator.MulNew(ctGuess3_m_GuessPow3tIn, pt_Frac_1_2) // option 1
		ctGuess3mGuessPow3tIn_t_1frac2 := evaluator.MultByConstNew(ctGuess3_m_GuessPow3tIn, 0.5) // option 2
		err = evaluator.Rescale(ctGuess3mGuessPow3tIn_t_1frac2, ctGuess3_m_GuessPow3tIn.Scale, ctGuess3mGuessPow3tIn_t_1frac2)
		if err != nil {
			ctOut = nil
			return
		}
		// auxio.Quick_check_vector(params, sk, ctGuess3mGuessPow3tIn_t_1frac2) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		// fmt.Printf("Step %d.7: update ctGuess = ctGuess3mGuessPow3tIn_t_1frac2 \n", i)
		// print(ctGuess3mGuessPow3tIn_t_1frac2.Level())
		ctGuess_duplicate = ctGuess3mGuessPow3tIn_t_1frac2.CopyNew()
		// print(ctGuess_duplicate.Level())
		// auxio.Quick_check_vector(params, sk, ctGuess_duplicate) // debug: check the ciphertext's plaintext value
	}
	ctOut = ctGuess_duplicate.CopyNew()
	return
}

// Debug Version of InvSqrtByNewton(*)
// this function is only for experiment purpose, since it will internally decrypt a ciphertext
func InvSqrtByNewton_dbg(params ckks.Parameters, rlk *rlwe.RelinearizationKey, pk *rlwe.PublicKey, sk *rlwe.SecretKey, ctIn *rlwe.Ciphertext, ctGuess *rlwe.Ciphertext, IterNum int) (ctOut *rlwe.Ciphertext, err error) {
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	encryptor := ckks.NewEncryptor(params, pk)
	decryptor := ckks.NewDecryptor(params, sk)
	encoder := ckks.NewEncoder(params)
	//fmt.Printf("Step 0.0: pre-compute pt_Int_3 = [3,3,...,3] and pt_Frac_1_2 = [1/2, 1/2,..., 1/2] \n")
	// ---------------- Acutal Operation ----------------
	// pt_Int_3 := auxio.Encode_single_int(params, 3, params.MaxLevel()) // option 1
	// pt_Frac_1_2 := auxio.Encode_single_float64(params, 0.5, params.MaxLevel(), params.DefaultScale()) // option 1
	// auxio.Quick_decode_vector(params, pt_Int_3)    // debug: check the plaintext value // option 1
	// auxio.Quick_decode_vector(params, pt_Frac_1_2) // debug: check the plaintext value // option 1
	// --------------------------------------------------
	ctGuess_duplicate := ctGuess.CopyNew()
	for i := 0; i < IterNum; i++ {
		if ctGuess_duplicate.Level() < 3 {
			ctGuess_duplicate = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctGuess_duplicate), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
		}
		//fmt.Printf("Round %d approximating InvSqrt by Newton method...\n", i)
		//fmt.Printf("Step %d.1: compute ctGuess_pow = ctGuess^2 \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess_pow2 := evaluator.MulRelinNew(ctGuess_duplicate, ctGuess_duplicate)
		err = evaluator.Rescale(ctGuess_pow2, ctGuess_duplicate.Scale, ctGuess_pow2)
		// print(ctGuess_pow2.Level())
		if err != nil {
			ctOut = nil
			return
		}
		//auxio.Quick_check_vector(params, sk, ctGuess_pow2) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		//fmt.Printf("Step %d.2: compute ctGuess_t_In = ctGuess * ctIn \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess_t_In := evaluator.MulRelinNew(ctGuess_duplicate, ctIn) // internally balances the level.
		//print(ctGuess_t_In.Level())
		err = evaluator.Rescale(ctGuess_t_In, ctGuess_duplicate.Scale, ctGuess_t_In)
		//print(ctGuess_t_In.Level())
		if err != nil {
			ctOut = nil
			return
		}
		//auxio.Quick_check_vector(params, sk, ctGuess_t_In) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		//fmt.Printf("Step %d.3: compute ctGuessPow3_t_In = ctGuess_pow2 * ctGuess_t_In \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuessPow3_t_In := evaluator.MulRelinNew(ctGuess_pow2, ctGuess_t_In)
		err = evaluator.Rescale(ctGuessPow3_t_In, ctGuess_pow2.Scale, ctGuessPow3_t_In)
		//print(ctGuessPow3_t_In.Level())
		if err != nil {
			ctOut = nil
			return
		}
		//auxio.Quick_check_vector(params, sk, ctGuessPow3_t_In) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		// fmt.Printf("Step %d.4: compute ctGuess_t_3 = ctGuess * 3 \n", i)
		// ---------------- Acutal Operation ----------------
		scale := rlwe.NewScale(ctGuess_duplicate.Scale.Float64() * ctGuess_t_In.Scale.Float64()) // have to synchronize the scale.
		pt_Int_3 := auxio.Encode_single_float64(params, 3.0, ctGuess.Level(), scale)
		ctGuess_t_3 := evaluator.MulNew(ctGuess_duplicate, pt_Int_3)
		err = evaluator.Rescale(ctGuess_t_3, ctGuessPow3_t_In.GetScale(), ctGuess_t_3)
		if err != nil {
			panic(err)
		}
		//auxio.Quick_check_vector(params, sk, ctGuess_t_3) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		//fmt.Printf("Step %d.5: compute ctGuess3_m_GuessPow3tIn = ctGuess_t_3 - ctGuessPow3_t_In \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess3_m_GuessPow3tIn := evaluator.SubNew(ctGuess_t_3, ctGuessPow3_t_In)
		//auxio.Quick_check_vector(params, sk, ctGuess3_m_GuessPow3tIn) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		//fmt.Printf("Step %d.6: compute ctGuess3mGuessPow3tIn_t_1frac2 = 1/2 * ctGuess3_m_GuessPow3tIn \n", i)
		// ---------------- Acutal Operation ----------------
		// ctGuess3mGuessPow3tIn_t_1frac2 := evaluator.MulNew(ctGuess3_m_GuessPow3tIn, pt_Frac_1_2) // option 1
		ctGuess3mGuessPow3tIn_t_1frac2 := evaluator.MultByConstNew(ctGuess3_m_GuessPow3tIn, 0.5) // option 2
		err = evaluator.Rescale(ctGuess3mGuessPow3tIn_t_1frac2, ctGuess3_m_GuessPow3tIn.Scale, ctGuess3mGuessPow3tIn_t_1frac2)
		if err != nil {
			ctOut = nil
			return
		}
		//auxio.Quick_check_vector(params, sk, ctGuess3mGuessPow3tIn_t_1frac2) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		//fmt.Printf("Step %d.7: update ctGuess = ctGuess3mGuessPow3tIn_t_1frac2 with level %d \n", i, ctGuess3mGuessPow3tIn_t_1frac2.Level())
		ctGuess_duplicate = ctGuess3mGuessPow3tIn_t_1frac2.CopyNew()
		// print(ctGuess_duplicate.Level())
		//auxio.Quick_check_vector(params, sk, ctGuess_duplicate) // debug: check the ciphertext's plaintext value
	}
	ctOut = ctGuess_duplicate.CopyNew()
	return
}

// Debug Version of InvSqrtByNewton(*)
// this function is only for experiment purpose, since it will internally decrypt a ciphertext
func InvSqrtByNewton_btpVer_dbg(params ckks.Parameters, rlk *rlwe.RelinearizationKey, pk *rlwe.PublicKey, sk *rlwe.SecretKey, ctIn *rlwe.Ciphertext, ctGuess *rlwe.Ciphertext, IterNum int, maxLevel int, reencryptionMode bool, Boostrapper ...*bootstrapping.Bootstrapper) (ctOut *rlwe.Ciphertext, btptimes int, err error) {
	d := 128 // only for dbg usage.
	encryptor := ckks.NewEncryptor(params, pk)
	decryptor := ckks.NewDecryptor(params, sk)
	encoder := ckks.NewEncoder(params)
	ctGuess_duplicate := ctGuess.CopyNew()
	var ctIn_duplicate *rlwe.Ciphertext
	var boostrapper *bootstrapping.Bootstrapper
	if len(Boostrapper) != 0 {
		boostrapper = Boostrapper[0]
	}
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	//fmt.Printf("Step 0.0: pre-compute pt_Int_3 = [3,3,...,3] and pt_Frac_1_2 = [1/2, 1/2,..., 1/2] \n")
	// ---------------- Acutal Operation ----------------
	// pt_Int_3 := auxio.Encode_single_int(params, 3, params.MaxLevel()) // option 1
	// pt_Frac_1_2 := auxio.Encode_single_float64(params, 0.5, params.MaxLevel(), params.DefaultScale()) // option 1
	// auxio.Quick_decode_vector(params, pt_Int_3)    // debug: check the plaintext value // option 1
	// auxio.Quick_decode_vector(params, pt_Frac_1_2) // debug: check the plaintext value // option 1
	// --------------------------------------------------
	if ctIn.Level() < maxLevel {
		if reencryptionMode {
			ctIn_duplicate = boostrapper.Bootstrap(ctIn)
		} else {
			ctIn_duplicate = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctIn), params.LogSlots()), maxLevel, params.DefaultScale(), params.LogSlots()))
		}
		btptimes++
	} else {
		ctIn_duplicate = ctIn.CopyNew()
	}
	for i := 0; i < IterNum; i++ {
		if ctGuess_duplicate.Level() < 3 {
			if reencryptionMode {
				ctGuess_duplicate = boostrapper.Bootstrap(ctGuess_duplicate)
			} else {
				ctGuess_duplicate = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctGuess_duplicate), params.LogSlots()), maxLevel, params.DefaultScale(), params.LogSlots()))
			}
			btptimes++
		}
		// fmt.Printf("Round %d approximating InvSqrt by Newton method...\n", i)
		fmt.Printf("ctGuess_duplicate has scale %f, level %d\n", math.Log2(ctGuess_duplicate.Scale.Float64()), ctGuess_duplicate.Level())
		auxio.Quick_check_matrix(params, sk, ctGuess_duplicate, d, d)
		//fmt.Printf("Step %d.1: compute ctGuess_pow = ctGuess^2 \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess_pow2 := evaluator.MulRelinNew(ctGuess_duplicate, ctGuess_duplicate)
		err = evaluator.Rescale(ctGuess_pow2, ctGuess_duplicate.Scale, ctGuess_pow2)
		fmt.Printf("ctGuess_pow2 has scale %f, level %d\n", math.Log2(ctGuess_pow2.Scale.Float64()), ctGuess_pow2.Level())
		auxio.Quick_check_matrix(params, sk, ctGuess_pow2, d, d)
		// print(ctGuess_pow2.Level())
		if err != nil {
			ctOut = nil
			return
		}
		//auxio.Quick_check_vector(params, sk, ctGuess_pow2) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		//fmt.Printf("Step %d.2: compute ctGuess_t_In = ctGuess * ctIn \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess_t_In := evaluator.MulRelinNew(ctGuess_duplicate, ctIn_duplicate) // internally balances the level.
		//print(ctGuess_t_In.Level())
		err = evaluator.Rescale(ctGuess_t_In, ctGuess_duplicate.Scale, ctGuess_t_In)
		fmt.Printf("ctGuess_t_In has scale %f, level %d\n", math.Log2(ctGuess_t_In.Scale.Float64()), ctGuess_t_In.Level())
		auxio.Quick_check_matrix(params, sk, ctGuess_t_In, d, d)
		//print(ctGuess_t_In.Level())
		if err != nil {
			ctOut = nil
			return
		}
		//auxio.Quick_check_vector(params, sk, ctGuess_t_In) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		//fmt.Printf("Step %d.3: compute ctGuessPow3_t_In = ctGuess_pow2 * ctGuess_t_In \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuessPow3_t_In := evaluator.MulRelinNew(ctGuess_pow2, ctGuess_t_In)
		err = evaluator.Rescale(ctGuessPow3_t_In, ctGuess_pow2.Scale, ctGuessPow3_t_In)
		fmt.Printf("ctGuessPow3_t_In has scale %f, level %d\n", math.Log2(ctGuessPow3_t_In.Scale.Float64()), ctGuessPow3_t_In.Level())
		auxio.Quick_check_matrix(params, sk, ctGuessPow3_t_In, d, d)
		//print(ctGuessPow3_t_In.Level())
		if err != nil {
			ctOut = nil
			return
		}
		//auxio.Quick_check_vector(params, sk, ctGuessPow3_t_In) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		// fmt.Printf("Step %d.4: compute ctGuess_t_3 = ctGuess * 3 \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess_t_3 := evaluator.MultByConstNew(ctGuess_duplicate, 3.0)
		evaluator.SetScale(ctGuess_t_3, ctGuessPow3_t_In.Scale)
		/*
			scale := rlwe.NewScale(ctGuess_duplicate.Scale.Float64() * ctGuess_t_In.Scale.Float64()) // have to synchronize the scale.
			pt_Int_3 := auxio.Encode_single_float64(params, 3.0, ctGuess_duplicate.Level(), scale)
			ctGuess_t_3 := evaluator.MulNew(ctGuess_duplicate, pt_Int_3)
			err = evaluator.Rescale(ctGuess_t_3, ctGuessPow3_t_In.GetScale(), ctGuess_t_3)
			fmt.Printf("ctGuess_t_3 has scale %f, level %d\n", math.Log2(ctGuess_t_3.Scale.Float64()), ctGuess_t_3.Level())
			auxio.Quick_check_matrix(params, sk, ctGuess_t_3, d, d)
			if err != nil {
				panic(err)
			}
		*/
		fmt.Printf("ctGuess_t_3 has scale %f, level %d\n", math.Log2(ctGuess_t_3.Scale.Float64()), ctGuess_t_3.Level())
		auxio.Quick_check_matrix(params, sk, ctGuess_t_3, d, d) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		//fmt.Printf("Step %d.5: compute ctGuess3_m_GuessPow3tIn = ctGuess_t_3 - ctGuessPow3_t_In \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess3_m_GuessPow3tIn := evaluator.SubNew(ctGuess_t_3, ctGuessPow3_t_In)
		fmt.Printf("ctGuess3_m_GuessPow3tIn has scale %f, level %d\n", math.Log2(ctGuess3_m_GuessPow3tIn.Scale.Float64()), ctGuess3_m_GuessPow3tIn.Level())
		auxio.Quick_check_matrix(params, sk, ctGuess3_m_GuessPow3tIn, d, d)
		//auxio.Quick_check_vector(params, sk, ctGuess3_m_GuessPow3tIn) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		//fmt.Printf("Step %d.6: compute ctGuess3mGuessPow3tIn_t_1frac2 = 1/2 * ctGuess3_m_GuessPow3tIn \n", i)
		// ---------------- Acutal Operation ----------------
		// ctGuess3mGuessPow3tIn_t_1frac2 := evaluator.MulNew(ctGuess3_m_GuessPow3tIn, pt_Frac_1_2) // option 1
		var ctGuess3mGuessPow3tIn_t_1frac2 *rlwe.Ciphertext
		if ctGuess3_m_GuessPow3tIn.Level() == 1 {
			scale2 := rlwe.NewScale(float64(params.RingQ().Modulus[ctGuess3_m_GuessPow3tIn.Level()]))
			scale3 := ctGuess3_m_GuessPow3tIn.Scale
			pt_1frac2 := auxio.Encode_single_float64(params, 0.5, ctGuess3_m_GuessPow3tIn.Level(), scale2.Mul(params.DefaultScale()).Div(scale3))
			ctGuess3mGuessPow3tIn_t_1frac2 = evaluator.MulNew(ctGuess3_m_GuessPow3tIn, pt_1frac2)
			err = evaluator.Rescale(ctGuess3mGuessPow3tIn_t_1frac2, params.DefaultScale(), ctGuess3mGuessPow3tIn_t_1frac2)
			if reencryptionMode {
				ctGuess3mGuessPow3tIn_t_1frac2 = boostrapper.Bootstrap(ctGuess3mGuessPow3tIn_t_1frac2)
			} else {
				ctGuess3mGuessPow3tIn_t_1frac2 = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctGuess3mGuessPow3tIn_t_1frac2), params.LogSlots()), params.MaxLevel(), params.DefaultScale(), params.LogSlots()))
			}
			btptimes++
		} else {
			ctGuess3mGuessPow3tIn_t_1frac2 = evaluator.MultByConstNew(ctGuess3_m_GuessPow3tIn, 0.5) // option 2
			err = evaluator.Rescale(ctGuess3mGuessPow3tIn_t_1frac2, ctGuess3_m_GuessPow3tIn.Scale, ctGuess3mGuessPow3tIn_t_1frac2)
		}
		fmt.Printf("ctGuess3mGuessPow3tIn_t_1frac2 has scale %f, level %d\n", math.Log2(ctGuess3mGuessPow3tIn_t_1frac2.Scale.Float64()), ctGuess3mGuessPow3tIn_t_1frac2.Level())
		auxio.Quick_check_matrix(params, sk, ctGuess3mGuessPow3tIn_t_1frac2, d, d)
		if err != nil {
			ctOut = nil
			return
		}
		//auxio.Quick_check_vector(params, sk, ctGuess3mGuessPow3tIn_t_1frac2) // debug: check the ciphertext's plaintext value
		// --------------------------------------------------
		//fmt.Printf("Step %d.7: update ctGuess = ctGuess3mGuessPow3tIn_t_1frac2 with level %d \n", i, ctGuess3mGuessPow3tIn_t_1frac2.Level())
		ctGuess_duplicate = ctGuess3mGuessPow3tIn_t_1frac2.CopyNew()
		// print(ctGuess_duplicate.Level())
		//auxio.Quick_check_vector(params, sk, ctGuess_duplicate) // debug: check the ciphertext's plaintext value
	}
	ctOut = ctGuess_duplicate.CopyNew()
	return
}

// Debug Version of InvSqrtByNewton(*)
// this function is only for experiment purpose, since it will internally decrypt a ciphertext
func InvSqrtByNewton_btpVer2_dbg(params ckks.Parameters, rlk *rlwe.RelinearizationKey, pk *rlwe.PublicKey, sk *rlwe.SecretKey, ctIn *rlwe.Ciphertext, ctGuess *rlwe.Ciphertext, IterNum int, maxLevel int, reencryptionMode bool, Boostrapper ...*bootstrapping.Bootstrapper) (ctOut *rlwe.Ciphertext, btptimes int, err error) {
	//d := 128 // only for dbg usage.
	encryptor := ckks.NewEncryptor(params, pk)
	decryptor := ckks.NewDecryptor(params, sk)
	encoder := ckks.NewEncoder(params)
	ctGuess_duplicate := ctGuess.CopyNew()
	var ctIn_duplicate *rlwe.Ciphertext
	var boostrapper *bootstrapping.Bootstrapper
	if len(Boostrapper) != 0 {
		boostrapper = Boostrapper[0]
	}
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	//fmt.Printf("Step 0.0: pre-compute pt_Int_3 = [3,3,...,3] and pt_Frac_1_2 = [1/2, 1/2,..., 1/2] \n")
	// ---------------- Acutal Operation ----------------
	// pt_Int_3 := auxio.Encode_single_int(params, 3, params.MaxLevel()) // option 1
	// pt_Frac_1_2 := auxio.Encode_single_float64(params, 0.5, params.MaxLevel(), params.DefaultScale()) // option 1
	// auxio.Quick_decode_vector(params, pt_Int_3)    // debug: check the plaintext value // option 1
	// auxio.Quick_decode_vector(params, pt_Frac_1_2) // debug: check the plaintext value // option 1
	// --------------------------------------------------
	if ctIn.Level() < maxLevel {
		if reencryptionMode {
			ctIn_duplicate = boostrapper.Bootstrap(ctIn)
		} else {
			ctIn_duplicate = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctIn), params.LogSlots()), maxLevel, params.DefaultScale(), params.LogSlots()))
		}
		btptimes++
	} else {
		ctIn_duplicate = ctIn.CopyNew()
	}
	for i := 0; i < IterNum; i++ {
		if ctGuess_duplicate.Level() < 3 {
			if reencryptionMode {
				ctGuess_duplicate = boostrapper.Bootstrap(ctGuess_duplicate)
			} else {
				ctGuess_duplicate = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctGuess_duplicate), params.LogSlots()), maxLevel, params.DefaultScale(), params.LogSlots()))
			}
			btptimes++
		}
		// fmt.Printf("Round %d approximating InvSqrt by Newton method...\n", i)
		//fmt.Printf("ctGuess_duplicate has scale %f, level %d\n", math.Log2(ctGuess_duplicate.Scale.Float64()), ctGuess_duplicate.Level())
		// auxio.Quick_check_matrix(params, sk, ctGuess_duplicate, d, d)
		//fmt.Printf("Step %d.1: compute ctGuess_pow = ctGuess^2 \n", i)
		// ---------------- Acutal Operation ----------------
		ctGuess_pow2 := evaluator.MulRelinNew(ctGuess_duplicate, ctGuess_duplicate)
		err = evaluator.Rescale(ctGuess_pow2, params.DefaultScale(), ctGuess_pow2)
		//fmt.Printf("ctGuess_pow2 has scale %f, level %d\n", math.Log2(ctGuess_pow2.Scale.Float64()), ctGuess_pow2.Level())
		//auxio.Quick_check_matrix(params, sk, ctGuess_pow2, d, d)
		if err != nil {
			ctOut = nil
			return
		}

		ctGuess1frac2 := evaluator.MultByConstNew(ctGuess_duplicate, -0.5)
		err = evaluator.Rescale(ctGuess1frac2, params.DefaultScale(), ctGuess1frac2)
		//fmt.Printf("ctGuess1frac2 has scale %f, level %d\n", math.Log2(ctGuess1frac2.Scale.Float64()), ctGuess1frac2.Level())
		//auxio.Quick_check_matrix(params, sk, ctGuess1frac2, d, d)
		if err != nil {
			ctOut = nil
			return
		}

		ctGuess1frac2tIn := evaluator.MulRelinNew(ctGuess1frac2, ctIn_duplicate)
		err = evaluator.Rescale(ctGuess1frac2tIn, params.DefaultScale(), ctGuess1frac2tIn)
		//fmt.Printf("ctGuess1frac2tIn has scale %f, level %d\n", math.Log2(ctGuess1frac2tIn.Scale.Float64()), ctGuess1frac2tIn.Level())
		//auxio.Quick_check_matrix(params, sk, ctGuess1frac2tIn, d, d)
		if err != nil {
			ctOut = nil
			return
		}

		evaluator.SetScale(ctGuess_pow2, params.DefaultScale().Mul(rlwe.NewScale(float64(params.RingQ().Modulus[ctGuess1frac2tIn.Level()]))).Div(ctGuess1frac2tIn.Scale))

		ctGuessPow31frac2tIn := evaluator.MulRelinNew(ctGuess1frac2tIn, ctGuess_pow2)
		err = evaluator.Rescale(ctGuessPow31frac2tIn, params.DefaultScale(), ctGuessPow31frac2tIn)
		fmt.Printf("ctGuessPow31frac2tIn has scale %f, level %d\n", math.Log2(ctGuessPow31frac2tIn.Scale.Float64()), ctGuessPow31frac2tIn.Level())
		//auxio.Quick_check_matrix(params, sk, ctGuessPow31frac2tIn, d, d)
		if err != nil {
			ctOut = nil
			return
		}

		ctGuess3frac2 := evaluator.MultByConstNew(ctGuess_duplicate, 1.5)
		err = evaluator.Rescale(ctGuess3frac2, params.DefaultScale(), ctGuess3frac2)
		evaluator.SetScale(ctGuess3frac2, ctGuessPow31frac2tIn.Scale)
		//fmt.Printf("ctGuess3frac2 has scale %f, level %d\n", math.Log2(ctGuess3frac2.Scale.Float64()), ctGuess3frac2.Level())
		//auxio.Quick_check_matrix(params, sk, ctGuess3frac2, d, d)
		if err != nil {
			ctOut = nil
			return
		}

		ctGuess_duplicate = evaluator.AddNew(ctGuess3frac2, ctGuessPow31frac2tIn)
	}
	ctOut = ctGuess_duplicate.CopyNew()
	return
}

func InvSqrtByNewton_btpVer3_dbg(params ckks.Parameters, rlk *rlwe.RelinearizationKey, pk *rlwe.PublicKey, sk *rlwe.SecretKey, ctIn *rlwe.Ciphertext, ctGuess *rlwe.Ciphertext, IterNum int, maxLevel int, reencryptionMode bool, Boostrapper ...*bootstrapping.Bootstrapper) (ctOut *rlwe.Ciphertext, btptimes int, err error) {
	// d := 128 // only for dbg
	encryptor := ckks.NewEncryptor(params, pk)
	decryptor := ckks.NewDecryptor(params, sk)
	encoder := ckks.NewEncoder(params)
	cty := ctGuess.CopyNew()
	var ctIn_dup *rlwe.Ciphertext
	var boostrapper *bootstrapping.Bootstrapper
	if len(Boostrapper) != 0 {
		boostrapper = Boostrapper[0]
	}
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	if ctIn.Level() < maxLevel {
		if reencryptionMode {
			ctIn_dup = boostrapper.Bootstrap(ctIn)
		} else {
			ctIn_dup = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctIn), params.LogSlots()), maxLevel, params.DefaultScale(), params.LogSlots()))
		}
		btptimes++
	} else {
		ctIn_dup = ctIn.CopyNew()
	}
	/*
		ctIn_pow2 := evaluator.MulRelinNew(ctIn_dup, ctIn_dup)
		err = evaluator.Rescale(ctIn_pow2, params.DefaultScale(), ctIn_pow2)
		if err != nil {
			ctOut = nil
			return
		}
		auxio.Quick_check_infos(ctIn_pow2, "ctIn_pow2")
		auxio.Quick_check_matrix(params, sk, ctIn_pow2, d, d)

		ctIn_pow4 := evaluator.MulRelinNew(ctIn_pow2, ctIn_pow2)
		err = evaluator.Rescale(ctIn_pow4, params.DefaultScale(), ctIn_pow4)
		if err != nil {
			ctOut = nil
			return
		}
		auxio.Quick_check_infos(ctIn_pow4, "ctIn_pow4")
		auxio.Quick_check_matrix(params, sk, ctIn_pow4, d, d)

		scale_InMod := make([]rlwe.Scale, 3)
		scale_InMod[0] = rlwe.NewScale(float64(params.RingQ().Modulus[ctIn_dup.Level()]))
		scale_InMod[1] = rlwe.NewScale(float64(params.RingQ().Modulus[ctIn_pow2.Level()]))
		scale_InMod[2] = rlwe.NewScale(float64(params.RingQ().Modulus[ctIn_pow4.Level()]))

		scale_In := make([]rlwe.Scale, 3)
		scale_In[0] = ctIn_dup.GetScale()
		scale_In[1] = ctIn_pow2.GetScale()
		scale_In[2] = ctIn_pow4.GetScale()
	*/

	for i := 0; i < IterNum; {
		if cty.Level() < 3 {
			if reencryptionMode {
				cty = boostrapper.Bootstrap(cty)
			} else {
				cty = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(cty), params.LogSlots()), maxLevel, params.DefaultScale(), params.LogSlots()))
			}
			btptimes++
		}
		if cty.Level() >= 4 && (i+2) < IterNum {
			/*
				ctGuess_pow2 := evaluator.MulRelinNew(cty, cty)
				err = evaluator.Rescale(ctGuess_pow2, params.DefaultScale(), ctGuess_pow2)
				if err != nil {
					ctOut = nil
					return
				}
				ctGuess_pow4 := evaluator.MulRelinNew(ctGuess_pow2, ctGuess_pow2)
				err = evaluator.Rescale(ctGuess_pow4, params.DefaultScale(), ctGuess_pow4)
				if err != nil {
					ctOut = nil
					return
				}
				ctGuess_pow8 := evaluator.MulRelinNew(ctGuess_pow4, ctGuess_pow4)
				err = evaluator.Rescale(ctGuess_pow8, params.DefaultScale(), ctGuess_pow8)
				if err != nil {
					ctOut = nil
					return
				}

					scale_GuessMod := make([]rlwe.Scale, 4)
					scale_GuessMod[0] = rlwe.NewScale(float64(params.RingQ().Modulus[cty.Level()]))
					scale_GuessMod[1] = rlwe.NewScale(float64(params.RingQ().Modulus[ctGuess_pow2.Level()]))
					scale_GuessMod[2] = rlwe.NewScale(float64(params.RingQ().Modulus[ctGuess_pow4.Level()]))
					scale_GuessMod[3] = rlwe.NewScale(float64(params.RingQ().Modulus[ctGuess_pow8.Level()]))

					scale_Guess := make([]rlwe.Scale, 4)
					scale_Guess[0] = cty.GetScale()
					scale_Guess[1] = ctGuess_pow2.GetScale()
					scale_Guess[2] = ctGuess_pow4.GetScale()
					scale_Guess[3] = ctGuess_pow8.GetScale()
			*/
			// Infact : y -> xy^3 -> x^2y^5....
			// y(0) -> -39/16xy^3 = -39/16y*xy^2(3) -> 27/16x^2y^5 = x^2y^2*y^2*27/16y(4) -> -9/16x^3y^7 = -9/16x
			// we need: xy, y^2, xy^2, x^2y^2,
			ctxy := evaluator.MulRelinNew(cty, ctIn_dup)
			evaluator.Rescale(ctxy, params.DefaultScale(), ctxy)
			//auxio.Quick_check_infos(ctxy, "ctxy")
			ctye2 := evaluator.MulRelinNew(cty, cty)
			evaluator.Rescale(ctye2, params.DefaultScale(), ctye2)
			//auxio.Quick_check_infos(ctye2, "ctye2")
			ctxye2 := evaluator.MulRelinNew(ctxy, cty)
			evaluator.Rescale(ctxye2, params.DefaultScale(), ctxye2)
			//auxio.Quick_check_infos(ctxye2, "ctxye2")
			ctxe2ye2 := evaluator.MulRelinNew(ctxy, ctxy)
			evaluator.Rescale(ctxe2ye2, params.DefaultScale(), ctxe2ye2)
			//auxio.Quick_check_infos(ctxe2ye2, "ctxe2ye2")
			ctye2x1over16 := evaluator.MultByConstNew(ctIn_dup, 1.0/16.0)
			evaluator.Rescale(ctye2x1over16, params.DefaultScale(), ctye2x1over16)
			//auxio.Quick_check_infos(ctye2x1over16, "ctye2x1over16 mult by const")
			evaluator.MulRelin(ctye2x1over16, ctye2, ctye2x1over16)
			evaluator.Rescale(ctye2x1over16, params.DefaultScale(), ctye2x1over16)
			//auxio.Quick_check_infos(ctye2x1over16, "ctye2x1over16 final")
			ctye2x9over16 := evaluator.MultByConstNew(ctye2x1over16, -9)
			// auxio.Quick_check_infos(ctye2x9over16, "ctye2x9over16")
			// (1/16*x^4*y^9) -> y^2(1), x^2(1) -> (1/16*x^2*y^7) -> y^2(1), x^2(1) -> (1/16*y^5) -> y^2(1), y^2(1), y(0), 1/16(0)
			// y^2*x^2(2), y^2*x^2(2), y^2(1), y^2(1), 1/16y(1)
			// y^2*x^2(2), y^2(1), y^2*x*1/16(2), y^2(1), xy(1) (two two combine)
			var ctTerm1 *rlwe.Ciphertext
			ctTerm1 = ctxy.CopyNew() // 12
			evaluator.MulRelin(ctTerm1, ctye2, ctTerm1)
			evaluator.Rescale(ctTerm1, params.DefaultScale(), ctTerm1) // 11
			evaluator.MulRelin(ctTerm1, ctye2x1over16, ctTerm1)
			evaluator.Rescale(ctTerm1, params.DefaultScale(), ctTerm1) // 9
			// auxio.Quick_check_infos(ctTerm1, "ctTerm1 current")
			// auxio.Quick_check_infos(ctye2, "ctye2 original")
			// auxio.Quick_check_infos(ctxe2ye2, "ctxe2ye2 original")
			ctye2tmp := ctye2.CopyNew()
			if ctxe2ye2.Level()-1 <= ctTerm1.Level() {
				evaluator.SetScale(ctye2tmp, params.DefaultScale().Mul(rlwe.NewScale(float64(params.RingQ().Modulus[ctxe2ye2.Level()]))).Div(ctxe2ye2.Scale).Mul(rlwe.NewScale(float64(params.RingQ().Modulus[ctxe2ye2.Level()-1]))).Div(ctTerm1.Scale))
			} else {
				evaluator.SetScale(ctye2tmp, params.DefaultScale().Mul(rlwe.NewScale(float64(params.RingQ().Modulus[ctxe2ye2.Level()]))).Div(ctxe2ye2.Scale).Mul(rlwe.NewScale(float64(params.RingQ().Modulus[ctTerm1.Level()]))).Div(ctTerm1.Scale))
			}
			// auxio.Quick_check_infos(ctye2tmp, "ctye2tmp")
			evaluator.MulRelin(ctye2tmp, ctxe2ye2, ctye2tmp)
			evaluator.Rescale(ctye2tmp, params.DefaultScale(), ctye2tmp)
			// auxio.Quick_check_infos(ctye2tmp, "ctye2tmp multiplied by xe2ye2")
			evaluator.MulRelin(ctye2tmp, ctTerm1, ctTerm1)
			evaluator.Rescale(ctTerm1, params.DefaultScale(), ctTerm1)
			// auxio.Quick_check_infos(ctTerm1, "ctTerm1")
			// auxio.Quick_check_matrix(params, sk, ctTerm1, d, d)

			// (9/16*x^3*y^7) -> y^2(1), x^2(1) -> (9/16*x*y^5) -> y^2(1), x(0) -> (9/16*y^3) -> y^2(1), 9/16y(1)
			// y^2*x*9/16(2), y^2*x(2), y^2(1), xy(1)
			var ctTerm2 *rlwe.Ciphertext
			ctTerm2 = ctxy.CopyNew() // 12
			evaluator.MulRelin(ctTerm2, ctye2, ctTerm2)
			evaluator.Rescale(ctTerm2, params.DefaultScale(), ctTerm2) // 11
			evaluator.MulRelin(ctTerm2, ctxye2, ctTerm2)
			evaluator.Rescale(ctTerm2, params.DefaultScale(), ctTerm2) // 10
			// fmt.Printf("Current level of ctTerm2:%d\n", ctTerm2.Level())
			// fmt.Printf("Original level of ctyex9over16:%d\n", ctye2x9over16.Level())
			if ctye2x9over16.Level() <= ctTerm2.Level() {
				evaluator.SetScale(ctye2x9over16, params.DefaultScale().Mul(rlwe.NewScale(params.RingQ().Modulus[ctye2x9over16.Level()-1])).Div(ctTerm2.Scale))
			} else {
				evaluator.SetScale(ctye2x9over16, params.DefaultScale().Mul(rlwe.NewScale(params.RingQ().Modulus[ctTerm2.Level()])).Div(ctTerm2.Scale))
			}
			// fmt.Printf("Estimated ctye2x9over16's scale: %f=%f+%f-%f\n", math.Log2(params.DefaultScale().Float64())+math.Log2(float64(params.RingQ().Modulus[ctTerm2.Level()]))-math.Log2(ctTerm2.Scale.Float64()), math.Log2(params.DefaultScale().Float64()), math.Log2(float64(params.RingQ().Modulus[ctTerm2.Level()])), math.Log2(ctTerm2.Scale.Float64()))
			// auxio.Quick_check_infos(ctye2x9over16, "ctye2x9over16 real info:")
			evaluator.MulRelin(ctTerm2, ctye2x9over16, ctTerm2)
			// auxio.Quick_check_infos(ctTerm2, "ctTerm2 mult by ctye2x9over16")
			evaluator.Rescale(ctTerm2, params.DefaultScale(), ctTerm2) // 3,2
			// auxio.Quick_check_infos(ctTerm2, "ctTerm2 after rescaling")
			// auxio.Quick_check_matrix(params, sk, ctTerm2, d, d)

			// 27/16x^2y^5 = x^2y^2*y^2*27/16y(3)
			var ctTerm3 *rlwe.Ciphertext
			ctTerm3 = cty.CopyNew()
			evaluator.MultByConst(ctTerm3, 27.0/16.0, ctTerm3)
			evaluator.Rescale(ctTerm3, params.DefaultScale(), ctTerm3)
			evaluator.MulRelin(ctTerm3, ctye2, ctTerm3)
			evaluator.Rescale(ctTerm3, params.DefaultScale(), ctTerm3)
			evaluator.MulRelin(ctTerm3, ctxe2ye2, ctTerm3)
			evaluator.Rescale(ctTerm3, params.DefaultScale(), ctTerm3)
			evaluator.SetScale(ctTerm3, params.DefaultScale())
			// auxio.Quick_check_infos(ctTerm3, "ctTerm3")
			// auxio.Quick_check_matrix(params, sk, ctTerm3, d, d)

			// -39/16xy^3 = -39/16y*xy^2(3)
			var ctTerm4 *rlwe.Ciphertext
			ctTerm4 = cty.CopyNew()
			evaluator.MultByConst(ctTerm4, -39.0/16.0, ctTerm4)
			evaluator.Rescale(ctTerm4, params.DefaultScale(), ctTerm4)
			evaluator.MulRelin(ctTerm4, ctxye2, ctTerm4)
			evaluator.Rescale(ctTerm4, params.DefaultScale(), ctTerm4)
			evaluator.SetScale(ctTerm4, params.DefaultScale())
			//auxio.Quick_check_infos(ctTerm4, "ctTerm4")
			// auxio.Quick_check_matrix(params, sk, ctTerm4, d, d)

			// 9/4 * y
			var ctTerm5 *rlwe.Ciphertext
			ctTerm5 = cty.CopyNew()
			evaluator.MultByConst(ctTerm5, 9.0/4.0, ctTerm5)
			evaluator.Rescale(ctTerm5, params.DefaultScale(), ctTerm5)
			evaluator.SetScale(ctTerm5, params.DefaultScale())
			// auxio.Quick_check_infos(ctTerm5, "ctTerm5")
			// auxio.Quick_check_matrix(params, sk, ctTerm5, d, d)

			// Sum
			evaluator.Add(ctTerm1, ctTerm2, cty)
			evaluator.Add(cty, ctTerm3, cty)
			evaluator.Add(cty, ctTerm4, cty)
			evaluator.Add(cty, ctTerm5, cty)
			// auxio.Quick_check_infos(cty, "cty")
			// auxio.Quick_check_matrix(params, sk, cty, d, d)
			i += 2
		} else if cty.Level() >= 3 {
			ctGuess_pow2 := evaluator.MulRelinNew(cty, cty)
			err = evaluator.Rescale(ctGuess_pow2, params.DefaultScale(), ctGuess_pow2)
			if err != nil {
				ctOut = nil
				return
			}
			ctGuess1frac2 := evaluator.MultByConstNew(cty, -0.5)
			err = evaluator.Rescale(ctGuess1frac2, params.DefaultScale(), ctGuess1frac2)
			//fmt.Printf("ctGuess1frac2 has scale %f, level %d\n", math.Log2(ctGuess1frac2.Scale.Float64()), ctGuess1frac2.Level())
			//auxio.Quick_check_matrix(params, sk, ctGuess1frac2, d, d)
			if err != nil {
				ctOut = nil
				return
			}
			ctGuess1frac2tIn := evaluator.MulRelinNew(ctGuess1frac2, ctIn_dup)
			err = evaluator.Rescale(ctGuess1frac2tIn, params.DefaultScale(), ctGuess1frac2tIn)
			//fmt.Printf("ctGuess1frac2tIn has scale %f, level %d\n", math.Log2(ctGuess1frac2tIn.Scale.Float64()), ctGuess1frac2tIn.Level())
			//auxio.Quick_check_matrix(params, sk, ctGuess1frac2tIn, d, d)
			if err != nil {
				ctOut = nil
				return
			}

			evaluator.SetScale(ctGuess_pow2, params.DefaultScale().Mul(rlwe.NewScale(float64(params.RingQ().Modulus[ctGuess1frac2tIn.Level()]))).Div(ctGuess1frac2tIn.Scale))

			ctGuessPow31frac2tIn := evaluator.MulRelinNew(ctGuess1frac2tIn, ctGuess_pow2)
			err = evaluator.Rescale(ctGuessPow31frac2tIn, params.DefaultScale(), ctGuessPow31frac2tIn)
			//fmt.Printf("ctGuessPow31frac2tIn has scale %f, level %d\n", math.Log2(ctGuessPow31frac2tIn.Scale.Float64()), ctGuessPow31frac2tIn.Level())
			//auxio.Quick_check_matrix(params, sk, ctGuessPow31frac2tIn, d, d)
			if err != nil {
				ctOut = nil
				return
			}
			ctGuess3frac2 := evaluator.MultByConstNew(cty, 1.5)
			err = evaluator.Rescale(ctGuess3frac2, params.DefaultScale(), ctGuess3frac2)
			evaluator.SetScale(ctGuess3frac2, ctGuessPow31frac2tIn.Scale)
			//fmt.Printf("ctGuess3frac2 has scale %f, level %d\n", math.Log2(ctGuess3frac2.Scale.Float64()), ctGuess3frac2.Level())
			//auxio.Quick_check_matrix(params, sk, ctGuess3frac2, d, d)
			if err != nil {
				ctOut = nil
				return
			}
			cty = evaluator.AddNew(ctGuess3frac2, ctGuessPow31frac2tIn)
			i++
		}

	}
	ctOut = cty.CopyNew()
	return
}

/*
func InvSqrtByGS_btpVer_dbg(params ckks.Parameters, rlk *rlwe.RelinearizationKey, pk *rlwe.PublicKey, sk *rlwe.SecretKey, ctIn *rlwe.Ciphertext, ctGuess *rlwe.Ciphertext, IterNum int, maxLevel int, reencryptionMode bool, Boostrapper ...*bootstrapping.Bootstrapper) (ctOut *rlwe.Ciphertext, btptimes int, err error) {
	// d := 128 // only for dbg
	encryptor := ckks.NewEncryptor(params, pk)
	decryptor := ckks.NewDecryptor(params, sk)
	encoder := ckks.NewEncoder(params)
	var boostrapper *bootstrapping.Bootstrapper
	if len(Boostrapper) != 0 {
		boostrapper = Boostrapper[0]
	}
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	var ctX, ctH,ctR *rlwe.Ciphertext
	var ctGuess_dup, ctIn_dup *rlwe.Ciphertext
	if ctGuess.Level() <= 2 {
		if reencryptionMode {
			ctGuess_dup = boostrapper.Bootstrap(ctGuess)
		} else {
			ctGuess_dup = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctGuess), params.LogSlots()), maxLevel, params.DefaultScale(), params.LogSlots()))
		}
		btptimes++
	}
	if ctIn.Level() <= 2 {
		if reencryptionMode {
			ctIn_dup = boostrapper.Bootstrap(ctIn)
		} else {
			ctIn_dup = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctIn), params.LogSlots()), maxLevel, params.DefaultScale(), params.LogSlots()))
		}
		btptimes++
	}
	ctX = evaluator.MulRelinNew(ctGuess_dup,ctIn_dup)
	evaluator.Rescale(ctX,params.DefaultScale(),ctX)
	ctH = evaluator.MultByConstNew(ctGuess_dup,0.5)
	evaluator.Rescale(ctH,params.DefaultScale(),ctH)
	ctR = evaluator.MulRelinNew(ctX,ctH)
	evaluator.Rescale(ctR,params.DefaultScale(),ctR)
	evaluator.MultByConst(ctR,-1,ctR)
	evaluator.AddConst(ctR,0.5,ctR)
	for i :=0;i<IterNum;i++{
		if ctX.Level() <= 2 {
			if reencryptionMode {
				ctX = boostrapper.Bootstrap(ctX)
			} else {
				ctX = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctX), params.LogSlots()), maxLevel, params.DefaultScale(), params.LogSlots()))
			}
			btptimes++
		}
		if ctH.Level() <= 2 {
			if reencryptionMode {
				ctH = boostrapper.Bootstrap(ctH)
			} else {
				ctH = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(ctH), params.LogSlots()), maxLevel, params.DefaultScale(), params.LogSlots()))
			}
			btptimes++
		}

	}

}

*/
