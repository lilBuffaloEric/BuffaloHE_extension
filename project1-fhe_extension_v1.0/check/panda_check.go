package check

import (
	"fmt"
	"time"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	auxio "project1-fhe_extension_v1.0/auxiliary_io"
)

func linear_approx(params ckks.Parameters, rlk *rlwe.RelinearizationKey, ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	a := -0.00019703
	b := 0.14777278
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	ctOut = evaluator.MultByConstNew(ctIn, a)
	evaluator.Rescale(ctOut, params.DefaultScale(), ctOut)
	evaluator.AddConst(ctOut, b, ctOut)
	return
}

func add_row_elements(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext, offset int) (ctOut *rlwe.Ciphertext) {
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	ctOut = ctIn.CopyNew()
	if offset == 1 {
		for i := 0; i < params.LogSlots(); i++ {
			temp := evaluator.RotateNew(ctOut, 1)
			evaluator.Add(ctOut, temp, ctOut)
		}
	} else {
		for i := 1; i < params.Slots()/256; i = i << 1 {
			temp := evaluator.RotateNew(ctOut, 1)
			evaluator.Add(ctOut, temp, ctOut)
		}

	}
	return

}

func expand_vector(params ckks.Parameters, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	var temp *rlwe.Ciphertext
	ctOut = ctIn.CopyNew()
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	for i := 1; i < 256; i = i << 1 {
		temp = evaluator.RotateNew(ctOut, 1)
		evaluator.Add(ctOut, temp, ctOut)
	}
	return
}

func partial_add(params ckks.Parameters, pk *rlwe.PublicKey, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	encoder := ckks.NewEncoder(params)
	encryptor := ckks.NewEncryptor(params, pk)
	ctOut = ctIn.CopyNew()
	temp_vec := make([]float64, params.Slots())
	temp := encryptor.EncryptNew(encoder.EncodeNew(temp_vec, 11, params.DefaultScale(), params.LogSlots()))
	for i := 1; i <= 256; i = i << 1 {

		temp2 := evaluator.RotateNew(temp, 1)
		s1 := evaluator.MulRelinNew(ctOut, temp2)
		evaluator.Rescale(s1, params.DefaultScale(), s1)
		evaluator.Rotate(s1, 1, s1)

		s2 := evaluator.MulRelinNew(ctOut, temp2)
		evaluator.Rescale(s1, params.DefaultScale(), s2)
		evaluator.Rotate(s2, 1, s2)

		evaluator.Add(s1, s2, ctOut)
		if i != 256 {
			evaluator.Rotate(temp2, 1, temp2)
			evaluator.MulRelin(temp, temp2, temp)
			evaluator.Rescale(temp, params.DefaultScale(), temp)
		}
	}
	ctTemp := ctOut.CopyNew()
	ctOut = expand_vector(params, rlk, galk, ctTemp)
	return
}

func sum_of_squares(params ckks.Parameters, sk *rlwe.SecretKey, pk *rlwe.PublicKey, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	norm := evaluator.MulRelinNew(ctIn, ctIn)
	evaluator.Rescale(norm, params.DefaultScale(), norm)
	norm = partial_add(params, pk, rlk, galk, norm.CopyNew())
	decryptor := ckks.NewDecryptor(params, sk)
	encryptor := ckks.NewEncryptor(params, pk)
	encoder := ckks.NewEncoder(params)
	norm = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(norm.CopyNew()), params.LogSlots()), 11, params.DefaultScale(), params.LogSlots()))
	ctOut = norm.CopyNew()
	return
}

func inv_norm_approx(params ckks.Parameters, sk *rlwe.SecretKey, pk *rlwe.PublicKey, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	norm := ctIn.CopyNew()
	x := norm.CopyNew()
	guess := linear_approx(params, rlk, norm)
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	div_vec := make([]float64, params.Slots())
	div1 := encoder.EncodeNew(div_vec, 11, params.DefaultScale(), params.LogSlots())
	div2 := encoder.EncodeNew(div_vec, 11, params.DefaultScale(), params.LogSlots())
	var sum_x *rlwe.Ciphertext
	r := 2
	for i := 0; i < r; i++ {
		square_x := evaluator.MulRelinNew(guess, guess)
		evaluator.Rescale(square_x, params.DefaultScale(), square_x)

		out_x := evaluator.MulRelinNew(guess, x)
		evaluator.Rescale(out_x, params.DefaultScale(), out_x)

		evaluator.MulRelin(square_x, div1, square_x)
		evaluator.Rescale(square_x, params.DefaultScale(), square_x)

		evaluator.MulRelin(square_x, out_x, square_x)
		evaluator.Rescale(square_x, params.DefaultScale(), square_x)

		sum_x = evaluator.MulRelinNew(guess, div2)
		evaluator.Rescale(sum_x, params.DefaultScale(), sum_x)

		evaluator.Add(sum_x, square_x, guess)
	}
	ctOut = guess.CopyNew()
	return

}

func goldschmidt(params ckks.Parameters, sk *rlwe.SecretKey, pk *rlwe.PublicKey, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, ctIn *rlwe.Ciphertext) (ctX *rlwe.Ciphertext, ctH *rlwe.Ciphertext) {
	encryptor := ckks.NewEncryptor(params, pk)
	decryptor := ckks.NewDecryptor(params, sk)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	encoder := ckks.NewEncoder(params)
	// temp := 0
	norm := sum_of_squares(params, sk, pk, rlk, galk, ctIn)
	y := inv_norm_approx(params, sk, pk, rlk, galk, norm)
	y = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(y.CopyNew()), params.LogSlots()), 11, params.DefaultScale(), params.LogSlots()))

	x := evaluator.MulRelinNew(norm, y)
	evaluator.Rescale(x, params.DefaultScale(), x)

	h := evaluator.MultByConstNew(y, 0.5)
	evaluator.Rescale(h, params.DefaultScale(), h)

	var r *rlwe.Ciphertext

	k := 4
	for i := 0; i < k; i++ {
		temp_r := evaluator.MulRelinNew(x, h)
		evaluator.Rescale(temp_r, params.DefaultScale(), temp_r)
		evaluator.MultByConst(temp_r, -1, temp_r)

		r = evaluator.AddConstNew(temp_r, 0.5)

		temp_x := evaluator.MulRelinNew(x, r)
		evaluator.Rescale(temp_x, params.DefaultScale(), temp_x)

		evaluator.Add(x, temp_x, x)

		temp_h := evaluator.MulRelinNew(h, r)
		evaluator.Rescale(temp_h, params.DefaultScale(), temp_h)

		evaluator.Add(h, temp_h, h)

		if x.Level() <= 2 {
			x = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(x.CopyNew()), params.LogSlots()), 11, params.DefaultScale(), params.LogSlots()))

		}
		if h.Level() <= 2 {

			h = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(h.CopyNew()), params.LogSlots()), 11, params.DefaultScale(), params.LogSlots()))
		}

	}
	evaluator.MultByConst(h, 2, h)
	return x, h
}

func vect_mat_product(params ckks.Parameters, sk *rlwe.SecretKey, pk *rlwe.PublicKey, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, Xj *rlwe.Ciphertext, r *rlwe.Ciphertext) (ctOut *rlwe.Ciphertext) {
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})
	x := evaluator.MulRelinNew(Xj, r)
	evaluator.Rescale(x, params.DefaultScale(), x)
	x = partial_add(params, pk, rlk, galk, x.CopyNew())
	evaluator.MulRelin(x, Xj, x)
	evaluator.Rescale(x, params.DefaultScale(), x)
	ctOut = x.CopyNew()
	return
}

func power_iteration(params ckks.Parameters, sk *rlwe.SecretKey, pk *rlwe.PublicKey, rlk *rlwe.RelinearizationKey, galk *rlwe.RotationKeySet, X []*rlwe.Ciphertext) (eigvec []*rlwe.Ciphertext, eigval []*rlwe.Ciphertext) {
	encryptor := ckks.NewEncryptor(params, pk)
	decryptor := ckks.NewDecryptor(params, sk)
	encoder := ckks.NewEncoder(params)
	evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk})
	iterations := 4
	r_vec := make([]float64, params.Slots())
	r := encryptor.EncryptNew(encoder.EncodeNew(r_vec, 11, params.DefaultScale(), params.LogSlots()))
	eigval = make([]*rlwe.Ciphertext, 0)
	eigvec = make([]*rlwe.Ciphertext, 0)
	for d := 0; d < 4; d++ {
		fmt.Printf("Computing %d-th eigenvec\n", d)
		var eig_val *rlwe.Ciphertext
		for k := 0; k < iterations; k++ {
			fmt.Printf("Computing %d-th PowerMet iteration\n", k)
			var S1 *rlwe.Ciphertext
			S1 = encryptor.EncryptZeroNew(11)
			var S2 *rlwe.Ciphertext
			S2 = encryptor.EncryptZeroNew(11)
			for j := range X {
				x := vect_mat_product(params, sk, pk, rlk, galk, X[j], r)
				evaluator.Add(S1, x, S1)
			}
			S1 = add_row_elements(params, rlk, galk, S1.CopyNew(), 256)
			S1 = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(S1.CopyNew()), params.LogSlots()), 11, params.DefaultScale(), params.LogSlots()))
			for j := range eigvec {
				x := vect_mat_product(params, sk, pk, rlk, galk, eigvec[j], r)
				evaluator.MulRelin(x, eigval[j], x)
				evaluator.Rescale(x, params.DefaultScale(), x)
				evaluator.Add(S2, x, S2)
			}
			S2 = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(S2.CopyNew()), params.LogSlots()), 11, params.DefaultScale(), params.LogSlots()))
			evaluator.Sub(S1, S2, S1)
			eig_val1, eig_inv := goldschmidt(params, sk, pk, rlk, galk, S1)
			evaluator.MulRelin(S1, eig_inv, S1)
			r = encryptor.EncryptNew(encoder.EncodeNew(encoder.Decode(decryptor.DecryptNew(S1.CopyNew()), params.LogSlots()), 11, params.DefaultScale(), params.LogSlots()))
			eig_val = eig_val1.CopyNew()
		}
		eigval = append(eigval, eig_val)
		eigvec = append(eigvec, r)
	}
	return

}

func PandaPCA_check(loadkeysFromDisk bool) {
	var err error
	var KeysPath = "Keys"
	// initialize encryption scheme params.
	var params ckks.Parameters
	var sk *rlwe.SecretKey
	var pk *rlwe.PublicKey
	var rlk *rlwe.RelinearizationKey
	var kgen rlwe.KeyGenerator

	if loadkeysFromDisk {
		params, sk, pk, rlk, err = LoadKeys_check()
		if err != nil {
			panic(err)
		}
		kgen = ckks.NewKeyGenerator(params)
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

	// inputLevel := 11

	// Encryptor
	encryptor := ckks.NewEncryptor(params, pk)

	// Decryptor
	// decryptor := ckks.NewDecryptor(params, sk)

	// encoder
	encoder := ckks.NewEncoder(params)

	// galk
	steps := make([]int, 1)
	steps[0] = 1
	galk := kgen.GenRotationKeysForRotations(steps, false, sk)

	// evaluator
	// evaluator := ckks.NewEvaluator(params, rlwe.EvaluationKey{Rlk: rlk, Rtks: galk})

	X := make([]*rlwe.Ciphertext, 4)
	vec := make([]float64, params.Slots())
	for i := range X {
		X[i] = encryptor.EncryptNew(encoder.EncodeNew(vec, 11, params.DefaultScale(), params.LogSlots()))
	}
	now := time.Now()
	power_iteration(params, sk, pk, rlk, galk, X)
	fmt.Printf("Panda done in %s\n", time.Since(now))
}
