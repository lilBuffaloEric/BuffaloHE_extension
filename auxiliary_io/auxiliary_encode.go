package auxiliary_io

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func Encode_single_int(params ckks.Parameters, value int, level int) (ptOut *rlwe.Plaintext) {
	values := make([]float64, params.Slots())
	encoder := ckks.NewEncoder(params)
	for i := range values {
		values[i] = float64(value)
	}
	scale := rlwe.NewScale(1)
	ptOut = encoder.EncodeNew(values, level, scale, params.LogSlots())
	return
}

func Encode_single_float64(params ckks.Parameters, value float64, level int, scale rlwe.Scale) (ptOut *rlwe.Plaintext) {
	values := make([]float64, params.Slots())
	encoder := ckks.NewEncoder(params)
	for i := range values {
		values[i] = value
	}
	ptOut = encoder.EncodeNew(values, level, scale, params.LogSlots())
	return
}

func DecryptDecode(params ckks.Parameters, sk *rlwe.SecretKey, ct *rlwe.Ciphertext) (values []complex128) {
	decryptor := ckks.NewDecryptor(params, sk)
	encoder := ckks.NewEncoder(params)
	pt := decryptor.DecryptNew(ct)
	values = encoder.Decode(pt, params.LogSlots())
	return
}

func Complex2Float(values []complex128) (realpart []float64) {
	realpart = make([]float64, len(values))
	for i := range values {
		realpart[i] = real(values[i])
	}
	return
}

func Quick_check_matrix(params ckks.Parameters, sk *rlwe.SecretKey, ct *rlwe.Ciphertext, row_num int, col_num int) {
	values := DecryptDecode(params, sk, ct)
	realpart := Complex2Float(values)
	Print_matrix_f64(realpart, row_num, col_num)
}

func Quick_check_matrix_full(params ckks.Parameters, sk *rlwe.SecretKey, ct *rlwe.Ciphertext, row_num int, col_num int) {
	values := DecryptDecode(params, sk, ct)
	realpart := Complex2Float(values)
	Print_matrix_f64_full(realpart, row_num, col_num)
}

func Quick_check_vector(params ckks.Parameters, sk *rlwe.SecretKey, ct *rlwe.Ciphertext) {
	values := DecryptDecode(params, sk, ct)
	realpart := Complex2Float(values)
	Print_vector_f64(realpart, len(realpart))
}

func Quick_decode_vector(params ckks.Parameters, pt *rlwe.Plaintext) {
	encoder := ckks.NewEncoder(params)
	values := encoder.Decode(pt, params.LogSlots())
	realpart := Complex2Float(values)
	Print_vector_f64(realpart, len(realpart))
}
