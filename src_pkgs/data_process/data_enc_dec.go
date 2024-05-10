package data_process

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// 编码和加密二维float64的slice, 若publicKey为nil则返回明文, 否则返回密文
func EncMatrixF64(matrix [][]float64, level int, scale rlwe.Scale, logSlots int, params ckks.Parameters, publicKey *rlwe.PublicKey) (ciphertext *rlwe.Ciphertext, plaintext *rlwe.Plaintext) {
	encoder := ckks.NewEncoder(params)
	plainVals := make([]float64, 1<<logSlots)
	for i, v := range matrix {
		copy(plainVals[i*len(matrix[0]):], v)
	}
	plaintext = encoder.EncodeNew(plainVals, level, scale, logSlots)
	if publicKey == nil {
		return nil, plaintext
	} else {
		encryptor := ckks.NewEncryptor(params, publicKey)
		return encryptor.EncryptNew(plaintext), nil
	}
}

// 编码和加密一维float64的slice, 并填充count次, 若publicKey为nil则返回明文, 否则返回密文
func EncRowVectorF64(vector []float64, count int, level int, scale rlwe.Scale, logSlots int, params ckks.Parameters, publicKey *rlwe.PublicKey) (ciphertext *rlwe.Ciphertext, plaintext *rlwe.Plaintext) {
	encoder := ckks.NewEncoder(params)

	plainVals := make([]float64, 1<<logSlots)
	for i := 0; i < count; i++ {
		copy(plainVals[i*len(vector):], vector)
	}
	plaintext = encoder.EncodeNew(plainVals, level, scale, logSlots)
	if publicKey == nil {
		return nil, plaintext
	} else {
		encryptor := ckks.NewEncryptor(params, publicKey)
		return encryptor.EncryptNew(plaintext), nil
	}
}

// 将密文解密为complex128的slice
func DecCiphertext(ciphertext *rlwe.Ciphertext, logSlots int, params ckks.Parameters, sk *rlwe.SecretKey) []complex128 {
	encoder := ckks.NewEncoder(params)
	decryptor := ckks.NewDecryptor(params, sk)
	plainVals := decryptor.DecryptNew(ciphertext)
	return encoder.Decode(plainVals, logSlots)
}
