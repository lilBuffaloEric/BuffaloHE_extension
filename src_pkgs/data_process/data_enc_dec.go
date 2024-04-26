package data_process

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

// 编码和加密二维float64的slice, 若encryptor为nil则返回明文, 否则返回密文
func EncMatrixF64(matrix [][]float64, level int, scale rlwe.Scale, logSlots int, encoder ckks.Encoder, encryptor rlwe.Encryptor) (ciphertext *rlwe.Ciphertext, plaintext *rlwe.Plaintext) {
	plainVals := make([]float64, 1<<logSlots)
	for i, v := range matrix {
		copy(plainVals[i*len(matrix[0]):], v)
	}
	plaintext = encoder.EncodeNew(plainVals, level, scale, logSlots)
	if encryptor == nil {
		return nil, plaintext
	}
	return encryptor.EncryptNew(plaintext), nil
}

// 编码和加密一维float64的slice, 并填充count次, 若encryptor为nil则返回明文, 否则返回密文
func EncRowVectorF64(vector []float64, count int, level int, scale rlwe.Scale, logSlots int, encoder ckks.Encoder, encryptor rlwe.Encryptor) (ciphertext *rlwe.Ciphertext, plaintext *rlwe.Plaintext) {
	plainVals := make([]float64, 1<<logSlots)
	for i := 0; i < count; i++ {
		copy(plainVals[i*len(vector):], vector)
	}
	plaintext = encoder.EncodeNew(plainVals, level, scale, logSlots)
	if encryptor == nil {
		return nil, plaintext
	}
	return encryptor.EncryptNew(plaintext), nil
}

// 将密文解密为complex128的slice
func DecCiphertext(ciphertext *rlwe.Ciphertext, decryptor rlwe.Decryptor, encoder ckks.Encoder, logSlots int) []complex128 {
	plainVals := decryptor.DecryptNew(ciphertext)
	return encoder.Decode(plainVals, logSlots)
}
