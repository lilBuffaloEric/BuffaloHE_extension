package auxiliary_io

/*
 * This part of package auxiliary_io mean to encapsulate the encoding,decoding and decryption
 * procedure. Specifically, decoding plaintext into a matrix is implented.
 */

import (
	"os"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

type Filetracker struct {
	Fp        string
	Off       int
	Hierarchy []int
}

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

func WriteCtVec2file(fp string, ctVec []*rlwe.Ciphertext, offset int) (n int, err error) {
	var file *os.File
	n = 0
	bytes := 0
	file, err = os.Open(fp)
	if err != nil {
		return 0, err
	}
	defer file.Close()
	MctVec := make([][]byte, len(ctVec))
	_, err = file.Seek(int64(offset), 0)
	if err != nil {
		return 0, err
	}
	for i := 0; i < len(ctVec); i++ {
		MctVec[i], err = ctVec[i].MarshalBinary()
		if err != nil {
			return
		}
		bytes, err = file.Write(MctVec[i])
		if err != nil {
			return
		}
		n += bytes
	}
	return

}

func WritePtVec2file(fp string, ptVec []*rlwe.Plaintext, offset int) (n int, err error) {
	var file *os.File
	n = 0
	bytes := 0
	file, err = os.Open(fp)
	if err != nil {
		return 0, err
	}
	defer file.Close()
	MptVec := make([][]byte, len(ptVec))
	_, err = file.Seek(int64(offset), 0)
	if err != nil {
		return 0, err
	}
	for i := 0; i < len(ptVec); i++ {
		MptVec[i], err = ptVec[i].MarshalBinary()
		if err != nil {
			return
		}
		bytes, err = file.Write(MptVec[i])
		if err != nil {
			return
		}
		n += bytes
	}
	return
}

func ReadCtVec4file(fp string, ctlen int, ctVec []*rlwe.Ciphertext, offset int) (n int, err error) {
	var file *os.File
	n = 0
	bytes := 0
	file, err = os.Open(fp)
	if err != nil {
		return 0, err
	}
	defer file.Close()
	_, err = file.Seek(int64(offset), 0)
	if err != nil {
		return 0, err
	}
	for i := 0; i < len(ctVec); i++ {
		buffer := make([]byte, ctlen)
		bytes, err = file.Read(buffer)
		n += bytes
		if err != nil {
			return
		}
		err = ctVec[i].UnmarshalBinary(buffer)
		if err != nil {
			return
		}
	}
	return
}

func ReadPtVec4file(fp string, ctlen int, ptVec []*rlwe.Plaintext, offset int) (n int, err error) {
	var file *os.File
	n = 0
	bytes := 0
	file, err = os.Open(fp)
	if err != nil {
		return 0, err
	}
	defer file.Close()
	_, err = file.Seek(int64(offset), 0)
	if err != nil {
		return 0, err
	}
	for i := 0; i < len(ptVec); i++ {
		buffer := make([]byte, ctlen)
		bytes, err = file.Read(buffer)
		n += bytes
		if err != nil {
			return
		}
		err = ptVec[i].UnmarshalBinary(buffer)
		if err != nil {
			return
		}
	}
	return
}

func WriteCt2file(fp string, ctIn *rlwe.Ciphertext, offset int) (n int, err error) {
	var file *os.File
	var MctIn []byte
	n = 0
	file, err = os.Open(fp)
	if err != nil {
		return 0, err
	}
	defer file.Close()
	MctIn, err = ctIn.MarshalBinary()
	if err != nil {
		return
	}
	n, err = file.WriteAt(MctIn, int64(offset))
	return
}

func ReadCt4file(fp string, ctIn *rlwe.Ciphertext, offset int, ctlen int) (n int, err error) {
	var file *os.File
	var MctIn []byte
	n = 0
	file, err = os.Open(fp)
	if err != nil {
		return 0, err
	}
	defer file.Close()
	MctIn = make([]byte, ctlen)
	n, err = file.ReadAt(MctIn, int64(offset))
	if err != nil {
		return
	}
	err = ctIn.UnmarshalBinary(MctIn)
	return
}
