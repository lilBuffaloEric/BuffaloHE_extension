package auxiliary_io

/*
 * This part of package auxiliary_io aims to
 * (1) encapsulate some encoding,decoding and decryption procedures. Specifically, decoding plaintext into a matrix is implented.
 * (2) Create a simple disk I/O structure: Filetracker, in order to handle storing and retrieving Marshalable object. Since most of the structure in Lattigo is marshalable,
 * 	   they can be store and read using our Filetracker.
 */

import (
	"encoding/binary"
	"errors"
	"fmt"
	"math"
	"os"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

type Filetracker struct {
	Fp        string
	Off       int
	Offidx    int
	Len       int
	Varylen   int
	Hierarchy []int
}

type Marshalable interface {
	MarshalBinary() ([]byte, error)
	UnmarshalBinary([]byte) error
	MarshalBinarySize() int
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

func Quick_decode_matrix_full(params ckks.Parameters, pt *rlwe.Plaintext, row_num int, col_num int) {
	encoder := ckks.NewEncoder(params)
	values := encoder.Decode(pt, params.LogSlots())
	realpart := Complex2Float(values)
	Print_matrix_f64_full(realpart, row_num, col_num)
}

func Quick_check_infos(ct *rlwe.Ciphertext, name string) {
	fmt.Printf("%s has scale %f, level %d\n", name, math.Log2(ct.Scale.Float64()), ct.Level())
}

func WriteMany2file(fp string, stuffs []Marshalable, offset int) (n int, err error) {
	var file *os.File
	n = 0
	bytes := 0
	file, err = os.OpenFile(fp, os.O_WRONLY|os.O_CREATE, os.ModeAppend|os.ModePerm)
	if err != nil {
		return 0, err
	}
	defer file.Close()
	buffer := make([][]byte, len(stuffs))
	_, err = file.Seek(int64(offset), 0)
	if err != nil {
		return 0, err
	}
	for i := 0; i < len(stuffs); i++ {
		buffer[i], err = stuffs[i].MarshalBinary()
		if err != nil {
			return 0, err
		}
		bytes, err = file.Write(buffer[i])
		n += bytes
		if err != nil {
			return
		}
	}
	return
}

func ReadMany4file(fp string, size []int, stuffs []Marshalable, offset int) (n int, err error) {
	var file *os.File
	n = 0
	bytes := 0
	file, err = os.Open(fp)
	if err != nil {
		return 0, err
	}
	defer file.Close()
	if len(size) != len(stuffs) {
		return 0, errors.New("given number of size dose not match the number of stuff")
	}
	_, err = file.Seek(int64(offset), 0)
	if err != nil {
		return 0, err
	}
	for i := 0; i < len(stuffs); i++ {
		buffer := make([]byte, size[i])
		bytes, err = file.Read(buffer)
		n += bytes
		if err != nil {
			return
		}
		err = stuffs[i].UnmarshalBinary(buffer)
		if err != nil {
			return
		}
	}
	return
}

func WriteOne2file(fp string, key Marshalable, offset int) (n int, err error) {

	var buffer []byte
	var file *os.File
	n = 0
	file, err = os.OpenFile(fp, os.O_WRONLY|os.O_CREATE, os.ModeAppend|os.ModePerm)
	if err != nil {
		return 0, err
	}
	defer file.Close()

	buffer, err = key.MarshalBinary()
	if err != nil {
		return 0, err
	}

	// 将数据拆分成较小的块，每次写入一小部分数据
	chunkSize := 65 * 1024 // 65KB
	position := 0
	bufferSize := len(buffer)
	written := 0
	for position < bufferSize {
		chunkEnd := position + chunkSize
		if chunkEnd > bufferSize {
			chunkEnd = bufferSize
		}
		chunk := buffer[position:chunkEnd]
		written, err = file.WriteAt(chunk, int64(offset+position))
		if err != nil {
			return
		}
		n += written
		if written != len(chunk) {
			return n, errors.New("failed to write complete chunk")
		}
		position = chunkEnd
	}
	return
}

func ReadOne4file(fp string, key Marshalable, offset int, keylen int) (n int, err error) {
	var buffer []byte
	var file *os.File
	n = 0
	file, err = os.Open(fp)
	if err != nil {
		return 0, err
	}
	defer file.Close()
	buffer = make([]byte, keylen)
	n, err = file.ReadAt(buffer, int64(offset))
	if err != nil {
		return
	}
	err = key.UnmarshalBinary(buffer)
	return
}

func (tracker *Filetracker) StoreFinish() (n int, err error) {
	n, err = WriteOne2file(tracker.Fp, tracker, tracker.Off)
	tracker.Off += n
	return
}

func (tracker *Filetracker) StoreUpdateOne(stuff Marshalable) (n int, err error) {
	if tracker.Varylen == 0 {
		if stuff.MarshalBinarySize() != tracker.Hierarchy[0] {
			err = errors.New("update bytes do not match the tracker's standard")
			return 0, err
		}
	}
	n, err = WriteOne2file(tracker.Fp, stuff, tracker.Off)
	tracker.Off += n
	if tracker.Varylen == 1 {
		tracker.Hierarchy = append(tracker.Hierarchy, n)
		tracker.Len++
	}
	if n != stuff.MarshalBinarySize() || err != nil {
		if err != nil {
			return n, err
		} else {
			return n, errors.New("written bytes do not match the expected MarshalBinarySize")
		}
	}
	return
}

func (tracker *Filetracker) StoreUpdateMany(stuffs []Marshalable) (n int, err error) {
	if tracker.Varylen == 0 {
		for i := 0; i < len(stuffs); i++ {
			if stuffs[i].MarshalBinarySize() != tracker.Hierarchy[0] {
				err = errors.New("update bytes do not match the tracker's standard")
				return 0, err
			}
		}
	}
	n, err = WriteMany2file(tracker.Fp, stuffs, tracker.Off)
	tracker.Off += n
	expectedSum := 0
	for i := 0; i < len(stuffs); i++ {
		expectedSum += stuffs[i].MarshalBinarySize()
	}
	if n != expectedSum || err != nil {
		if err != nil {
			return n, err
		} else {
			return n, errors.New("written bytes do not match the expected MarshalBinarySize")
		}
	}
	if tracker.Varylen == 1 {
		for i := 0; i < len(stuffs); i++ {
			tracker.Hierarchy = append(tracker.Hierarchy, stuffs[i].MarshalBinarySize())
			tracker.Len++
		}
	}
	return
}

func (tracker *Filetracker) ReadInit() (n int, err error) {
	var file *os.File
	n = 0
	file, err = os.Open(tracker.Fp)
	if err != nil {
		return 0, err
	}
	defer file.Close()
	file.Seek(-8, 2)
	buffer := make([]byte, 8)
	file.Read(buffer)
	tracker.Len = int(binary.BigEndian.Uint64(buffer))
	n = 8 + 8 + 8*tracker.Len
	buffer = make([]byte, n)
	file.Seek(int64(-n), 2)
	file.Read(buffer)
	tracker.Off = 0
	tracker.Offidx = 0
	err = tracker.UnmarshalBinary(buffer)
	return
}

func (tracker *Filetracker) ReadUpdateOne(stuff Marshalable) (n int, err error) {
	n = 0
	// version 2
	if tracker.Varylen == 1 {
		n, err = ReadOne4file(tracker.Fp, stuff, tracker.Off, tracker.Hierarchy[tracker.Offidx])
	} else {
		n, err = ReadOne4file(tracker.Fp, stuff, tracker.Off, tracker.Hierarchy[0])
	}
	if err != nil {
		return
	}
	tracker.Off += n
	if tracker.Varylen == 1 {
		tracker.Offidx++
	}
	return

	/* version 1
	n, err = ReadOne4file(tracker.Fp, stuff, tracker.Off, tracker.Hierarchy[0])
	if err != nil {
		return
	}
	tracker.Off += n
	if tracker.Varylen == 1 {
		tracker.Hierarchy = tracker.Hierarchy[1:]
	}
	return
	*/
}

func (tracker *Filetracker) Locate(idxLocation int) {
	tracker.Off = 0
	if tracker.Varylen == 1 {
		for i := 0; i < idxLocation; i++ {
			tracker.Off += tracker.Hierarchy[i]
		}
	} else {
		tracker.Off += idxLocation * tracker.Hierarchy[0]
	}
	tracker.Offidx = idxLocation
	return
}

// untested...
func (tracker *Filetracker) ReadUpdateMany(stuffs []Marshalable) (n int, err error) {
	n = 0
	sizes := make([]int, len(stuffs))
	if tracker.Varylen == 1 {
		for i := 0; i < len(stuffs); i++ {
			sizes[i] = tracker.Hierarchy[i]
		}
	}
	if tracker.Varylen == 0 {
		for i := 0; i < len(stuffs); i++ {
			sizes[i] = tracker.Hierarchy[0]
		}
	}
	n, err = ReadMany4file(tracker.Fp, sizes, stuffs, tracker.Off)
	if err != nil {
		return
	}
	tracker.Off += n
	if tracker.Varylen == 1 {
		tracker.Hierarchy = tracker.Hierarchy[len(stuffs):]
	}
	return
}

func NewTracker(fp string, EntrySize int) (tracker *Filetracker) {
	if EntrySize <= 0 {
		tracker = &Filetracker{Fp: fp, Varylen: 1, Hierarchy: make([]int, 0), Off: 0, Len: 0}
	} else {
		tracker = &Filetracker{Fp: fp, Varylen: 0, Hierarchy: make([]int, 1), Off: 0, Len: 0}
		tracker.Hierarchy[0] = EntrySize
	}
	file, err := os.OpenFile(fp, os.O_WRONLY|os.O_TRUNC|os.O_CREATE, 0666)
	if err != nil {
		fmt.Println("Failed to open file:", err)
		return
	}
	defer file.Close()
	return
}

func NewTracker4File(fp string) (tracker *Filetracker) {
	tracker = &Filetracker{Fp: fp, Varylen: 0, Off: 0, Len: 0}
	tracker.ReadInit()
	return
}

func (tracker *Filetracker) MarshalBinary() (data []byte, err error) {
	data = make([]byte, tracker.MarshalBinarySize())
	ptr := 0
	for i := 0; i < len(tracker.Hierarchy); i++ {
		binary.BigEndian.PutUint64(data[ptr:], uint64(tracker.Hierarchy[i]))
		ptr += 8
	}
	binary.BigEndian.PutUint64(data[ptr:], uint64(tracker.Varylen))
	ptr += 8
	binary.BigEndian.PutUint64(data[ptr:], uint64(tracker.Len))
	return
}

func (tracker *Filetracker) UnmarshalBinary(data []byte) (err error) {
	tracker.Hierarchy = make([]int, tracker.Len)
	for i := 0; i < tracker.Len; i++ {
		tracker.Hierarchy[i] = int(binary.BigEndian.Uint64(data))
		data = data[8:]
	}
	tracker.Varylen = int(binary.BigEndian.Uint64(data))
	data = data[8:]
	tracker.Len = int(binary.BigEndian.Uint64(data))
	return
}

func (tracker *Filetracker) MarshalBinarySize() (datalen int) {
	datalen = 8 + 8 + 8*len(tracker.Hierarchy)
	return
}
