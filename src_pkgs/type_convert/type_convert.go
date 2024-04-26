package type_convert

// 将complex128的slice转换为float64的slice
func Complex128ToFloat64(complexSlice []complex128) []float64 {
	floatSlice := make([]float64, len(complexSlice))
	for i, v := range complexSlice {
		floatSlice[i] = real(v)
	}
	return floatSlice
}
