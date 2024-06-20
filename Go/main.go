package main

import (
	"bytes"
	"fmt"
	"log"
	"math"
	"os"
	"strings"

	"gonum.org/v1/gonum/mat"
)

const (
	C = 299792458 // 光速 m/s
	// C = 3e8 // 光速 m/s

	a         = 20 * 1e-3 // 长 mm
	b         = 10 * 1e-3 // 宽 mm
	h         = 1 * 1e-3  // 间隔 mm
	tolerance = 1e-5      // 允许误差
	maxIter   = 1e5
	ezInit    = 1   // Ez内点初始
	kcInit    = 0.5 // kc初始
)

func main() {
	Ez := newEz(a, b, h, ezInit)
	n := Ez.calcTM11(kcInit)
	fmt.Printf("迭代次数: %d\n", n)

	fmt.Print(Ez)
	writeCsv(Ez, "Ez")

	kcTheory, lambdaTheory, freqTheory := theory(a, b)
	kc := Ez.kc()
	lambda, freq := compute(kc)

	fmt.Printf("理论 kc: %.6f\n", kcTheory)
	fmt.Printf("计算 kc: %.6f\n", kc)
	fmt.Printf("理论波长: %.6f m\n", lambdaTheory)
	fmt.Printf("计算波长: %.6f m\n", lambda)
	fmt.Printf("理论频率: %.6f Hz\n", freqTheory)
	fmt.Printf("计算频率: %.6f Hz\n", freq)
}

type Ez struct {
	a, b, h float64
	*mat.Dense
}

func newEz(a, b, h, init float64) *Ez {
	nx := int(a/h) + 1
	ny := int(b/h) + 1
	U := mat.NewDense(nx, ny, nil)

	for i := 1; i < nx-1; i++ {
		for j := 1; j < ny-1; j++ {
			U.Set(i, j, init)
		}
	}

	return &Ez{
		a: a,
		b: b,
		h: h,

		Dense: U,
	}
}

func (Ez *Ez) update(kc float64) {
	r, c := Ez.Dims()
	for i := 1; i < r-1; i++ {
		for j := 1; j < c-1; j++ {
			E2 := Ez.At(i+1, j)
			E3 := Ez.At(i-1, j)
			E4 := Ez.At(i, j+1)
			E5 := Ez.At(i, j-1)
			E1 := (E2 + E3 + E4 + E5) / (4 - kc*kc*Ez.h*Ez.h)
			Ez.Set(i, j, E1)
		}
	}
}

func (Ez *Ez) kc() float64 {
	r, c := Ez.Dims()
	kckcSum := 0.0

	for i := 1; i < r-1; i++ {
		for j := 1; j < c-1; j++ {
			E2 := Ez.At(i+1, j)
			E3 := Ez.At(i-1, j)
			E4 := Ez.At(i, j+1)
			E5 := Ez.At(i, j-1)
			E1 := Ez.At(i, j)
			term := -(E2 + E3 + E4 + E5 - 4*E1) / (E1 * h * h)
			kckcSum += term
		}
	}

	kckcAvg := kckcSum / float64((r-2)*(c-2))

	return math.Sqrt(kckcAvg)
}

func (Ez *Ez) calcTM11(kcInit float64) int {
	r, c := Ez.Dims()

	kc := kcInit
	kcOld := kc + 1
	EzOld := mat.DenseCopyOf(Ez)
	EzDiff := mat.NewDense(r, c, nil)

	n := 0
	for {
		n++
		if n >= maxIter {
			break
		}

		EzOld.Copy(Ez)
		Ez.update(kc)

		if n%10 == 0 {
			kcOld = kc
			kc = Ez.kc()
		}

		EzDiff.Sub(Ez, EzOld)
		EzDiffNorm := mat.Norm(EzDiff, 2)
		kcDiff := math.Abs(kc - kcOld)
		if EzDiffNorm < tolerance || kcDiff < tolerance {
			break
		}
	}

	return n
}

func (Ez *Ez) String() string {
	sb := new(strings.Builder)
	sb.WriteString("Ez:\n")
	r, c := Ez.Dims()
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			sb.WriteString(fmt.Sprintf("%.6f ", Ez.At(i, j)))
		}
		sb.WriteByte('\n')
	}
	return sb.String()
}

func theory(a, b float64) (kc, lambda, freq float64) {
	kc = math.Sqrt(math.Pow(math.Pi/a, 2) + math.Pow(math.Pi/b, 2))
	lambda = 2 * math.Pi / kc
	freq = C / lambda
	return
}

func compute(kc float64) (lambda, freq float64) {
	lambda = 2 * math.Pi / kc
	freq = C / lambda
	return
}

// matrixToCsv 将矩阵转换为逗号分隔符格式
func matrixToCsv(m mat.Matrix) []byte {
	var (
		r, c = m.Dims()
		buf  bytes.Buffer
	)

	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			if j > 0 {
				buf.WriteString(",")
			}
			fmt.Fprintf(&buf, "%v", m.At(i, j))
		}
		buf.WriteString("\n")
	}

	return buf.Bytes()
}

// writeCsv 写入 csv 文件
func writeCsv(m mat.Matrix, name string) {
	filePath := name + ".csv"
	csvFile, err := os.Create(filePath)
	if err != nil {
		log.Fatal(err)
	}
	defer csvFile.Close()
	csvFile.Write(matrixToCsv(m))
}
