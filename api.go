package gobigwig
import("fmt"
	"math"
	"runtime"
	"sync"
)

type Bigwig_file_out struct {
	bf_fp *bigWigFile_t
	Info FileInfo_bw_out
}

type FileInfo_bw_out struct {
	Version uint16
	NLevels uint16
	FieldCount uint16
	DefinedFieldCount uint16
	Bufsize uint32
	Extensionoffset uint64
	NBasesCovered uint64
	MinVal float64
	MaxVal float64
	SumData float64
	SumSquared float64
}

// openBigWig 打开 BigWig 文件的辅助函数
func OpenBigWig(fname string) (*Bigwig_file_out, error) {
	// 1. 检查是否是 BigWig 文件
	isBw, err := bwisBigWig(fname)
	if err != nil {
		return nil, fmt.Errorf("检查文件格式失败: %w", err)
	}
	if !isBw {
		return nil, fmt.Errorf("不是有效的 BigWig 文件")
	}
	// 2. 打开文件
	url, err := Open(fname)
	if err != nil {
		return nil, fmt.Errorf("打开文件失败: %w", err)
	}
	fp := &bigWigFile_t{
		URL:     url,
		IsWrite: false,
		Type:    0, // 0 = BigWig
	}
	// 3. 读取文件头
	if err := bwHdrRead(fp); err != nil {
		url.Close()
		return nil, fmt.Errorf("读取文件头失败: %w", err)
	}
	// 4. 读取染色体列表
	cl, err := bwReadchromList(fp)
	if err != nil {
		url.Close()
		return nil, fmt.Errorf("读取染色体列表失败: %w", err)
	}
	fp.Cl = cl
	// 5. 读取索引
	idx := bwReadIndex(fp, 0)
	if idx == nil {
		url.Close()
		return nil, fmt.Errorf("读取索引失败")
	}
	fp.Idx = idx
	fbo:=FileInfo_bw_out{
		Version: fp.Hdr.version,
		NLevels: fp.Hdr.nLevels,
		FieldCount: fp.Hdr.fieldCount,
		DefinedFieldCount: fp.Hdr.definedFieldCount,
		Bufsize: fp.Hdr.bufsize,
		Extensionoffset: fp.Hdr.extensionoffset,
		NBasesCovered: fp.Hdr.NBasesCovered,
		MinVal: fp.Hdr.MinVal,
		MaxVal: fp.Hdr.MaxVal,
		SumData: fp.Hdr.SumData,
		SumSquared: fp.Hdr.SumSquared,
	}

	return &Bigwig_file_out{
		bf_fp: fp,
		Info: fbo,
	}, nil
}



// closeBigWig 关闭 BigWig 文件
func CloseBigWig(fp *Bigwig_file_out) {
	if fp.bf_fp != nil && fp.bf_fp.URL != nil {
		fp.bf_fp.URL.Close()
	}
}


func (fp *Bigwig_file_out)ReadBigWigSignal(chrom string, start int, end int) []float32{
	start_uint32:=uint32(start)
	end_uint32:=uint32(end)
	blocksPerIteration := uint32(10) // 每次处理10个块
	iter := bwOverlappingIntervalsIterator(fp.bf_fp, chrom, start_uint32, end_uint32, blocksPerIteration)
	if iter == nil {
		fmt.Println("创建迭代器失败")
		return nil
	}
	totalIntervals := uint32(0)
	iterCount := 0
	output_float32:=[]float32{}
	// 迭代所有数据块
	for iter.Data != nil {
		iterCount++
		intervals := iter.Intervals
		if intervals != nil {
			totalIntervals += intervals.L
		}
		output_float32 = append(output_float32,intervals.Value... )
		iter = bwIteratorNext(iter)
	}
	return output_float32
}


func (fp *Bigwig_file_out)Getmeta_hdr(){
	fmt.Println("\n--- 文件头信息 ---")
	getmeta_hdr(fp.bf_fp)
	fmt.Println("\n--- 染色体信息 ---")
}


func (fp *Bigwig_file_out) GetVersion() uint16 {
	return fp.Info.Version
}

func (fp *Bigwig_file_out) GetNLevels() uint16 {
	return fp.Info.NLevels
}

func (fp *Bigwig_file_out) GetFieldCount() uint16 {
	return fp.Info.FieldCount
}

func (fp *Bigwig_file_out) GetDefinedFieldCount() uint16 {
	return fp.Info.DefinedFieldCount
}

func (fp *Bigwig_file_out) GetBufsize() uint32 {
	return fp.Info.Bufsize
}

func (fp *Bigwig_file_out) GetExtensionOffset() uint64 {
	return fp.Info.Extensionoffset
}

func (fp *Bigwig_file_out) GetNBasesCovered() uint64 {
	return fp.Info.NBasesCovered
}

func (fp *Bigwig_file_out) GetMinVal() float64 {
	return fp.Info.MinVal
}

func (fp *Bigwig_file_out) GetMaxVal() float64 {
	return fp.Info.MaxVal
}

func (fp *Bigwig_file_out) GetSumData() float64 {
	return fp.Info.SumData
}

func (fp *Bigwig_file_out) GetSumSquared() float64 {
	return fp.Info.SumSquared
}

// func (fp *Bigwig_file_out) PrintLevels_detail(){
// 	if fp.bf_fp.Hdr.ZoomHdrs != nil && len(fp.bf_fp.Hdr.ZoomHdrs) > 0 {
// 		zh := fp.bf_fp.Hdr.ZoomHdrs[0] // 取第一个 zoom header
// 		for i := 0; i < len(zh.Level); i++ {
// 			fmt.Printf("Zoom Level %d: level=%d\n",
// 				i, zh.Level[i])
// 		}
// 	}
// }


// PrintZoomInfo 打印zoom level信息（调试用）
func (fp *Bigwig_file_out)PrintZoomInfo() {
	if fp.bf_fp.Hdr == nil || len(fp.bf_fp.Hdr.ZoomHdrs) == 0 {
		fmt.Println("No zoom levels available")
		return
	}

	zhdr := fp.bf_fp.Hdr.ZoomHdrs[0]
	fmt.Println("=== Zoom Levels ===")
	for i, level := range zhdr.Level {
		fmt.Printf("Level %d: reduction=%d, indexOffset=%d, dataOffset=%d\n",
			i, level, zhdr.IndexOffset[i], zhdr.DataOffset[i])
	}
}

// ZoomSelector 定义了选择 zoom level 的函数类型
type ZoomSelector func(zhdr *bwZoomHdr_t, desiredReduction uint32) int

// BWOptions_Zoom 表示 zoom 层级选择和取值的参数
type BWOptions_Zoom struct {
	NumBins        int          // 输出分辨率（输出多少个bin）
	SummaryType    string       // 求值方式，如 "mean" / "max"
	IndexZoomModel ZoomSelector // zoom选择策略函数
}


// GetZoom 获取指定区间的 zoom 数据
// 参数：
//   chrom —— 染色体名称
//   start, end —— 区间范围
//   numBins —— 目标输出bin数
//   useClosest —— 是否使用“最接近”zoom选择模式
//   desiredReduction —— 目标缩放比例
// GetZoomValues 获取指定区间的 zoom 数据（并行版，自动替换 NaN → 0）
func (fp *Bigwig_file_out) GetZoomValues(
	chrom string,
	start int,
	end int,
	numBins int,
	useClosest bool,
	desiredReduction int,
) []float32 {

	opts := BWOptions_Zoom{
		NumBins:     numBins,
		SummaryType: "mean",
	}

	if useClosest {
		opts.IndexZoomModel = bwGetBestZoomClosest
	} else {
		opts.IndexZoomModel = bwSelectBestZoomLevel
	}

	if len(fp.bf_fp.Hdr.ZoomHdrs) == 0 {
		fmt.Println("no zoom headers available")
		return nil
	}

	zhdr := fp.bf_fp.Hdr.ZoomHdrs[0]
	zoomIdx := opts.IndexZoomModel(zhdr, uint32(desiredReduction))
	if zoomIdx < 0 {
		fmt.Printf("no suitable zoom level found for desiredReduction=%d\n", desiredReduction)
		return nil
	}

	values, err := bwGetValuesFromZoom(
		fp.bf_fp, zoomIdx, chrom,
		uint32(start), uint32(end),
		opts.NumBins, opts.SummaryType,
	)
	if err != nil {
		fmt.Printf("failed to read zoom data: %v\n", err)
		return nil
	}

	// 并行替换 NaN 为 0
	n := len(values)
	if n == 0 {
		return values
	}

	numCPU := runtime.NumCPU()
	chunkSize := (n + numCPU - 1) / numCPU
	var wg sync.WaitGroup

	for i := 0; i < numCPU; i++ {
		startIdx := i * chunkSize
		endIdx := startIdx + chunkSize
		if endIdx > n {
			endIdx = n
		}

		if startIdx >= n {
			break
		}
		wg.Add(1)
		go func(s, e int) {
			defer wg.Done()
			for j := s; j < e; j++ {
				if math.IsNaN(float64(values[j])) {
					values[j] = 0
				}
			}
		}(startIdx, endIdx)
	}

	wg.Wait()
	return values
}


func (fp *Bigwig_file_out)GetChromosomesLen() (uint64, error) {
	if fp.bf_fp == nil || fp.bf_fp.Cl == nil {
		return 0, fmt.Errorf("invalid BigWig file or chromosome list is nil")
	}
	// fmt.Println("=== Chromosome Information ===")
	var totalLength uint64 = 0
	for i := 0; i < len(fp.bf_fp.Cl.Chrom); i++ {
		_= fp.bf_fp.Cl.Chrom[i]
		length := uint32(0)
		if i < len(fp.bf_fp.Cl.Len) {
			length = fp.bf_fp.Cl.Len[i]
		}
		totalLength += uint64(length)
		// fmt.Printf("Chrom %d: %-10s Length: %d\n", i, chrom, length)
	}
	// fmt.Printf("Total Chromosomes: %d\n", len(bw.Cl.Chrom))
	// fmt.Printf("Total Genome Length: %d\n", totalLength)

	return totalLength, nil
}

// ShowChromosomes 从BigWig文件中提取染色体信息并返回一个映射（染色体名称到长度）
// 如果发生错误，返回nil和相应的错误信息
func (fp *Bigwig_file_out)ShowChromosomes() (map[string]uint64, error) {
    // 检查输入参数有效性
    if fp.bf_fp == nil {
        return nil, fmt.Errorf("bigWigFile_t指针为nil")
    }
    if fp.bf_fp.Cl == nil {
        return nil, fmt.Errorf("染色体列表(Cl)为nil")
    }
    if fp.bf_fp.Cl.Chrom == nil {
        return nil, fmt.Errorf("染色体名称数组(Chrom)为nil")
    }
    
    // 创建用于存储染色体信息的映射
    // 预设容量为35，兼顾标准染色体和可能的扩展需求
    chromInfo := make(map[string]uint64, 35)
    
    // 遍历所有染色体并收集信息
    for i, chrom := range fp.bf_fp.Cl.Chrom {
        // 检查染色体名称有效性
        if chrom == "" {
            return nil, fmt.Errorf("第%d个染色体名称为空", i)
        }
        
        // 获取染色体长度，处理长度数组可能较短的情况
        var length uint64
        if i < len(fp.bf_fp.Cl.Len) {
            length = uint64(fp.bf_fp.Cl.Len[i])
        } else {
            // 对于没有长度信息的染色体，使用0并记录警告
            length = 0
            // 可以根据需要选择是否输出警告日志
            // log.Printf("警告: 第%d个染色体%s没有长度信息", i, chrom)
        }
        
        // 检查是否有重复的染色体名称
        if _, exists := chromInfo[chrom]; exists {
            return nil, fmt.Errorf("发现重复的染色体名称: %s", chrom)
        }
        
        chromInfo[chrom] = length
    }
    
    // 验证染色体数量
    if len(chromInfo) == 0 {
        return nil, fmt.Errorf("未找到任何染色体信息")
    }
    
    return chromInfo, nil
}

