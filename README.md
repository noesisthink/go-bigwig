---

# gobigwig

[![Go Reference](https://pkg.go.dev/badge/github.com/noesisthink/gobigwig.svg)](https://pkg.go.dev/github.com/noesisthink/gobigwig)

`gobigwig` is a high-performance Go library for parsing BigWig files. It supports cross-platform usage and allows reading genomic signal values, file metadata, chromosome information, and zoom (summary) data efficiently.

It is suitable for genomic data analysis and visualization workflows requiring fast BigWig file processing.

---

## Features

* Cross-platform: Windows / Linux / macOS
* Fast reading of BigWig file headers
* Retrieve chromosome list and signal values
* Support for zoom-level (summary) data
* Rich API for accessing file metadata

---

## Installation

```bash
go get github.com/noesisthink/gobigwig
```

---

## Quick Start

```go
package main

import (
    "fmt"
    "github.com/noesisthink/gobigwig"
)

func main() {
    bw, err := gobigwig.OpenBigWig("example.bw")
    if err != nil {
        panic(err)
    }
    defer gobigwig.CloseBigWig(bw)

    // Print file header metadata
    bw.Getmeta_hdr()

    // Read signal values for chr1 from 1000 to 5000
    signal := bw.ReadBigWigSignal("chr1", 1000, 5000)
    fmt.Println("Signal length:", len(signal))

    // Get chromosome information
    chroms, _ := bw.ShowChromosomes()
    fmt.Println("Chromosomes:", chroms)

    // Get zoom-level summary data
    zoomValues := bw.GetZoomValues("chr1", 0, 100000, 20, true, 10)
    fmt.Println("Zoom values:", zoomValues)
}
```

---

## API Overview

| Feature         | Method                                                                                        | Return Type                 | Description                          |
| --------------- | --------------------------------------------------------------------------------------------- | --------------------------- | ------------------------------------ |
| File operations | `OpenBigWig(fname string)`                                                                    | `(*Bigwig_file_out, error)` | Open a BigWig file                   |
| File operations | `CloseBigWig(fp *Bigwig_file_out)`                                                            | `void`                      | Close the BigWig file                |
| Signal reading  | `ReadBigWigSignal(chrom string, start int, end int)`                                          | `[]float32`                 | Retrieve signal values in a range    |
| Zoom data       | `GetZoomValues(chrom string, start, end, numBins int, useClosest bool, desiredReduction int)` | `[]float32`                 | Retrieve zoom-level summary data     |
| Metadata        | `GetVersion()`                                                                                | `uint16`                    | File version                         |
| Metadata        | `GetNLevels()`                                                                                | `uint16`                    | Number of zoom levels                |
| Metadata        | `GetFieldCount()`                                                                             | `uint16`                    | Total number of fields               |
| Metadata        | `GetDefinedFieldCount()`                                                                      | `uint16`                    | Number of defined fields             |
| Metadata        | `GetBufsize()`                                                                                | `uint32`                    | Buffer size                          |
| Metadata        | `GetExtensionOffset()`                                                                        | `uint64`                    | Extension offset                     |
| Metadata        | `GetNBasesCovered()`                                                                          | `uint64`                    | Number of bases covered              |
| Metadata        | `GetMinVal()`                                                                                 | `float64`                   | Minimum value                        |
| Metadata        | `GetMaxVal()`                                                                                 | `float64`                   | Maximum value                        |
| Metadata        | `GetSumData()`                                                                                | `float64`                   | Sum of data                          |
| Metadata        | `GetSumSquared()`                                                                             | `float64`                   | Sum of squared data                  |
| Chromosome info | `ShowChromosomes()`                                                                           | `map[string]uint64, error`  | Returns chromosome name → length map |
| Chromosome info | `GetChromosomesLen()`                                                                         | `uint64, error`             | Total genome length                  |
| Debugging       | `PrintZoomInfo()`                                                                             | `void`                      | Print zoom-level info                |
| Debugging       | `Getmeta_hdr()`                                                                               | `void`                      | Print file header & chromosome info  |

---

## Zoom Example (Compressed Multi-Chromosome Bar Chart)

```
chr1: 0-1,000,000 bp
+------------------------------------------------+
| ▇▇▆▆▅▅▃▃▂▂▁▁ ▇▇▆▅▃▂▁▁ ▇▆▆▅▃▂▁ ▇▇▇▅▃▂▁ |
+------------------------------------------------+

chr2: 0-500,000 bp
+------------------------------+
| ▆▆▅▃▂▁▁ ▇▇▆▅▃▂▁ ▆▆▅▃▂▁ |
+------------------------------+

chr3: 0-750,000 bp
+----------------------------------+
| ▇▇▇▆▆▅▃▂▁ ▆▆▅▃▂▁ ▇▇▆▆▅▃▂▁ |
+----------------------------------+
```

* Each character represents the average signal intensity of a bin (`▇` high, `▁` low)
* Useful for quick visualization of large regions in genomic data

---

## Example: Retrieve Chromosome Info and Print First 5

```go
bw, _ := gobigwig.OpenBigWig("example.bw")
defer gobigwig.CloseBigWig(bw)

chroms, _ := bw.ShowChromosomes()
count := 0
for chrom, length := range chroms {
    fmt.Printf("%s: %d\n", chrom, length)
    count++
    if count >= 5 {
        break
    }
}
```

---

## Contributing

Issues and Pull Requests are welcome. If you want to improve parsing performance or add new features, please fork the repository and work on a new branch.

---

## License

MIT License © 2025 noesisthink

---

Do you want me to do that next?

