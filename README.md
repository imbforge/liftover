# liftover
A tool to convert genome coordinates and genome annotation files between assemblies.
Underneath, it uses the fantastic [`rtracklayer`](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html) package from Bioconductor, and the [chain files](https://hgdownload.soe.ucsc.edu/downloads.html#liftover) supplied by UCSC.

## Current limitations
* it assumes a UCSC genome. Even if you install new chains based on Ensembl, there's hardcoded a `seqlevelsStyle(x) <- "UCSC"` somewhere in the code.
* poor handling of overlapping ranges. Formats like BigWig don't allow overlapping ranges. We're simply skipping them, which is really very poor handling. Be careful with that!
