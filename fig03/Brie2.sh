## brie2 installed by follow https://brie.readthedocs.io/en/latest/install.html
## SE.filtered.gff3.gz downloaded from https://sourceforge.net/projects/brie-rna/files/annotation/mouse/gencode.vM17/
mkdir brie2
brie-count -a SE.filtered.gff3.gz -S sample.txt -o brie2 -p 38

brie-quant -i brie2/brie_count.h5ad -o brie2/brie_quant_cell.h5ad -c sample.txt --interceptMode gene --LRTindex=All

## file brie2/brie_quant_cell_L5IT.brie_ident.tsv generated for downstream analysis
