
## read GEO data
library( GEOquery )
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072*6 )
meta =getGEO(filename = "GSE244173_series_matrix.txt.gz", AnnotGPL = F, getGPL = F)
meta = pData(meta)
head( meta )

### read in expr matrix
xxx = Read10X( "./GSE244173_RAW/GSM7808275_CG-3Dose-F-3" )


### 
class( xxx ) 
# xxx = as( xxx, "dgCMatrix") # adjust the class to dgCMatrix if necessary

object.size( xxx ) # check mem
table( xxx@x < 0 ) # 
table( xxx@x == 0 ) #


hist( xxx@x, breaks = 100 ) # 
hist( log2( xxx@x), breaks = 100 ) # 

genexx = dimnames( xxx )[[1]]
head( genexx )
cellxx = dimnames( xxx )[[2]]
head( cellxx )

metadata = data.frame( row.names = cellxx ) %>% 
  mutate( . , sample = "xxx" ) %>%
  mutate( individual = "C199E" ) %>% 
  mutate( vaccine = "pre" )


### Seurat object
yyy = CreateSeuratObject( counts = xxx, # expr matrix
                          meta.data = metadata, # meta
                          project = "xxx_1", # 
                          min.cells = 5, #
                          min.features = 100 #
                          )


### merge Seurat objects
xxx2 =  Read10X( "./GSE244173_RAW/GSM7808276_CG-3Dose-Af-3" ) 
dim( xxx2 )
cellxx2 = dimnames( xxx2 )[[2]]
metadata2 = data.frame( row.names = cellxx2 ) %>% 
  mutate( sanmple = "xxx2" ) %>%
  mutate( individual = "C199E" ) %>% 
  mutate( vaccine = "post" )
#
yyy2 = CreateSeuratObject( counts = xxx2, 
                          meta.data = metadata2,
                          project = "xxx_2",
                          min.cells = 100,
                          min.features = 100
)


# merging
sc_int = merge( yyy, yyy2 )
# for more samples, to merge, sc_int = merge( yyy, list( yyy2, yyy3 ) )


