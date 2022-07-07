library(data.table)
data.table::setDTthreads(threads = 2)
library(stringi)
library(ggplot2)
library(bettermc)
library(magrittr)
library(qqman)

# Patients and folders
# ----
selectedFolders = c(
  "folderName1",
  "folderName2"
)
selectedPatients = c(
  "sampleName1", "sampleName2"
)
# ----

# Annotation is generated and exported
if (T) {
  library(annotatr)
  library(AnnotationHub)
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", path = .libPaths()[1]) 
  # This package is necessary, but not automatically installed
  hg38Annots = builtin_annotations()[
    stri_detect_fixed(builtin_annotations(), "hg38") # all hg38 are built, although not used
  ]
  hg38Annots = hg38Annots[hg38Annots != "hg38_cpg_inter" &
      hg38Annots != "hg38_enhancers_fantom"] # We remove the "regions that are not CpG islands" annotation
  # and enhancers_fantom
  hg38Anno = build_annotations(genome = "hg38", annotations = hg38Annots)
  hg38Anno = hg38Anno[GenomicRanges::mcols(hg38Anno)$type != "hg38_cpg_inter"] # The "regions that are not CpG islands" has to be removed again
  saveRDS(hg38Anno, "/path/to/annotationsMulti.rds")
}

# Do the annotation
if (T) {
hg38Anno = readRDS("/path/to/annotationsMulti.rds")
rmResults = annotatr::read_regions( # rmResults: results from RepeatMasker, as a data.table
  con = "/path/to/ins/data/all.merged.ins.85.min3.rm.bed",
  format = "bed",
  genome = "hg38",
  rename_score = "meaningless", # This column is necessary in BED format, but it has no meaning here
  extraCols = c(
    seqId = "character",
    sw_score = "numeric",
    repeat.class = "character",
    repeat.subclass = "character",
    repeat.start = "numeric",
    repeat.end = "numeric",
    repeat.left = "numeric",
    repeat.strains = "character",
    repeat.divergence_percentage = "numeric",
    repeat.deletion_percentage = "numeric",
    repeat.insertion_percentage = "numeric",
    vcf_alt = "character",
    vcf_info = "character",
    setNames(rep("character", times = length(sampleNames)), paste0("id", sampleNames))
  )
) %>% data.frame() %>% data.table()
} # Reading RepeatMasker's results

if (T) {
  # Generating the usable datasets: min3 and trusty, lax and strict critera, respectively
  genoFields = c("GT", "PSV", "LN", "DR", "ST", "QV", "TY", "ID", "RAL", "AAL", "CO", "SC")
  genoFieldsOfInt = c("PSV", "DR") # Used for counting occurrence and genotyping
  pos = which(genoFields %in% genoFieldsOfInt) # THe position of "PSV" and "DR" in the genotype columns
  
  for (col in colnames(rmResults[1, .SD, .SDcols = patterns( 
    paste0("^(" , paste0(paste0("id", selectedPatients %>% stri_replace_all_fixed(., pattern = "-", replacement = ".")), collapse = "|"), ")$")
    )]) ) {
    # create columns for PSV and DR for each patient
    rmResults[, paste0(col, "_", genoFieldsOfInt) := tstrsplit(.SD[[col]], ":", fixed = T, keep = pos)]
    # keep second value of DR, which is the one supporting the INS
    # Variant callers dont consistently report the support for the reference allele
    rmResults[, paste0(col, "_DR") := tstrsplit(.SD[[paste0(col, "_DR")]], ",", fixed = T, keep = 2)]
    # Fix NAs
    rmResults[, 
        paste0(col, "_DR") := 
          fifelse(.SD[[paste0(col, "_DR")]] == ".", 0, as.numeric(.SD[[paste0(col, "_DR")]]) )
        ] # This will warn about NA by coercion, which are not included in the data.table
  }
  
  # Take the DR columns, create a bool vector, make it numeric, assign as column
  rmResults[, SUPP_VEC_min3 := apply(.SD, 2, `>=`, y = 3) %>% apply(., 2, as.integer) %>% apply(., 1, paste0, collapse = ""), .SDcols = grep(pattern = "_DR$", colnames(rmResults))]
  # For every INS, we have created a bool column containing the min3 criterion
  
  rmResults[, SUPP_min3 := stri_count_fixed(SUPP_VEC_min3, "1")] # Count SUPPORT with the min3 criteria
  
  # Initialize the column with only the two callers criteria
  # Since a data.table is a list of columns (like a data.frame)
  # apply with margin = 2 iterates rows instead of columns, and the opposite for margin = 1
  # because it is "transposed"
  rmResults[
    , SUPP_VEC_trusty := apply(.SD, 2, `==`, y = "11") %>% 
      apply(., 2, as.integer) %>% 
      apply(., 1, paste0, collapse = ""), 
    .SDcols = patterns("_PSV$")]
  # The value for the column is provisional, and missing the RE >= 3 part
  
  # This should be in R, but alas
  # ("00", "11") => "00"
  # ("01", "11") => "01"
  # ("11", "11") => "11"
  bitwiseAndStr = function (a, b) { 
    lapply (1:length(a), function (i) {
      a2 = unlist( strsplit(a[[i]], split = "", fixed = T) )
      b2 = unlist( strsplit(b[[i]], split = "", fixed = T) )
      
      paste0(ifelse(a2 == b2 & a2 == "1", "1", "0"), collapse = "")
    }
    )
  }
  
  # Now we do both callers AND min3 -> strict criteria
  correctTrusty = unlist( bitwiseAndStr(rmResults$SUPP_VEC_min3, rmResults$SUPP_VEC_trusty) )
  rmResults[, SUPP_VEC_trusty := correctTrusty]
  rm(correctTrusty)
  rmResults[, SUPP_trusty := stri_count_fixed(SUPP_VEC_trusty, "1")]
  
  # We will also need the INS' length
  rmResults[, SVLEN := as.numeric(gsub(".*SVLEN=([^;]+);.*", "\\1", vcf_info))]
  
  # A vector column of supports from callers, for genotyping
  supports = rmResults[, .SD, .SDcols = grep("DR$", colnames(rmResults))]
  supportColumn = supports %>% as.matrix() %>% t() %>% split(rep(1:ncol(.), each = nrow(.)))
  rmResults[, support := supportColumn]
  rm(supports, supportColumn) ; gc()
  
  # Create vector bool columns
  suppMin3 = rmResults$SUPP_VEC_min3
  rmResults[, SUPP_mask_min3 := suppMin3 %>% 
              lapply(strsplit, split = "") %>% 
              lapply(function(x)as.logical(as.numeric(unlist(x))))]
  suppTrusty = rmResults$SUPP_VEC_trusty
  rmResults[, SUPP_mask_trusty := suppTrusty %>% 
              lapply(strsplit, split = "") %>% 
              lapply(function(x)as.logical(as.numeric(unlist(x))))]
  
  rmResults[, # split level2-level3 and replace NA for level 3 with "-"
    c("repeat.subclass", "repeat.family") := stri_match_all_regex(repeat.subclass, "([^-]+)(?:-([^-]+))?$") %>% rapply(., f = `[`, ... = 2:3, how = "r") %>% 
        lapply(function (x) {x[is.na(x)] = "-"; x}) %>% data.table::transpose() ]
  
  # Remove unnecessary columns
  rmResults[, grep("^id", colnames(rmResults)) := NULL]
  rmResults[, vcf_info := NULL]
  
  # Now lets annotate
  annotatedFilt = data.table(data.frame(
      annotatr::annotate_regions(regions = GenomicRanges::makeGRangesFromDataFrame(
        rmResults, 
        seqnames.field = "seqnames", 
        keep.extra.columns = T),
                       annotations = hg38Anno,
                       ignore.strand = F,
                       quiet = T
      )
    ))
  annotatedFilt[, annot.type := stri_replace_first_regex(annot.type, "hg38_(.*)", "$1")] # "hg38_genes_introns" => "genes_introns"
  # Lets separate "genes_introns" into 2 columns
  text = stri_match_all_regex(annotatedFilt$annot.type, "([^-]*)_(.*)") %>% 
    lapply(., function(x) {x[2:3]})
  annotatedFilt[, `:=`(annot.class = text %>% lapply(`[`, 1) %>% unlist(),
                 annot.subclass = text %>% lapply(`[`, 2) %>% unlist()
                 )]
  annotatedFilt[, annot.type := NULL] # Redundant
  annotatedFilt[, repeat.class := as.factor(repeat.class)]
  rm(text); gc()
  
  annotatedFilt[, c("annot.strand", "annot.seqnames", "annot.width", "annot.id", "annot.tx_id") := NULL]
  # We only save the ones that pass at least the min3 criteria for at least one occurrence
  saveRDS(annotatedFilt[SUPP_min3 > 0], "/path/to/annotatedFiltered.rds")
  # The ones that dont could get save like this:
  # saveRDS(annotatedFilt[SUPP_min3 == 0], "/path/to/annotatedFilteredLowRe.rds")
}

# 
# Genotyping
# 
# ----
chrs = paste0("chr", c(1:22, "X", "Y")) # Human chromosomes used to filter out Unknown and Random fragments
# chrM and Epstein-Barr are removed, too, but they had no presence in our study
class1 = c("Retroposon", "SINE", "LINE", "LTR") # Families of class I elements
class2 = c("DNA") # Class II elements

repeats = readRDS("/path/to/annotatedFiltered.rds")
repeats[, repeat.percentage := (repeat.end - repeat.start + 1) / SVLEN] # Percentage of INS covered by annotation
# We filter out the ones that are not >=85% repetitive elements (includes TEs and non TEs),
# the ones that are not TEs
# and the ones not in chr1-22, X or Y
repeats = repeats[repeat.percentage >= 0.85 & repeat.class %in% c(class1, class2) & seqnames %in% chrs] 
ins = repeats[, .SD[1], by = seqId] # Only one annotation per INS, a "unique" version of the above
ins[, grep("^annot\\.", colnames(ins)) := NULL] # Remove the annotation data (incomplete in this table)
allIns = readRDS("/path/to/annotatedFiltered.rds")[ # An unfiltered version, with only one row per INS (so no annotation)
  , repeat.percentage := (repeat.end - repeat.start + 1) / SVLEN
][, .SD[1], by = seqId]

if (T) { # Calculate coverage on affected sites, on all patients, with Mosdepth, for genotyping 
  mosdepthBed = ins[seqnames %in% chrs, c("seqnames", "start", "end", "seqId")]
  mosdepthBed[, `:=`(
    start = pmax(start - 61, 0), # The extra 1 adjusts to BED coordinates
    end = end + 60
  )]
  # mosdepthBed %>% nrow()
  fwrite(mosdepthBed, file = "/path/to/insertions.bed", sep = "\t", col.names = F, scipen=50)
  
  system("bash slurm/mosdepth.slurm config/config.sh INS")
  
  seqids = mosdepthBed[order(seqId)]$seqId
  coveragesForGenotyping = mclapply(selectedPatients, function (x) {
    fread(paste0("/path/to/", x, ".minimap.ins.regions.bed.gz"),
                  col.names = c("chr", "start", "end", "seqId", "coverage"))[order(seqId)]$coverage
  }, mc.cores = 1)
  coverageForGenotyping = data.table(
    seqId = seqids, 
    coverage = lapply(asplit(simplify2array(coveragesForGenotyping), 1), as.numeric)
    )
}
colnamesoi = c("seqnames", "start", "end", "seqId", grep("^SUPP", colnames(ins), value = T), "support")

probRatio = function (q1, q2) {
  if (q1 / q2 > 0) {
    log = log10(q1 / q2)
    if (is.infinite(log) || is.nan(log)) return (0)
    else return(log)
    }
  else return(0)
}
library(Rcpp)
sourceCpp("/path/to/insSuport.cpp")

genotype = function (supportVal, coverageVal, 
                     chances = c(0.05, 0.5, 0.95), 
                     genotypes = c(0, 1, 2), 
                     genoNames = c("homRef", "het", "homAlt"),
                     maxNorm = 250) {
  if(supportVal > coverageVal) {
    coverageVal = supportVal 
  }
  maxValue = max(supportVal, coverageVal)
  if (maxValue > maxNorm) {
    factorNorm = 250 / maxValue
    supportVal = supportVal * factorNorm
    coverageVal = coverageVal * factorNorm
  }
  
  probs = lapply(chances, dbinomSnifflesCpp, x = supportVal, size = coverageVal) %>% unlist() %>% setNames(genoNames) %>% sort(decreasing = T)
  probs = probs / sum(probs)
  q1 = probs[1]
  q2 = probs[2]
  qual = min(60, floor(-10 * probRatio(q2, q1)))
  geno = q1 %>% attr("name")
  genoResult = (if (geno == "homRef") 0 else if (geno == "het") 1 else 2)
  c(genoResult, qual)
}
genotypeVec = Vectorize(genotype, vectorize.args = c("supportVal", "coverageVal"))

genotypes = mclapply(1:length(selectedPatients), function (i) {
  x = selectedPatients[i]
  support = unlist(lapply(ins[seqnames %in% chrs][order(seqId)]$support, `[`, ... = i))
  coverage = round(fread(paste0("/path/to/", x, ".minimap.ins.regions.bed.gz"),
                col.names = c("chr", "start", "end", "seqId", "coverage"))[order(seqId)]$coverage)
  result = t(genotypeVec(support, coverage))
  result
}, mc.cores = 10)

genotypesValues = mclapply(1:length(selectedPatients), function (i) {genotypes[[i]][, 1]}, mc.cores = 16)
genotypesColumn = genotypesValues %>% simplify2array() %>% asplit(1) %>% lapply(as.numeric)
if (seqids %>% unique() %>% length() != seqids %>% length()) {warning("Duplicated IDs in the genotyping!")}
genotyped = data.table(seqId = seqids, genotype = genotypesColumn)
genotyped = genotyped[ins, on = "seqId", nomatch = NULL]

genotyped = genotyped[repeat.percentage >= 0.85 &
                        repeat.class %in% c(class1, class2)]
genotyped[, SUPP_geno := genotype %>% lapply(`>=`, y = 1) %>% lapply(as.numeric) %>% lapply(sum) %>% unlist()]

# False positives of trusty
hasAlleleleMask = lapply(genotyped$genotype, `>=`, y = 1) # For each INS and patient, bool for (0/1 or 1/1)
trustyWithoutAllelesMask = mapply(function (x, y) {!x & y},
         hasAlleleleMask, genotyped$SUPP_mask_trusty,
         SIMPLIFY = F
         )
trustyWithoutAlleles = trustyWithoutAllelesMask %>% # was detected by trusty but genotyped as 0/0
   lapply(as.numeric) %>% lapply(sum) # sum of the cases where it happens per INS
falpostrusty = trustyWithoutAlleles >= 1
fpt = genotyped[falpostrusty]
allFalsePositivesIds = fpt[trustyWithoutAlleles[falpostrusty] %>% unlist() == fpt$SUPP_trusty]$seqId

observedMafTable = copy(genotyped)
observedMafTable[, maf := genotype %>% lapply(function(x) sum(x) / (length(selectedPatients) * 2 ) ) %>% unlist() ] # humans are diploid, thus the length() * 2

# We remove the ones that were genotyped as 0/0 on all cases
genotyped[!seqId %in% allFalsePositivesIds]
repeats = repeats[!seqId %in% allFalsePositivesIds]
ins = ins[!seqId %in% allFalsePositivesIds]
# ----

# 
# Genes
# 
# ----
# Lets prepare a list of the genes affected by our INS
genes = repeats[SUPP_trusty > 0 & annot.subclass %in% c("cds", "introns", "promoters") & !is.na(annot.symbol), 
        .SD[1], 
        by = c("seqId", "annot.symbol")][order(annot.symbol)]
genes = genes[, .(
  uniques = seqId %>% unique() %>% length(),
  effectLax = sum(SUPP_min3),
  effectStrict = sum(SUPP_trusty)
  ), by = annot.symbol][order(-effectStrict)]

repetitiveRef = fread("/path/to/repeatsReferenceTENoOverlap.bed", select = c(1:3,5,6), col.names = c("seqnames", "start", "end", "name", "fam"))
setkey(repetitiveRef, seqnames, start, end)

# Annotation from genes: we need both promoter and 3UTR
genesFromAnnotation = data.table(data.frame(readRDS("/path/to/annotationsMulti.rds")))[startsWith(type, "hg38_genes")]
# It turns out not all genes (including many in our study) have 3UTR
genesCoords = genesFromAnnotation[!is.na(gene_id) & !is.na(symbol) & 
  type != "hg38_genes_1to5kb",
                    .(
                      seqnames = seqnames[1],
                      start = min(start),
                      end = max(end)
                    ), by = symbol]
genesCoords[, length := end - start + 1]
setcolorder(genesCoords, c("seqnames", "start", "end", "symbol", "length"))
setkey(genesCoords, seqnames, start, end)
# Lets calculate the overlap
overlap = foverlaps(repetitiveRef, genesCoords, nomatch = NULL)
overlap[, `:=`(
  i.start = pmax(i.start, start),
  i.end = pmin(i.end, end)
)][, repeatSpan := i.end - i.start + 1]

# This is finally the result of genes with % covered by repeats.
geneOverlapRepeat = overlap[, 
        .(seqnames = seqnames[1],
          start = start[1], end = end[1],
          length = length[1],
          repeatSpan = sum(repeatSpan)), 
        by = symbol][, repeatPercent := repeatSpan / length * 100]
data.table::setnames(geneOverlapRepeat, "symbol", "annot.symbol")

goiRepeatPercent = genes[geneOverlapRepeat, on = "annot.symbol"][!is.na(uniques)][order(-uniques)]
gnoiRepeatPercent = genes[geneOverlapRepeat, on = "annot.symbol"][is.na(uniques)][order(-uniques)]

t = t.test(goiRepeatPercent$repeatPercent, gnoiRepeatPercent$repeatPercent)
t2 = t.test(goiRepeatPercent[uniques == 1]$repeatPercent, goiRepeatPercent[uniques > 1]$repeatPercent)
# ----

# This is also used latter, liberally
roi = repeats[repeat.class %in% c("Retroposon", "SINE", "LINE", "DNA", "LTR") & SUPP_trusty > 0 & repeat.percentage >= 0.85] # repeats of interest for the article
ioi = roi[, .SD[1], by = seqId]

# Insertions' targets
# ----
# A smaller version of ioi
ioi2 = ioi[,
    .(seqId, seqnames, start, end, name, repeat.class, repeat.subclass)]
setkey(ioi2, seqnames, start, end)

# This is a very inefficient step
# which takes a significant amount of time
# to do something practically useless
if (T) {
fams = repetitiveRef$fam %>%
  stri_split_fixed(",") %>% 
  mclapply(function (x) {
    sub(pattern = "([^/]+)/.*", replacement = "\\1", perl = T, x = x) %>% 
      unique() %>% 
      paste0(collapse = ",")
  }, mc.cores = 20)
subfams =
  repetitiveRef$fam %>%
  stri_split_fixed(",") %>% 
  mclapply(function (x) {
    sub(pattern = "[^/]+/(.*)", replacement = "\\1", perl = T, x = x) %>% 
      unique() %>% 
      paste0(collapse = ",")
  }, mc.cores = 20)
repetitiveRef[, c("fam", "subfam") := .(fams, subfams)]
saveRDS(fams, "/path/to/refProcessedFam.rds")
saveRDS(subfams, "/path/to/refProcessedSubfam.rds")
} else {
  repetitiveRef[
    , c("fam", "subfam") := .(
    readRDS("/path/to/refProcessedFam.rds"), 
    readRDS("/path/to/refProcessedSubfam.rds"))
    ]
}
  
o = foverlaps(repetitiveRef, ioi2, nomatch = NULL)
o[, c("start", "end", "seqnames", "i.start", "i.end") := NULL]

o2 = o[, .(
    inRepeats = .SD %>% nrow(),
    sameFam = .SD[stri_detect_regex(str = fam, pattern = paste0("(^|,)",repeat.class,"(,|$)"))] %>% nrow(),
    sameSubfam = .SD[stri_detect_regex(
      str = subfam, 
      pattern = paste0("(^|,)", repeat.subclass,"(\\?|(-[^,]*))?(,|$)"))] %>% nrow(),
    sameElement = .SD[stri_detect_regex(str = i.name, pattern = paste0("(^|,)",name,"(,|$)"))] %>% nrow() # Since name can be "Name1,Name2":
  ), by = c("repeat.class", "repeat.subclass")
  ][
    ioi2[, .(total = .N), by = repeat.subclass], on = "repeat.subclass", nomatch = NULL # join to get total count
    ] %T>% setcolorder(
      c("repeat.class", "repeat.subclass", "total", "inRepeats", "sameFam", "sameSubfam", "sameElement")
      ) %>% .[
    , .(
      class = repeat.class,
      subclass = repeat.subclass,
      total,
      inRepeats = inRepeats,
      sameFam = sameFam,
      sameSubfam = sameSubfam,
      sameElement = sameElement
    )
    ] %T>% setorder(class, subclass)
# ----

# 
# Comparisons against other datasets
# 
# ----
vcf = ioi[observedMafTable[, c("seqId", "SUPP_geno", "maf")], on = "seqId", nomatch = NULL]
vcf = vcf[ins[, c("seqId", "vcf_alt")], on = "seqId", nomatch = NULL]
vcf$seqId %>% length() == vcf$seqId %>% unique() %>% length()
vcf = vcf[, .(
  `#CHROM` = seqnames,
  POS = start,
  ID = seqId,
  REF = "N",
  ALT = vcf_alt,
  QUAL = ".",
  FILTER = ".",
  INFO = paste(
    "SVTYPE=INS",
    paste0("METYPE=", repeat.subclass),
    paste0("SVLEN=", SVLEN),
    paste0("END=", start),
    paste0("SUPP_strict=", SUPP_trusty),
    paste0("SUPP_lax=", SUPP_min3),
    paste0("SUPP_geno=", SUPP_geno),
    paste0("MAF=", maf),
    sep = ";")
)]
vcf[stri_detect_fixed(INFO, "METYPE=Alu;")]
fwrite(vcf, "/path/to/all.me.txt", sep = "\t")
system("bash /path/to/sets/prepareOurMe.sh")
system("bash /path/to/runComparisons.sh")

vsIndigen = fread("/path/to/us.vs.indigen.vcf.gz",
                   skip = "#CHROM")
vsIndigen[
  ,  c( "len.us", "pos.us") := US %>% stri_split_fixed(pattern = ":") %>% lapply(`[`, y = c(3, 11)) %>% transpose()
  ][]
vsIndigen[
  ,  c( "len.indigen", "pos.indigen") := Indigen %>% stri_split_fixed(pattern = ":") %>% lapply(`[`, y = c(3, 11)) %>% transpose()
  ]
vsIndigen[, SUPP := sub(x = INFO, pattern = ".*SUPP=([^;\t]+);.*", replacement = "\\1", perl = T) %>% as.numeric()]
vsIndigen[ # This will introduce NA by coercion: "." -> NA
    , `:=`(
    len.us = as.numeric(len.us),
    pos.us = sub(pos.us, pattern = ".*_([^_]+)_[^_]+$", replacement = "\\1") %>% as.numeric(),
    len.indigen = as.numeric(len.indigen),
    pos.indigen = sub(pos.indigen, pattern = ".*_([^_]+)_[^_]+$", replacement = "\\1") %>% as.numeric()
  )
]
vsIndigen[, seqId := sub(pattern = ".*:(seq[[:digit:]]+|\\.):.*", replacement = "\\1", x = US)]
vsIndigen[SUPP == 2, `:=`(lendiff = abs(len.us - len.indigen), posdiff = abs(pos.us - pos.indigen))] 

vsChaisson = fread("/path/to/us.vs.hgsvc.vcf.gz",
                   skip = "#CHROM")
vsChaisson[
  ,  c( "len.us", "pos.us") := US %>% stri_split_fixed(pattern = ":") %>% lapply(`[`, y = c(3, 11)) %>% transpose()
  ][]
vsChaisson[
  ,  c( "len.indigen", "pos.indigen") := HGSVC %>% stri_split_fixed(pattern = ":") %>% lapply(`[`, y = c(3, 11)) %>% transpose()
  ]
vsChaisson[, SUPP := sub(x = INFO, pattern = ".*SUPP=([^;\t]+);.*", replacement = "\\1", perl = T) %>% as.numeric()]
vsChaisson[ # This will introduce NA by coercion: "." -> NA
    , `:=`(
    len.us = as.numeric(len.us),
    pos.us = sub(pos.us, pattern = ".*_([^_]+)_[^_]+$", replacement = "\\1") %>% as.numeric(),
    len.indigen = as.numeric(len.indigen),
    pos.indigen = sub(pos.indigen, pattern = ".*_([^_]+)_[^_]+$", replacement = "\\1") %>% as.numeric()
  )
]
vsChaisson[, seqId := sub(pattern = ".*:(seq[[:digit:]]+|\\.):.*", replacement = "\\1", x = US)]
vsChaisson[SUPP == 2, `:=`(lendiff = abs(len.us - len.indigen), posdiff = abs(pos.us - pos.indigen))] 
# ----

# 
# Deletions
# 
deletedMeMin3 = fread("/path/to/all.merged2.min3.deletedMe.vcf.gz", skip = "#CHROM")

# This removes "LTR?" elements and similar ones
deletedMeMin3 = deletedMeMin3[!INFO %like% "ME=[^;]*\\?[^;]*" & !INFO %like% "MEFAM=[^;]*\\?[^;]*"]
deletedMeMin3 = deletedMeMin3[`#CHROM` %in% chrs]

deletedMeMin3[ # Ids from survivor output are not guaranted to be unique, and indeed they are not
  , `:=`(
    survivorId = ID,
    ID = paste0("surv", 1:nrow(deletedMeMin3)),
    svlen = INFO %>% sub(".*SVLEN=-([0-9]+).*", "\\1", x = .) %>% as.numeric()
  )
]
delMeMin3Genotypes = deletedMeMin3[
  , .SD,
  .SDcols = colnames(deletedMeMin3)[c(3, 10:(10 + length(selectedPatients) ) )]
]

# Add columns
# ----
# The read supports gets lost in the intersample merge, but we since the file is min3,
# so all instances are min3, and if PSV != NaN, RE >= 3 -> min3 for that instance, since it existed
# then, if PSV == "11" -> trusty, since min3 is already checked
genoFields = deletedMeMin3[1, 9] %>% unname() %>% stri_split_fixed(pattern = ":") %>% unlist()
genoFieldsOfInt = c("PSV")
pos = which(genoFields %in% genoFieldsOfInt)
for (col in colnames(deletedMeMin3[1, .SD, .SDcols = patterns("^Sample")]) ) {
  deletedMeMin3[, paste0(col, "_", genoFieldsOfInt) := tstrsplit(.SD[[col]], ":", fixed = T, keep = pos)]
}
deletedMeMin3[, SUPP_VEC_min3 := apply(.SD, 2, `!=`, y = "NaN") %>% apply(., 2, as.integer) %>% apply(., 1, paste0, collapse = ""), .SDcols = grep(pattern = "_PSV$", colnames(deletedMeMin3))]
deletedMeMin3[, SUPP_min3 := stri_count_fixed(SUPP_VEC_min3, "1")]
deletedMeMin3[, SUPP_VEC_trusty := apply(.SD, 2, `==`, y = "11") %>% apply(., 2, as.integer) %>% apply(., 1, paste0, collapse = ""), .SDcols = patterns("_PSV$")]
deletedMeMin3[, SUPP_trusty := stri_count_fixed(SUPP_VEC_trusty, "1")]
deletedMeMin3[, SVLEN := as.numeric(gsub(".*SVLEN=([^;]+);.*", "\\1", INFO)) %>% abs()]
deletedMeMin3[, grep("^(Sample|FORMAT)", colnames(deletedMeMin3)) := NULL]
deletedMeMin3 = deletedMeMin3[
  , .(
    seqnames = `#CHROM`,
    start = POS,
    end = sub(pattern = ".*;END=([^;]+).*", replacement = "\\1", x = INFO, perl = T) %>% as.numeric(),
    id = ID,
    SUPP_min3,
    SUPP_trusty,
    svlen,
    # info = INFO,
    name = sub(pattern = ".*ME=([^;]+).*", replacement = "\\1", x = INFO, perl = T),
    repeat.class = sub(pattern = ".*MEFAM=([^;]+).*", replacement = "\\1", x = INFO, perl = T),
    repeat.coords = sub(pattern = ".*MECOORDS=([^;]+).*", replacement = "\\1", x = INFO, perl = T)
  )
][, c("repeat.class", "repeat.subclass") := repeat.class %>% 
  stri_match_all_regex("([^/]+)(?:/(.+))?$") %>% 
  rapply(., f = `[`, ... = 2:3, how = "r") %>% transpose()][
    , repeat.subclass := fifelse(is.na(repeat.subclass), "-", repeat.subclass)
  ][
    , c("repeat.subclass", "repeat.family") := 
        stri_match_all_regex(repeat.subclass, "([^-]+|-)(?:-(.+))?$") %>% 
        rapply(., f = `[`, ... = 2:3, how = "r") %>% transpose()
  ][, repeat.family := fifelse(is.na(repeat.family), "-", repeat.family)]
# ----

mosdepthDelBed = deletedMeMin3[, c("seqnames", "start", "end", "id")]
fwrite(mosdepthDelBed, file = "/path/to/teDeletions.bed", sep = "\t", col.names = F, scipen=50)

system("bash slurm/mosdepth.slurm config/config.sh DEL")

delMeMin3Ids = delMeMin3Genotypes[, .(ID)]
# will create a column of vectors
# position 1 is id por P1 in their file
delMeMin3Ids[, ids := delMeMin3Genotypes[, 2:(ncol(delMeMin3Genotypes) - 1)] %>%
  apply(2, function (x) {
    x %>% 
      stri_split_fixed(pattern = ":") %>% lapply(`[`, y = 8) %>% 
      stri_split_fixed(pattern = ",") %>% lapply(`[`, y = 1)
  }) %>% simplify2array() %>% asplit(1) %>% lapply(unlist) %>% lapply(unname)]

supports = mclapply(1:length(selectedPatients), function (i) {
  patient = selectedPatients[i]
  folder = folders[selectedPatientsMask][i]
  
  file = paste0("/path/to/", folder, "/variant_calling/GRCh38/",
                patient, ".minimap.hg38.merged2.min3.vcf")
  vcf = fread(file, skip = "#CHROM")
  setnames(
    vcf, 
    c("#CHROM", "ID", "POS"), 
    c("seqnames", "id", "start"))
  vcf[
    , `:=`(
      svtype = gsub(".*SVTYPE=([^;]+).*", "\\1", vcf$INFO, perl = T)
    )]
  idTable = delMeMin3Ids[ # A table with ids of merge and Pi's
    , .(
      ID = ID,
      ids = ids %>% lapply(`[`, y = i) %>% unlist()
    )
  ][order(ID)] # alphanumerical sort to be coherent with other tables
  mergedIdsPresent = idTable[ids != "NaN", ID] # The merged ids of deletions Pi has
  patientIds = idTable$ids %>% .[. != "NaN"] %>% unlist()  # The ids from Pi on their own file
  supportVector = vcf[id %in% patientIds, c(10, 11)] %>% # We take the Genotype columns from both callers in Pi file
    apply(2, function (x) {
      x %>% 
        stri_split_fixed(pattern = ":") %>% 
        lapply(`[`, y = 4) %>% # We select the field with read support
        stri_split_fixed(pattern = ",") %>% 
        lapply(`[`, y = 2) %>% unlist() %>% # and we get the support for the SV
        as.numeric()
      }) %>% 
    apply(1, max) # Then we get the support from the caller that returns the highest
  result = merge(
    idTable, 
    data.table(ID = mergedIdsPresent, support = supportVector),
    by = "ID",
    all.x = T)
  result[, `:=`(p = paste0("P", i), support = fifelse(is.na(support), 0, support))]
  result
}, mc.cores = 5) # each element comes out sorted alphanumerically by ID (format: "survX")
supports = supports %>% rbindlist()
supports = supports %>% dcast(ID ~ p, value.var = "support");
setcolorder(supports, colnames(supports) %>% stri_sort(numeric = T))
supportsCol = supports[, -"ID"] %>% simplify2array() %>% asplit(1) %>% lapply(unname) # still sorted

seqids = mosdepthDelBed[order(id)]$id # its sorted alphanumerically
coveragesForGenotypingDel = mclapply(selectedPatients, function (x) {
  fread(paste0("/path/to/", x, ".minimap.del.regions.bed.gz"),
                col.names = c("chr", "start", "end", "seqId", "coverage"))[order(seqId)]$coverage # also sorted alphanumerically
}, mc.cores = 1)
# coveragesForGenotypingDel[[13]][1:5]
# seqids %>% length()
# lapply(asplit(simplify2array(coveragesForGenotypingDel), 1), as.numeric) %>% length()
# supportsCol %>% length()
coverageForGenotypingDel = data.table(
  seqId = seqids, # has been sorted alphanumerically
  coverage = lapply(asplit(simplify2array(coveragesForGenotypingDel), 1), as.numeric),
  support = supportsCol
)
coverageForGenotypingDel[, totalReads := mapply(`+`, coverage, support, SIMPLIFY = F)]
genotypesDel = mclapply(1:length(selectedPatients), function (i) {
  # x = selectedPatients[i]
  
  support = coverageForGenotypingDel$support %>% lapply(`[`, y = i) %>% unlist()
  reads = coverageForGenotypingDel$totalReads %>% lapply(`[`, y = i) %>% unlist()
  
  result = t(genotypeVec(support, reads))
  result
}, mc.cores = 10)

genotypesValuesDel = mclapply(1:length(selectedPatients), function (i) {genotypesDel[[i]][, 1]}, mc.cores = 16)
genotypesDelColumn = genotypesValuesDel %>% simplify2array() %>% asplit(1) %>% lapply(as.numeric)
genotypedDel = deletedMeMin3[data.table(id = seqids, genotype = genotypesDelColumn), on = "id"]
alleleN = (genotypedDel[1, genotype] %>% unlist() %>% length() ) * 2
genotypedDel[, maf := (genotype %>% lapply(sum) %>% unlist() ) / alleleN]

genotypesValuesDel = mclapply(1:length(selectedPatients), function (i) {genotypesDel[[i]][, 1]}, mc.cores = 16)
genotypesDelColumn = genotypesValuesDel %>% simplify2array() %>% asplit(1) %>% lapply(as.numeric)
genotypedDel = deletedMeMin3[data.table(id = seqids, genotype = genotypesDelColumn), on = "id"]
alleleN = (genotypedDel[1, genotype] %>% unlist() %>% length() ) * 2
genotypedDel[, maf := (genotype %>% lapply(sum) %>% unlist() ) / alleleN]

# Annotation for deletions
delAnno = annotatr::annotate_regions(regions = GenomicRanges::makeGRangesFromDataFrame(
      deletedMeMin3[!id %in% genotypedDel[maf == 0]$id], 
      seqnames.field = "seqnames", 
      keep.extra.columns = T),
     annotations = hg38Anno,
     ignore.strand = F,
     quiet = T
    )
delAnno = delAnno %>% as.data.frame(stringsAsFactors = F) %>% data.table()
delAnno[, annot.type := stri_replace_first_regex(annot.type, "hg38_(.*)", "$1")]
delAnno = delAnno[stri_detect_fixed(annot.type, pattern = "genes_")]
text = stri_match_all_regex(delAnno$annot.type, "([^-]*)_(.*)") %>% 
  lapply(., function(x) {x[2:3]})
delAnno[, `:=`(annot.class = text %>% lapply(`[`, 1) %>% unlist(),
               annot.subclass = text %>% lapply(`[`, 2) %>% unlist()
               )]
delAnno[, annot.type := NULL] # Redundant
delAnno[, repeat.class := as.factor(repeat.class)]
rm(text); gc()
delAnno[, c("annot.strand", "annot.seqnames", "annot.width", "annot.id", "annot.tx_id") := NULL]

saveRDS(covs, "/path/to/covs.rds")
saveRDS(covDist, "/path/to/ins/covDist.rds")
saveRDS(ins, "/path/to/ins/ins.rds")
saveRDS(allIns, "/path/to/ins/allIns.rds")
saveRDS(repeats, "/path/to/ins/repeats.rds")
saveRDS(genes, "/path/to/ins/genes.rds")
saveRDS(goiRepeatPercent, "/path/to/ins/goi.rds")
saveRDS(gnoiRepeatPercent, "/path/to/ins/gnoi.rds")
saveRDS(t, "/path/to/ins/t.rds")
saveRDS(t2, "/path/to/ins/t2.rds")
geneAnno = data.table(as.data.frame(hg38Anno))
saveRDS(observedMafTable, "/path/to/ins/observedMafTable.rds")
saveRDS(geneAnno, "/path/to/ins/geneAnno.rds")
saveRDS(genotyped, "/path/to/ins/genotyped.rds")
saveRDS(o2, "/path/to/ins/o2.rds")
saveRDS(vsIndigen, "/path/to/ins/vsIndigen.rds")
saveRDS(vsChaisson, "/path/to/ins/vsChaisson.rds")
saveRDS(deletedMeMin3, "/path/to/ins/deletedMeMin3.rds")
saveRDS(genotypedDel, "/path/to/ins/genotypedDel.rds")
saveRDS(delAnno, "/path/to/ins/delAnno.rds")