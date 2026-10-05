rm(list=ls())
library(dplyr)
library(tidyr)
library(Seurat)
library(patchwork)
library(data.table)
library(Matrix)
library(hdf5r)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggsignif)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(GEOquery)
library(VennDiagram)
library(grid)
library(tibble)
base_dir <- "E:/sc_RNA_seq/PGK1_PD1"
paths <- list(
  GSE120575 = file.path(base_dir, "mel", "gse120575"),
  GSE243013 = file.path(base_dir, "NSCLC", "GSE243013"),
  GSE205506 = file.path(base_dir, "CRC", "GSE205506"),
  GSE206325 = file.path(base_dir, "HCC", "GSE206325")
)
first_existing <- function(x) {
  x <- x[file.exists(x)]
  if (length(x) == 0) return(NA_character_)
  x[1]
}
find_col <- function(df, patterns) {
  nms <- colnames(df)
  hit <- unique(unlist(lapply(patterns, function(p) grep(p, nms, ignore.case = TRUE, value = TRUE))))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}
normalize_response <- function(x) {
  y <- trimws(as.character(x))
  z <- tolower(y)
  out <- rep(NA_character_, length(y))
  out[grepl("non[- _]?responder|non[- _]?response|non[- _]?mpr|non[- _]?pcr|\\bnr\\b|stable disease|progressive disease|\\bsd\\b|\\bpd\\b", z)] <- "NR"
  out[is.na(out) & grepl("responder|complete response|partial response|major pathological response|pathological complete response|\\bmpr\\b|\\bpcr\\b|\\bcr\\b|\\bpr\\b|^r$", z)] <- "R"
  out
}
standardize_metadata <- function(obj, dataset, cancer) {
  md <- obj@meta.data
  md$dataset <- dataset
  md$cancer <- cancer
  patient_col <- find_col(md, c("patient", "subject", "donor", "case", "individual"))
  sample_col <- find_col(md, c("sample", "orig.ident", "library", "specimen"))
  response_col <- find_col(md, c("response", "responder", "outcome", "efficacy", "pathological", "pathologic", "mpr", "pcr"))
  if (!is.na(patient_col)) md$patient_id <- as.character(md[[patient_col]]) else if (!is.na(sample_col)) md$patient_id <- as.character(md[[sample_col]]) else md$patient_id <- rownames(md)
  if (!is.na(sample_col)) md$sample_id <- as.character(md[[sample_col]]) else md$sample_id <- md$patient_id
  if (!is.na(response_col)) {
    md$response_raw <- as.character(md[[response_col]])
    md$response <- normalize_response(md$response_raw)
  } else {
    md$response_raw <- NA_character_
    md$response <- NA_character_
  }
  obj@meta.data <- md
  obj
}
load_GSE120575 <- function(path) {
  data <- read.table(file.path(path, "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"), sep = "\\t", header = TRUE, quote = "", fill = TRUE, check.names = FALSE)
  data <- data[-1, -1]
  patient_id <- fread(file.path(path, "GSE120575_patient_ID_single_cells.txt.gz"), sep = "\\t", header = TRUE)
  colnames(patient_id) <- as.character(patient_id[19,])
  patient_id <- as.data.frame(patient_id[20:16310, 1:11])
  rownames(patient_id) <- patient_id$title
  obj <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200, project = "GSE120575")
  common_cells <- intersect(colnames(obj), rownames(patient_id))
  obj <- subset(obj, cells = common_cells)
  obj <- AddMetaData(obj, patient_id[common_cells, , drop = FALSE])
  new_colnames <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "sample_name", "title", "source_name", "organism", "patient_id", "response", "therapy", "molecule_description", "processed_data_file", "raw_file")
  if (ncol(obj@meta.data) == length(new_colnames)) colnames(obj@meta.data) <- new_colnames
  if ("therapy" %in% colnames(obj@meta.data)) obj <- subset(obj, subset = therapy == "anti-PD1")
  if ("response" %in% colnames(obj@meta.data)) obj$response <- dplyr::case_when(obj$response == "Responder" ~ "R", obj$response == "Non-responder" ~ "NR", TRUE ~ as.character(obj$response))
  obj$dataset <- "GSE120575"
  obj$cancer <- "Melanoma"
  if (!"sample_id" %in% colnames(obj@meta.data)) obj$sample_id <- obj$patient_id
  obj$response_raw <- obj$response
  obj
}
load_GSE243013 <- function(path) {
  mtx_file <- first_existing(c(file.path(path, "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz"), list.files(path, pattern = "NSCLC_immune_scRNA_counts.*mtx", full.names = TRUE, recursive = TRUE)))
  gene_file <- first_existing(c(file.path(path, "GSE243013_genes.csv.gz"), list.files(path, pattern = "genes.*csv", full.names = TRUE, recursive = TRUE)))
  barcode_file <- first_existing(c(file.path(path, "GSE243013_barcodes.csv.gz"), list.files(path, pattern = "barcodes.*csv", full.names = TRUE, recursive = TRUE)))
  meta_file <- first_existing(c(file.path(path, "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz"), list.files(path, pattern = "NSCLC_immune_scRNA_metadata.*csv", full.names = TRUE, recursive = TRUE)))
  mat <- readMM(mtx_file)
  genes <- fread(gene_file, header = FALSE)
  barcodes <- fread(barcode_file, header = FALSE)
  meta <- fread(meta_file)
  gene_names <- make.unique(as.character(genes[[ncol(genes)]]))
  cell_names <- as.character(barcodes[[ncol(barcodes)]])
  if (nrow(mat) == length(cell_names) && ncol(mat) == length(gene_names)) mat <- t(mat)
  rownames(mat) <- gene_names
  colnames(mat) <- cell_names
  obj <- CreateSeuratObject(counts = mat, min.cells = 3, min.features = 200, project = "GSE243013")
  meta <- as.data.frame(meta)
  barcode_col <- find_col(meta, c("barcode", "cell"))
  if (!is.na(barcode_col)) rownames(meta) <- make.unique(as.character(meta[[barcode_col]])) else if (nrow(meta) == ncol(obj)) rownames(meta) <- colnames(obj)
  common <- intersect(colnames(obj), rownames(meta))
  if (length(common) > 0) {
    obj <- subset(obj, cells = common)
    obj <- AddMetaData(obj, meta[common, , drop = FALSE])
  }
  standardize_metadata(obj, "GSE243013", "NSCLC")
}
load_GSE205506 <- function(path) {
  mtx_files <- list.files(path, pattern = "matrix\\.mtx\\.gz$", full.names = TRUE, recursive = TRUE)
  if (length(mtx_files) == 0) stop("GSE205506目录中未找到matrix.mtx.gz文件")
  geo_meta <- tryCatch({
    g <- getGEO("GSE205506", GSEMatrix = TRUE, getGPL = FALSE)
    do.call(rbind, lapply(g, pData))
  }, error = function(e) NULL)
  objs <- lapply(mtx_files, function(mtx) {
    feature <- sub("matrix\\.mtx\\.gz$", "features.tsv.gz", mtx)
    barcode <- sub("matrix\\.mtx\\.gz$", "barcodes.tsv.gz", mtx)
    if (!file.exists(feature) || !file.exists(barcode)) return(NULL)
    sample_name <- sub("_matrix\\.mtx\\.gz$", "", basename(mtx))
    sample_name <- sub("_v[23]$", "", sample_name)
    mat <- ReadMtx(mtx = mtx, features = feature, cells = barcode, feature.column = 2, unique.features = TRUE)
    obj <- CreateSeuratObject(mat, min.cells = 3, min.features = 200, project = sample_name)
    obj$sample_id <- sample_name
    gsm <- sub("_.*$", "", sample_name)
    obj$geo_accession <- gsm
    if (!is.null(geo_meta) && "geo_accession" %in% colnames(geo_meta) && gsm %in% geo_meta$geo_accession) {
      one <- geo_meta[match(gsm, geo_meta$geo_accession), , drop = FALSE]
      for (nm in colnames(one)) obj[[paste0("geo_", make.names(nm))]] <- as.character(one[[nm]][1])
    }
    obj
  })
  objs <- Filter(Negate(is.null), objs)
  if (length(objs) == 1) obj <- objs[[1]] else obj <- merge(objs[[1]], y = objs[-1], add.cell.ids = paste0("S", seq_along(objs)), project = "GSE205506")
  standardize_metadata(obj, "GSE205506", "CRC")
}
load_GSE206325 <- function(path) {
  rda_file <- first_existing(c(file.path(path, "GSE206325_data_SingleCellExperiment_object.Rda.gz"), list.files(path, pattern = "SingleCellExperiment.*Rda", full.names = TRUE, recursive = TRUE)))
  env <- new.env(parent = emptyenv())
  con <- if (grepl("\\.gz$", rda_file)) gzfile(rda_file, "rb") else rda_file
  load(con, envir = env)
  if (inherits(con, "connection")) close(con)
  candidates <- mget(ls(env), envir = env)
  sce <- candidates[[which(vapply(candidates, function(x) inherits(x, "SingleCellExperiment"), logical(1)))[1]]]
  if ("counts" %in% assayNames(sce)) {
    obj <- CreateSeuratObject(counts = assay(sce, "counts"), meta.data = as.data.frame(colData(sce)), min.cells = 3, min.features = 200, project = "GSE206325")
  } else {
    obj <- as.Seurat(sce)
  }
  ann_file <- first_existing(c(file.path(path, "GSE206325_sample_annots_Liver_Treated_patients.csv.gz"), list.files(path, pattern = "sample_annots.*csv", full.names = TRUE, recursive = TRUE)))
  if (!is.na(ann_file)) {
    ann <- as.data.frame(fread(ann_file))
    obj_col <- find_col(obj@meta.data, c("sample", "orig.ident", "patient"))
    ann_col <- find_col(ann, c("sample", "orig.ident", "patient"))
    if (!is.na(obj_col) && !is.na(ann_col)) {
      idx <- match(as.character(obj@meta.data[[obj_col]]), as.character(ann[[ann_col]]))
      add <- ann[idx, , drop = FALSE]
      rownames(add) <- colnames(obj)
      colnames(add) <- paste0("ann_", make.names(colnames(add)))
      obj <- AddMetaData(obj, add)
    }
  }
  standardize_metadata(obj, "GSE206325", "HCC")
}
seurat_list <- list(
  GSE120575 = load_GSE120575(paths$GSE120575),
  GSE243013 = load_GSE243013(paths$GSE243013),
  GSE205506 = load_GSE205506(paths$GSE205506),
  GSE206325 = load_GSE206325(paths$GSE206325)
)
cancer_names <- c(GSE120575 = "Melanoma", GSE243013 = "NSCLC", GSE205506 = "CRC", GSE206325 = "HCC")
run_analysis <- function(all, dataset_name) {
  output_dir <- file.path(base_dir, "results", paste0(dataset_name, "_", cancer_names[[dataset_name]]))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(output_dir)
  all <- subset(all, subset = !is.na(response) & response %in% c("R", "NR"))
  all <- NormalizeData(all)
  all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
  all <- ScaleData(all)
  all <- RunPCA(all, features = VariableFeatures(all))
  pdf("ElbowPlot.pdf", width = 6, height = 5)
  print(ElbowPlot(all))
  dev.off()
  all <- FindNeighbors(all, dims = 1:15)
  all <- FindClusters(all, resolution = 0.5)
  all <- RunUMAP(all, dims = 1:15)
  p1 <- DimPlot(all, reduction = "umap", label = TRUE, label.size = 5)
  p2 <- DimPlot(all, reduction = "umap", label = TRUE, label.size = 5, group.by = "patient_id")
  ggsave("UMAP_clusters.pdf", p1, width = 7, height = 6)
  ggsave("UMAP_patient.pdf", p2, width = 8, height = 6)
  ref2 <- BlueprintEncodeData()
  all <- JoinLayers(all)
  test_data <- GetAssayData(all, layer = "data")
  singler_results2 <- SingleR(test = test_data, ref = ref2, labels = ref2$label.main, clusters = NULL, assay.type.test = "logcounts")
  all$singler_labels2 <- singler_results2$labels
  saveRDS(all, paste0(dataset_name, "_all.rds"))
  myloidcell<-subset(all,subset = singler_labels2=="Monocytes"|singler_labels2=="Macrophages"|singler_labels2=="DC")
  expr_data <- FetchData(
    myloidcell,
    vars = c("CD274", "PDCD1LG2", "response", "singler_labels2","patient_id")
  )
  expr_data$response <- factor(
    expr_data$response,
    levels = c("R", "NR")
  )
  plot_data <- melt(
    expr_data,
    id.vars = c("response", "singler_labels2"),
    measure.vars = c("CD274", "PDCD1LG2"),
    variable.name = "Gene",
    value.name = "Expression"
  )
  ggplot(plot_data, aes(x = singler_labels2, y = Expression, fill = response)) +
    geom_violin(
      position = position_dodge(0.8),
      scale = "width",
      width = 0.7,
      trim = TRUE,
      alpha = 0.8
    ) +
    geom_boxplot(
      width = 0.15,
      position = position_dodge(0.8),
      outlier.shape = NA,
      alpha = 0.8,
      color = "black"
    ) +
    facet_wrap(~ Gene, ncol = 2, scales = "free_y") +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12, color = "black")
    ) +
    stat_compare_means(
      aes(group = response),
      method = "wilcox.test",
      label = "p.signif",
      label.y = max(plot_data$Expression) * 1.15,
      size = 5,
      vjust = 0.5,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c("R" = "#E64B35", "NR" = "#4DBBD5")) +
    labs(
      x = "Myeloid Cell Type",
      y = "Expression Level",
      title = "PD-L1/PD-L2 Expression in Myeloid Subsets"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.title = element_text(face = "bold", size = 14),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      panel.border = element_blank()
    )
  ggsave("PDL1_L2_myeloid_expression.pdf", width = 6, height = 5, dpi = 300)
  avg_expr <- FetchData(
    myloidcell,
    vars = c("CD274", "PDCD1LG2", "response", "singler_labels2", "patient_id")
  ) %>%
    group_by(patient_id, singler_labels2, response) %>%
    summarise(
      CD274_avg = mean(CD274, na.rm = TRUE),
      PDCD1LG2_avg = mean(PDCD1LG2, na.rm = TRUE)
    ) %>%
    melt(
      id.vars = c("patient_id", "singler_labels2", "response"),
      variable.name = "Gene",
      value.name = "Expression"
    )
  avg_expr$response <- factor(
    avg_expr$response,
    levels = c("R", "NR")
  )
  ggplot(avg_expr, aes(x = singler_labels2, y = Expression, fill = response)) +
    geom_violin(
      width = 0.8,
      position = position_dodge(0.8),
      alpha = 0.3,
      trim = TRUE,
      scale = "width",
    )+
    geom_boxplot(
      width = 0.2,
      position = position_dodge(0.8),
      outlier.shape = NA,
      alpha = 0.8,
      color = "black"
    ) +
    geom_point(
      aes(group = response),
      position = position_jitterdodge(
        jitter.width = 0.15,
        dodge.width = 0.8
      ),
      shape = 21,
      color = "black",
      size = 2.5,
      alpha = 0.7,
      show.legend = FALSE
    ) +
    facet_wrap(~ Gene, ncol = 2, scales = "free_y",
               labeller = labeller(Gene = c("CD274_avg" = "PD-L1 (CD274)",
                                            "PDCD1LG2_avg" = "PD-L2 (PDCD1LG2)"))) +
    stat_compare_means(
      aes(group = response),
      method = "wilcox.test",
      label = "p.signif",
      size = 5,
      vjust = 0.5,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c("R" = "#E64B35", "NR" = "#4DBBD5")) +
    labs(
      x = "Myeloid Cell Type",
      y = "Average Expression per Sample",
      title = "Sample-level PD-L1/PD-L2 Expression in Myeloid Subsets"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.title = element_text(face = "bold", size = 14),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12)
    )
  ggsave("PDL1_L2_myeloid_expression_sample_avg.pdf", width = 6, height = 5, dpi = 300)
  avg_by_sample_overall <- FetchData(
    myloidcell,
    vars = c("CD274", "PDCD1LG2", "response", "patient_id")
  ) %>%
    group_by(patient_id, response) %>%
    summarise(
      CD274_avg = mean(CD274, na.rm = TRUE),
      PDCD1LG2_avg = mean(PDCD1LG2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    reshape2::melt(
      id.vars = c("patient_id", "response"),
      variable.name = "Gene",
      value.name = "Expression"
    )
  avg_by_sample_overall$Gene <- dplyr::recode(
    avg_by_sample_overall$Gene,
    "CD274_avg" = "PD-L1 (CD274)",
    "PDCD1LG2_avg" = "PD-L2 (PDCD1LG2)"
  )
  avg_by_sample_overall$response <- factor(avg_by_sample_overall$response, levels = c("R", "NR"))
  library(ggplot2)
  library(ggpubr)
  p_overall <- ggplot(avg_by_sample_overall, aes(x = response, y = Expression, fill = response)) +
    geom_violin(width = 0.9, alpha = 0.3, trim = TRUE) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.9, color = "black") +
    geom_point(
      position = position_jitter(width = 0.12, height = 0, seed = 1),
      shape = 21, color = "black", size = 2.5, alpha = 0.8, show.legend = FALSE
    ) +
    facet_wrap(~ Gene, ncol = 2, scales = "free_y") +
    stat_compare_means(
      method = "wilcox.test", label = "p.signif",
      size = 5, vjust = 0.5, show.legend = FALSE
    ) +
    scale_fill_manual(values = c("R" = "#E64B35", "NR" = "#4DBBD5")) +
    labs(
      x = "Response", y = "Average Expression per Sample (Myeloid total)"
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 12),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
    )
  print(p_overall)
  ggsave("PDL1_L2_myeloid_overall_sample_avg.pdf", p_overall, width = 6, height = 5, dpi = 300)
  exhaustion_related_NVE_U<-c(
    "ABCC2", "ABTB2", "ADAMTS2", "ADAR", "ADM", "ADORA2A", "AIM2", "AMPD3",
    "ANKLE2", "APOBEC3A", "APOL1", "ARAP2", "ARHGAP27", "ASAP2", "B4GALT5",
    "BAALC", "BAZ1A", "BHLHE22", "BIRC2", "BLVRA", "BLZF1", "LINC01588",
    "C15orf32", "SHFL", "C3orf38", "C4orf46", "CPLANE1", "SNHG15", "CCL4",
    "CCL5", "NOCT", "CD40", "CD83", "CDC14C", "CDKN2B", "CELA1", "CFLAR",
    "CLEC2D", "CPNE8", "CSRNP1", "CXCL9", "CXorf65", "DAPP1", "RIGI",
    "DDX60L", "DGKH", "DMRT2", "DNAAF1", "DUSP5", "EGFL7", "EGFLAM", "EHD1",
    "EHD4", "EPSTI1", "EREG", "ERRFI1", "ESPL1", "FAM174B", "FBP2", "FBXO40",
    "FCGR1BP", "FKSG41", "FUT4", "FXYD6", "GCH1", "GJB2", "GNGT1", "GPR183",
    "GPR32", "GRHL1", "GTPBP1", "HERC5", "HERC6", "HSF4", "IDO2", "IFI35",
    "IFIH1", "IL18R1", "IL18RAP", "IL19", "IL6", "IL7R", "IMPG1", "ISX",
    "KCNJ2", "WDR97", "KRT17", "KRT80", "LAG3", "SIPA1L1", "LINC01785",
    "MIR3945HG", "LYN", "MAP2K3", "MAP3K8", "MAP4K4", "MASP1", "MCTP1",
    "MEST", "MIR155HG", "MX2", "MXD1", "MYO1A", "NAMPT", "NCF1C", "NEDD4L",
    "NINJ1", "NPY1R", "NXF1", "OAS1", "OAS3", "OASL", "OTUD4", "P2RX4",
    "P2RX7", "PALM3", "TENT4A", "PARP14", "PELI1", "PI4K2B", "PIK3R3",
    "PIK3R5", "PIM3", "PIWIL4", "PLEK", "PLK3", "PML", "POPDC2", "PPFIBP2",
    "PPP1R15A", "PRG3", "HELZ2", "PRSS23", "PRX", "PSD3", "PTGER2", "PTGS2",
    "PTK2B", "PTPRF", "PVR", "NECTIN3", "RAB24", "RBBP6", "RBMXL1", "REN",
    "RIPK1", "RIPK2", "RNF144B", "RNF19A", "RNF213", "SAMD9L", "SCGB2A1",
    "SDC4", "SERPINB1", "SERPINE1", "SKIL", "SLAMF1", "SLAMF7", "SLC1A2",
    "SLC1A3", "SNN", "SORBS1", "SPAG1", "SPATA6", "STARD5", "STAT1", "STAT5A",
    "STX11", "TAAR5", "TDRD7", "TNS2", "TNIP2", "TNIP3", "TP53INP2", "TRIL",
    "TRIM22", "TSHB", "TSIX", "TSPO2", "TTN", "TYSND1", "UBE2L6", "USP18",
    "UST", "VCAN", "DNAI3", "WTAP", "XRN1", "ZBTB10", "ZC3HAV1", "ZNFX1")
  exhaustion_related_NVE_D<-c("AAK1", "ABCA1", "AICDA", "ALDH2", "ALG2", "ANKRD50", "ANXA2", "AP3M1", "ARID3B", "ATF2",
                              "ATG13", "ATG9A", "BMPR2", "BOD1L", "BRAF", "MPC2", "ODR4", "GLMP", "CYSTM1", "CABP4",
                              "CALU", "CENATAC", "CCNI", "CD38", "CD5L", "CDK2AP2", "CDKL2", "CHCHD10", "CHD2", "CHPT1",
                              "CLINT1", "CLK4", "CLN8", "CMPK1", "COG3", "COPA", "COPS3", "CR1L", "CRK", "CTSB",
                              "CUTA", "CYFIP2", "DDX17", "INTS6L", "DERL1", "DHRS1", "DIP2A", "DMWD", "DNPEP", "DPY19L3",
                              "DUSP3", "E2F2", "EBP", "ECH1", "ECI2", "EHD4", "ENDOU", "ENPP1", "EPM2AIP1", "ERGIC3",
                              "ERO1B", "ESYT2", "EVI5", "FAM117A", "FAM3C", "FAM53C", "FAM78A", "FARP2", "FBXO2", "FBXO4",
                              "FILIP1", "FLCN", "FOXO3", "FUCA2", "GALNT2", "GATA4", "GLB1", "GLT8D1", "GPATCH2", "GPCPD1",
                              "GRINA", "GRIP2", "GSK3B", "GSTT2", "GZF1", "HCFC2", "HINT3", "H2BC3", "HIVEP3", "HNRNPH3",
                              "HPSE", "IDH1", "IL18", "ING2", "INPP5K", "IQCB1", "ITPR3", "KDM5A", "AREL1", "KHDC4",
                              "LBH", "LHPP", "LIMS4", "LMBRD1", "LPIN2", "MAP3K1", "MBNL2", "MBTPS1", "MCCC2", "MIA3",
                              "MLLT3", "MNT", "MZB1", "N4BP2", "NCKAP1L", "NIPAL3", "NSF", "NUP210", "OTUD5", "PAIP2",
                              "PAQR3", "PDIA6", "PDK1", "PHF3", "PHYHD1", "PLOD2", "PNRC1", "PPP1R15A", "PPP2R5B", "SLC66A2",
                              "PRDX2", "PRDX5", "PRELID1", "PRR3", "PRRG4", "PSME4", "PSRC1", "PTGR1", "PTPRCAP", "RAB11A",
                              "RABGAP1L", "RAI1", "RAPGEF1", "RASA4", "RDH10", "RFFL", "RHOBTB1", "RPN2", "SAMHD1", "SCP2",
                              "SDF2", "SELENOM", "SELENOP", "SEPTIN7", "SERHL2", "SETX", "SF3B1", "SLAMF7", "SLC25A24", "SLC35A2",
                              "SLC44A1", "SLC48A1", "SMPDL3A", "SNAP23", "STK25", "SVIL", "SYS1", "TAPT1", "TAX1BP1", "TCP11L2",
                              "TECR", "TGFBR2", "THRA", "TJP2", "TM9SF1", "TMCO1", "TMED10", "TMED2", "TMED3", "NEMP1",
                              "TMEM9", "TOB1", "TPST1", "TRIM35", "TSPAN14", "TTC9C", "TXNIP", "UBAC2", "UBE2Q2", "UBN1",
                              "UHRF2", "UNC50", "USP22", "VCPIP1", "XRCC1", "YIPF4", "ZC3H11A", "ZMYND8", "ZNF420", "ZNF790")
  exhaustion_related_MemVE_U<-c(
    "ADAR", "AIDA", "ALS2", "ANKRA2", "ANKRD42", "ANXA1", "AP3M2", "ARHGAP25",
    "ARHGAP30", "ARHGEF6", "ARPP21", "ATP10A", "ATP8A1", "ATP8A2", "ATP8B4",
    "ADGRB2", "BATF2", "BCAT1", "BCL7B", "BMP2", "VSTM4", "NOXRED1", "CFAP126",
    "SUCO", "COA5", "OARD1", "LAMTOR4", "CARD6", "CASP1", "WASHC3", "CCDC90B",
    "CD244", "CD40", "CD83", "CD84", "CD99L2", "CEBPD", "CEP120", "CHMP4B",
    "CNN3", "CORO2A", "CSF2", "CSRNP1", "CTSG", "CTSS", "CXCR6", "TASL",
    "CYBB", "CYP4F3", "DAAM1", "DAPP1", "DAXX", "DCAF8", "DENND4A", "DOK1",
    "DSC2", "DUSP2", "DYNLRB1", "EPPK1", "ETNK1", "F2R", "UBALD2", "CALHM6",
    "FASLG", "FBXO4", "FBXO7", "FGF13", "FNBP1", "GBP4", "GBP7", "GCNT1",
    "GGPS1", "GIMAP6", "ADGRD1", "GPR18", "GPR65", "HHEX", "HINT3", "HMGN3",
    "HSD11B1", "HSH2D", "HSPA1A", "IFT25", "ICAM1", "ICOSLG", "IFIT2", "IFNG",
    "IL12RB1", "IL18", "IL18RAP", "IQSEC1", "IRF1", "IRF8", "ITGAL", "ITGB5",
    "JUN", "KIF9", "LAT2", "LCK", "LEAP2", "LFNG", "LGALS9B", "LIMD2", "LRRC15",
    "LRRC8E", "LRRN4", "LTB", "LYST", "MAF", "MAP3K5", "MCOLN3", "MICAL1",
    "MITD1", "MLLT6", "MMAA", "MPO", "MPV17L", "MS4A6A", "MSRB2", "MTUS1",
    "MYBPC2", "MYCBP2", "NCEH1", "NCF2", "NFU1", "NIT2", "NLRC5", "NQO2",
    "NT5C3A", "OAS3", "OCIAD1", "PARP14", "PCNX2", "PGLYRP2", "PHF23", "PHLPP1",
    "PIAS1", "GSAP", "PLEK", "PLSCR3", "HELZ2", "PRTN3", "PSME2", "PTPRZ1",
    "RAB5C", "IFT22", "RASGEF1B", "RGS17", "RNASEL", "RNF145", "RNF157",
    "RNF31", "RRAD", "SAMHD1", "SCLY", "SENP8", "SERPINB1", "SERPINB6",
    "SERPINB9", "SERPINE2", "SERPINI1", "SESN1", "SETDB1", "SGK3", "SLC15A3",
    "SLC6A12", "SLITRK6", "SOCS1", "SPI1", "DENND2B", "STAT2", "STIM2",
    "TAPBPL", "TBC1D21", "TBC1D22B", "TBK1", "TDRD7", "TFEC", "TLR1", "TLR4",
    "TLR7", "TMED6", "TMEM140", "TMEM17", "TNNI2", "TOR3A", "TRIM14", "TRIM21",
    "TRIM26", "TRIM5", "TSPAN14", "TTYH3", "USP43", "USP8", "VPS54", "WDR37",
    "ZFYVE26", "ZNF513", "ZSWIM2"
  )
  exhaustion_related_MemVE_D<-c("ABCB11", "ACADSB", "ACOX1", "ACTR8", "ADD1", "AGFG2", "ALDOC", "ANXA5", "APLP2", "ARFGAP3",
                                "ARFIP2", "ARMCX2", "ATF6", "ATP8B2", "BDP1", "BPGM", "MPC2", "FUOM", "MYDGF", "CIMIP5",
                                "YAE1", "CA2", "CABP4", "CACNG8", "CALU", "CAPN2", "CASP3", "CCDC141", "CCDC88A", "CD2BP2",
                                "CD48", "CD81", "CDC42SE2", "CDK2AP2", "CFI", "CGRRF1", "CKAP4", "CLINT1", "CLPTM1L", "CMAS",
                                "CPD", "CPNE3", "CRADD", "CREB3L2", "CREG1", "CRIP1", "DCAF12", "DCAF6", "DCUN1D4", "INTS6L",
                                "DGKA", "DHDDS", "DKKL1", "DNASE2", "DNPEP", "DPY19L3", "DRAM2", "RIPPLY3", "DYRK3", "ECHDC1",
                                "ECI2", "EDEM1", "EIF2AK3", "EPCAM", "ERLEC1", "ERP44", "ETHE1", "ETV5", "PABIR1", "MVB12A",
                                "FBXO2", "FER", "FKBP11", "GALNT2", "GKAP1", "GNAQ", "GNPNAT1", "GRIN2A", "GSTT1", "HADHB",
                                "HEXIM2", "HM13", "HMOX2", "HNRNPLL", "ICA1", "IDH2", "IFT80", "IGF1R", "IKZF3", "IL10RB",
                                "IL13RA1", "IL18", "IL6R", "INS", "IRF4", "ITFG1", "JTB", "KDELR1", "GARRE1", "KIFC2",
                                "KRTCAP2", "LHPP", "LPCAT3", "MAGED1", "MAN2A1", "MBTPS1", "MCCC2", "MCFD2", "MEGF8", "MGAT1",
                                "MRPL57", "MRPS14", "MTOR", "MYADM", "MYLK", "NARF", "NDFIP1", "NDUFA1", "NLRP6", "OAT",
                                "OSBP", "P2RX4", "PANK3", "PBK", "PDIA4", "PDIA6", "PEX2", "PGK1", "PGRMC1", "JADE1",
                                "PHYHD1", "PON3", "PPP1R15A", "PPP1R15B", "PRDX2", "PRDX4", "PRR3", "PRRC2A", "HACD2",
                                "RAB3D", "RAI1", "RAP1A", "RASL12", "RDH10", "RELN", "REXO2", "RTN4RL1", "S100PBP", "SAA4",
                                "SBF2", "SEC23B", "SELENOF", "SEPSECS", "SKIL", "SLC44A1", "SMAD4", "SMAD5", "SNX21", "SOX13",
                                "SPCS1", "SPCS3", "SPINT2", "SRSF9", "SSR4", "ST14", "STAT5B", "STX5", "SYNRG", "TEX14",
                                "TM2D1", "TM2D3", "TMBIM4", "TMED10", "TMEM150A", "TMEM165", "TMEM167A", "TMEM182", "TMEM25",
                                "TMEM52", "TNFRSF17", "TNS3", "TRABD", "TRAM1", "TRAPPC2B", "TSPAN13", "TXLNG", "TXNDC5",
                                "UAP1", "UBE2A", "UBN1", "UGDH", "USP22", "VCL", "TCEAL9", "WDR27", "ZCWPW2", "ZKSCAN1",
                                "ZNF148", "ZNF496", "ZSCAN21")
  ex_up<-intersect(exhaustion_related_NVE_D, exhaustion_related_MemVE_D)
  ex_down<-intersect(exhaustion_related_MemVE_U, exhaustion_related_NVE_U)
  ex_up<-intersect(exhaustion_related_NVE_D, exhaustion_related_MemVE_D)
  ex_down<-intersect(exhaustion_related_MemVE_U, exhaustion_related_NVE_U)
  CD8<-subset(all, subset= singler_labels2=="CD8+ T-cells")
  CD8 <- NormalizeData(CD8)
  CD8 <- FindVariableFeatures(CD8, selection.method = "vst", nfeatures = 2000)
  CD8<- ScaleData(CD8)
  CD8<- RunPCA(CD8, features = VariableFeatures(object = CD8))
  table(CD8$patient_id)
  pre_cells <- grep("Pre", CD8@meta.data$patient_id, ignore.case = TRUE, value = TRUE)
  keep_cells <- !grepl("Pre", CD8@meta.data$patient_id, ignore.case = TRUE)
  CD8post <- CD8[, keep_cells]
  table(CD8post$patient_id)
  saveRDS(CD8,file = "CD8.rds")
  library(VennDiagram)
  library(grid)
  venn.diagram(
    x = list(V1 = exhaustion_related_NVE_D, V2 = exhaustion_related_MemVE_D),
    filename = "ex_up.png",
    fill = c("#87CEEB", "#FA8072"),
    alpha = 0.7,
    cat.col = c("#87CEEB", "#FA8072"),
    cat.cex = 1.2,
    margin = 0.1,
    cex = 1.5,
    category.names = c("Exhausted vs Naive UP", "Exhausted vs Memory UP"),
    output = TRUE
  )
  venn.diagram(
    x = list(V1 = exhaustion_related_NVE_U, V2 = exhaustion_related_MemVE_U),
    filename = "ex_down.png",
    fill = c("#C13B2A", "#E3B130"),
    alpha = 0.7,
    cat.col = c("#C13B2A", "#E3B130"),
    cat.cex = 1.2,
    margin = 0.1,
    cex = 1.5,
    category.names = c("Exhausted vs Naive DOWN", "Exhausted vs Memory DOWN"),
    output = TRUE
  )
  glycolysis_genes <- c(
    "HK1", "HK2", "HK3", "GCK",
    "GPI",
    "PFKL", "PFKM", "PFKP",
    "ALDOA", "ALDOB", "ALDOC",
    "TPI1",
    "GAPDH",
    "PGK1",
    "PGAM1", "PGAM2",
    "ENO1", "ENO2", "ENO3",
    "PKLR", "PKM",
    "LDHA", "LDHB", "LDHC", "LDHD"
  )
  exhaustion_related_genes<-c("PDCD1")
  library(tibble)
  avg_expr <- FetchData(
    CD8,
    vars = c(glycolysis_genes,exhaustion_related_genes,"response", "patient_id"),
    slot = "data"
  ) %>%
    pivot_longer(
      cols = all_of(glycolysis_genes),
      names_to = "Gene",
      values_to = "Expression"
    ) %>%
    group_by(patient_id, response, Gene) %>%
    summarise(
      Expression = mean(Expression, na.rm = TRUE),
      .groups = "drop"
    )
  avg_expr$response <- factor(
    avg_expr$response,
    levels = c("R", "NR")
  )
  ggplot(avg_expr, aes(x = response, y = Expression, fill = response)) +
    geom_violin(
      width = 0.8,
      position = position_dodge(0.8),
      alpha = 0.3,
      trim = TRUE,
      scale = "width",
    )+
    geom_boxplot(
      width = 0.15,
      position = position_dodge(0.8),
      outlier.shape = NA,
      alpha = 0.8,
      color = "black"
    ) +
    geom_point(
      aes(group = response),
      position = position_jitterdodge(
        jitter.width = 0.15,
        dodge.width = 0.8
      ),
      shape = 21,
      size = 2.5,
      color = "black",
      alpha = 0.7,
      show.legend = FALSE
    ) +
    facet_wrap(
      ~ Gene,
      ncol = 5,
      scales = "free_y",
      labeller = labeller(Gene = label_value)
    ) +
    stat_compare_means(
      aes(group = response),
      method = "wilcox.test", paired = FALSE,
      method.args = list(alternative = "less"),
      size = 3,
    ) +
    scale_fill_manual(values = c("R" = "#E64B35", "NR" = "#4DBBD5")) +
    labs(
      x = "Treatment Response",
      y = "Normalized Expression",
      title = "Glycolysis Pathway Enzyme Expression by Response Group"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10, margin = margin(3,0,3,0))
    )
  ggsave("glycolysis_expression in CD8+T.pdf", width = 10, height = 12.5, dpi = 300)
  pgk1_data_single <- avg_expr %>%
    filter(Gene == "PGK1")
  ggplot(pgk1_data_single, aes(x = response, y = Expression, fill = response)) +
    geom_violin(
      width = 0.8,
      position = position_dodge(0.8),
      alpha = 0.3,
      trim = TRUE,
      scale = "width",
    )+
    geom_boxplot(
      width = 0.2,
      position = position_dodge(0.8),
      outlier.shape = NA,
      alpha = 0.8,
      color = "black"
    ) +
    geom_point(
      aes(group = response),
      position = position_jitterdodge(
        jitter.width = 0.15,
        dodge.width = 0.8
      ),
      shape = 21,
      size = 2.5,
      color = "black",
      alpha = 0.7,
      show.legend = FALSE
    ) +scale_y_continuous(limits = c(0, max(pgk1_data_single$Expression) * 1.2), expand = c(0, 0))+
    stat_compare_means(
      aes(group = response),
      method = "wilcox.test", paired = FALSE,
      method.args = list(alternative = "greater"),
      size = 6,
      label.y = max(pgk1_data_single$Expression) * 1.1
    ) +
    scale_fill_manual(values = c("R" = "#E64B35", "NR" = "#4DBBD5")) +
    labs(
      x = "Treatment Response",
      y = "PGK1 Expression Level",
      title = "PGK1 Expression by Response Group"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, size = 12, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      legend.position = "none"
    )
  ggsave("PGK1_expression in CD8+T.pdf", width = 4, height = 5, dpi = 300)
  avg_expr_ex <- FetchData(
    CD8,
    vars = c(exhaustion_related_genes,"response", "patient_id"),
    slot = "data"
  ) %>%
    pivot_longer(
      cols = all_of(exhaustion_related_genes),
      names_to = "Gene",
      values_to = "Expression"
    ) %>%
    group_by(patient_id, response, Gene) %>%
    summarise(
      Expression = mean(Expression, na.rm = TRUE),
      .groups = "drop"
    )
  avg_expr_ex$response <- factor(
    avg_expr_ex$response,
    levels = c("R", "NR")
  )
  pd1_data_single <- avg_expr_ex %>%
    filter(Gene == "PDCD1")
  ggplot(pd1_data_single, aes(x = response, y = Expression, fill = response)) +
    geom_violin(
      width = 0.8,
      position = position_dodge(0.8),
      alpha = 0.3,
      trim = TRUE,
      scale = "width",
    )+
    geom_boxplot(
      width = 0.2,
      position = position_dodge(0.8),
      outlier.shape = NA,
      alpha = 0.8,
      color = "black"
    ) +
    geom_point(
      aes(group = response),
      position = position_jitterdodge(
        jitter.width = 0.15,
        dodge.width = 0.8
      ),
      shape = 21,
      size = 2.5,
      color = "black",
      alpha = 0.7,
      show.legend = FALSE
    ) +scale_y_continuous(limits = c(0, max(pd1_data_single$Expression) * 1.2), expand = c(0, 0))+
    stat_compare_means(
      aes(group = response),
      method = "wilcox.test",
      size = 6,
      label.y = max(pd1_data_single$Expression) * 1.1
    ) +
    scale_fill_manual(values = c("R" = "#E64B35", "NR" = "#4DBBD5")) +
    labs(
      x = "Treatment Response",
      y = "PDCD1 Expression Level",
      title = "PDCD1 Expression by Response Group"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, size = 12, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      legend.position = "none"
    )
  ggsave("Pd1_expression in CD8+T.pdf", width = 4, height = 5, dpi = 300)
  exhaustion_set <- list(Exhaustion = c("CXCL13","HAVCR2","PDCD1","TIGIT","LAG3","CTLA4","LAYN","RBPJ","VCAM1","GZMB","TOX","MYO7A"))
  CD8 <- AddModuleScore(CD8, features = exhaustion_set, ctrl = 100, name = "Exhaustion_")
  CD8$Exhaustion_Score <- CD8@meta.data$Exhaustion_1
  exhaustion_plot_data <- data.frame(
    Exhaustion_Score = CD8$Exhaustion_Score,
    response = CD8$response
  )
  exhaustion_plot_data$response <- factor(
    exhaustion_plot_data$response,
    levels = c("R", "NR")
  )
  ggplot(exhaustion_plot_data, aes(x = response, y = Exhaustion_Score, fill = response)) +
    geom_violin(
      width = 0.8,
      position = position_dodge(0.8),
      alpha = 0.3,
      trim = TRUE,
      scale = "width",
    ) +
    geom_boxplot(
      width = 0.2,
      position = position_dodge(0.8),
      outlier.shape = NA,
      alpha = 0.8,
      color = "black"
    ) +
    scale_y_continuous(
      limits = c(min(exhaustion_plot_data$Exhaustion_Score)*0.9,
                 max(exhaustion_plot_data$Exhaustion_Score)*1.3),
      expand = c(0, 0)
    ) +
    stat_compare_means(
      aes(group = response),
      method = "wilcox.test",
      size = 6,
      vjust = -0.5,
      label.y = max(exhaustion_plot_data$Exhaustion_Score)*1.2
    ) +
    scale_fill_manual(values = c("R" = "#E64B35", "NR" = "#4DBBD5")) +
    labs(
      x = "Treatment Response",
      y = "T Cell Exhaustion Score",
      title = "T Cell Exhaustion by Response Group"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, size = 12, face = "bold", hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      legend.position = "none"
    )
  ggsave("exhaustion score in CD8+T.pdf", width = 4, height = 5, dpi = 300)
  expr_data <- FetchData(CD8, vars =c("PGK1","PDCD1"), slot = "data")
  head(expr_data)
  expr_data$cell <- rownames(expr_data)
  pd1_positive <- expr_data %>% filter(PDCD1 > 0)
  pgk1_nonzero <- pd1_positive %>% filter(PGK1 > 0)
  pgk1_cutoff <- quantile(pgk1_nonzero$PGK1, probs = c(0.3, 0.7))
  expr_data <- expr_data %>% mutate(
    group = case_when(
      PDCD1 == 0  ~ "PD1neg",
      PDCD1 > 0 & PGK1 > 0 & PGK1 <= pgk1_cutoff[1] ~ "PD1pos_PGK1low",
      PDCD1 > 0 & PGK1 > 0 & PGK1 >= pgk1_cutoff[2] ~ "PD1pos_PGK1high",
      TRUE ~ "Other"
    )
  )
  CD8$PD1_pgk1 <- expr_data$group
  keep<-c("PD1neg","PD1pos_PGK1low","PD1pos_PGK1high")
  CD8_PD1_PGK1<-subset(CD8, subset= PD1_pgk1%in%keep)
  exhaustion_plot_data <- data.frame(
    Exhaustion_Score = CD8_PD1_PGK1$Exhaustion_Score,
    PD1_pgk1 = CD8_PD1_PGK1$PD1_pgk1
  )
  exhaustion_plot_data$PD1_pgk1 <- factor(
    exhaustion_plot_data$PD1_pgk1,
    levels = c("PD1neg","PD1pos_PGK1low","PD1pos_PGK1high")
  )
  ggplot(exhaustion_plot_data, aes(x = PD1_pgk1, y = Exhaustion_Score, fill = PD1_pgk1)) +
    geom_boxplot(width = 0.15, position = position_dodge(0.8), outlier.shape = NA, alpha = 0.8, color = "black") +
    geom_violin(
      width = 0.8,
      position = position_dodge(0.8),
      alpha = 0.3,
      trim = TRUE,
      scale = "width",
    ) +
    geom_signif(
      comparisons = list(
        c("PD1neg", "PD1pos_PGK1low"),
        c("PD1neg", "PD1pos_PGK1high"),
        c("PD1pos_PGK1low", "PD1pos_PGK1high")
      ),
      test = "wilcox.test",
      label = "p.format",
      map_signif_level=FALSE,
      textsize = 5,
      tip_length = 0.01,
      vjust = 0.1,
      y_position = c(
        max(exhaustion_plot_data$Exhaustion_Score) * 1.0,
        max(exhaustion_plot_data$Exhaustion_Score) * 1.12,
        max(exhaustion_plot_data$Exhaustion_Score) * 0.9
      )
    ) +
    scale_y_continuous(
      limits = c(min(exhaustion_plot_data$Exhaustion_Score) * 1.2,
                 max(exhaustion_plot_data$Exhaustion_Score) * 1.25)
    ) +
    scale_fill_manual(values = c("PD1neg" = "white", "PD1pos_PGK1low" = "#4DBBD5","PD1pos_PGK1high" = "#E64B35")) +
    labs(
      x = "Treatment Response",
      y = "T Cell Exhaustion Score",
      title = "T Cell Exhaustion by PD1_pgk1 Group"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, size = 12, face = "bold",hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      legend.position = "none"
    )
  ggsave("exhaustion score in PGK1_CD8+T.pdf", width = 5, height = 5, dpi = 300)
  saveRDS(CD8, paste0(dataset_name, "_CD8_final.rds"))
  invisible(list(all = all, CD8 = CD8))
}
results <- lapply(names(seurat_list), function(ds) {
  message("Running ", ds, " (", cancer_names[[ds]], ")")
  tryCatch(run_analysis(seurat_list[[ds]], ds), error = function(e) {
    message(ds, " failed: ", conditionMessage(e))
    NULL
  })
})
names(results) <- names(seurat_list)
saveRDS(results, file.path(base_dir, "results", "four_cohort_analysis_results.rds"))
