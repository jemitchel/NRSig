library(annotate)
library(jsonlite)
library("hgu133plus2.db")
library(ggplot2)

# Loads gene symbol annotations for all Affy probes
x <- hgu133plus2SYMBOL
mapped_probes <- mappedkeys(x)
glist <- as.list(x[mapped_probes]) # list of probes mapped to genes

# opens the tf-target probe data source
all_targets <- fromJSON("./data/NR_target_probes.json", flatten = TRUE)

# loads dataframe of serum starved MCF7 samples
data <- readRDS(file = "./data/serum-starved.rds")

# computes mean and standard deviation for serum starved data
data_av_std <- data.frame(
  Mean = rowMeans(data, na.rm = TRUE),
  Std = apply(data, 1, sd, na.rm = TRUE)
)
av <- data_av_std[, 1]
std <- data_av_std[, 2]

# all 48 NRs (some not tested since have < 15 target genes in TTRUST)
NRs <- c(
  "ESR1", "AR", "NR3C1", "NR3C2", "PGR",
  "NR4A1", "NR4A2", "NR4A3", "NR5A1", "NR5A2", "NR6A1",
  "NR0B1", "NR0B2", "THRA", "THRB", "RARA", "RARB", "RARG",
  "PPARA", "PPARD", "PPARG", "NR1D1", "NR1D2", "RORA", "RORB",
  "RORC", "NR1H4", "NR1H3", "NR1H2", "VDR", "NR1I2",
  "NR1I3", "HNF4A", "HNF4G", "RXRA", "RXRB", "RXRG", "NR2C1",
  "NR2C2", "NR2E1", "NR2E3", "NR2F1", "NR2F2", "NR2F6",
  "ESR2", "ESRRA", "ESRRB", "ESRRG"
)


#' DiffProbes2Genes
#'
#' Converts z scores from probe level -> gene level
#'
#' @param z_probe a dataframe of z scores for probes of either the current NR
#'  target or all other probes
#' @return 2 element list:
#' (1) list of max z scores for all represented genes and 
#' (2) list of probe that gave max z score of each gene
#'
DiffProbes2Genes <- function(z_probe) {
  z_genes <- list()
  mappings <- list()
  for (i in 1:nrow(z_probe)) {
    id <- rownames(z_probe)[i]
    if (id %in% mapped_probes) { # if the probe maps to a gene symbol, continue
      sym <- glist[[id]]
      score <- z_probe[i, 1]
      if (sym %in% names(z_genes)) {
        z_genes[[sym]] <- c(z_genes[[sym]], score)
        mappings[[sym]] <- c(mappings[[sym]], id)
      } else {
        z_genes[[sym]] <- c(score)
        mappings[[sym]] <- c(id)
      }
    }
  }

  newz <- list()
  newz_mappings <- list()
  for (i in 1:length(z_genes)) {
    ndx_max <- which.max(abs(z_genes[[i]]))
    newz[[names(z_genes)[i]]] <- z_genes[[i]][ndx_max]
    newz_mappings[[names(mappings)[i]]] <- mappings[[i]][ndx_max]
  }

  return(list(newz, newz_mappings))
}

#' GetNumAboveThreshold
#'
#' This function counts number of genes that are above a z-score threshold of 2
#'
#' @param zscores list of gene level z scores for either the itaration's
#'  NR targets or background
#' @return the number of genes with absolute value z score larger than 2
GetNumAboveThreshold <- function(zscores) {
  count <- 0
  thresh <- 2 # this is the zscore threshold
  for (i in 1:length(zscores)) {
    if (abs(zscores[[i]]) > thresh) {
      count <- count + 1
    }
  }
  return(count)
}

#' MakeBoxPlots
#'
#' Produces a figure of boxplots for distributions for all NR target
#' serum starved expression levels as well as the input's gene levels
#'
#' @param test_data dataframe of replicate averaged preprocessed input
#' @param gene_z list of max z scores for each target gene
#' @param mappings list of Affy probes corresponding to each target gene
#' @param NR character of the current NR being analyzed
#' @return a ggplot object of the figure
#' 
MakeBoxPlots <- function(test_data, gene_z, mappings, NR) {

  # preps serum-starved data for plotting
  distrib <- data[unlist(mappings), ]
  genelst <- rep(glist[rownames(distrib)], each = ncol(distrib))
  distrib <- c(t(distrib)) # linearizes data
  distrib <- cbind(distrib, genelst)
  distrib <- data.frame(distrib)
  colnames(distrib) <- c("expr", "gene")
  distrib <- as.data.frame(lapply(distrib, unlist))

  # preps user input data and results for plotting
  ordered_z <- as.list(gene_z[order(abs(unlist(gene_z)), decreasing = FALSE)])
  in_data <- test_data[unlist(mappings), , drop = FALSE]
  genelst <- glist[rownames(in_data)]
  in_data <- cbind(in_data, unlist(genelst))
  colnames(in_data) <- c("expr", "gene")
  in_data <- cbind(in_data, unlist(gene_z[as.character(in_data$gene)]))
  colnames(in_data)[3] <- "zsc"
  in_data <- in_data[match(names(ordered_z), in_data$gene), ]
  in_data$gene <- factor(in_data$gene, levels = in_data$gene) # orders the factor
  in_data$zsc <- round(in_data$zsc, digits = 2)

  # makes data to build the legend
  bplot_legend <- data.frame(
    c(26.5, 27),
    c(
      as.character(in_data[nrow(in_data) - 1, 2]),
      as.character(in_data[nrow(in_data) - 1, 2])
    )
  )
  colnames(bplot_legend) <- c("vals", "loc")
  in_plot_legend <- data.frame(
    c(26.5),
    c(as.character(in_data[nrow(in_data) - 2, 2]))
  )
  colnames(in_plot_legend) <- c("vals", "loc")
  text_legend <- data.frame(
    c("Serum-Starved", "Input Sample"),
    c(
      as.character(in_data[nrow(in_data) - 1, 2]),
      as.character(in_data[nrow(in_data) - 2, 2])
    )
  )
  colnames(text_legend) <- c("vals", "loc")

  # makes the plot
  p <- ggplot(data = in_data, aes(x = as.factor(gene), y = expr)) +
    geom_point(color = "red", size = 3) +
    geom_text(
      data = in_data, aes(x = as.factor(gene), y = 20, label = zsc),
      vjust = 0
    ) +
    geom_boxplot(
      data = distrib, aes(x = as.factor(gene), y = expr),
      width = .3
    ) +
    geom_point(
      data = in_data, aes(x = as.factor(gene), y = expr),
      color = "red", size = 3
    ) +
    coord_flip() +
    geom_line(
      data = bplot_legend, aes(x = as.factor(loc), y = vals),
      color = "black"
    ) +
    geom_point(
      data = in_plot_legend, aes(x = as.factor(loc), y = vals),
      color = "red", size = 3
    ) +
    geom_text(data = text_legend, aes(
      x = as.factor(loc), y = 24,
      label = vals
    )) +
    scale_y_continuous(breaks = seq(0, 15, 3)) +
    labs(
      title = paste(NR, " Target Expression", sep = ""),
      subtitle = "Z-Score", x = "Target Gene", y = "fRMA Expression"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, vjust = 2, face = "bold"),
      plot.subtitle = element_text(
        hjust = .7, vjust = -1,
        face = "bold"
      )
    ) +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )

  return(p)
}

#' UpdateTargets
#'
#' Modifies the list of NR target gene probes if the user has opted to remove
#' certain probes from the analysis (for example if they identified 
#' cross-hybridized probes in a xenograft study)
#'
#' @param data the dataframe of serum starved MCF7 gene expression
#' @return an updated list of NR target gene probes
#' 
UpdateTargets <- function(data) {
  for (i in 1:length(NRs)) {
    NR <- NRs[i]
    rems <- c()
    if (NR %in% names(all_targets)) { # if it hasnt alreay been removed...
      tlist <- all_targets[[NR]]
      for (j in 1:length(tlist)) {
        probe <- tlist[j]
        if (!(probe %in% rownames(data))) {
          rems <- c(rems, probe)
        }
      }

      if (length(rems) > 0) {
        for (j in 1:length(rems)) {
          tlist <- tlist[tlist != rems[j]]
        }
      }

      asGenes <- unique(glist[tlist])
      if (length(asGenes) < 15) {
        all_targets[[NR]] <- NULL
      } else {
        all_targets[[NR]] <- tlist
      }
    }
  }
  return(all_targets)
}


#' CalcEnrich
#'
#' This is the main function for data processing. It calculates differentially
#' expressed probes and NR enrichment by hypergeometric test. Also corrects P
#' values for multiple comparisons.
#'
#' @param test_samples dataframe of preprocessed input data
#' @param cross_probes dataframe of probes to remove (NULL if none)
#' @return 4 element list:
#' (1) list of list of NR target expression plots with number of targets
#' (2) dataframe of p values for overrepresentation of NR target sets
#' (3) dataframe of differential expression z scores
#' (4) dataframe of the preprocessed input used for calculations
#' 
CalcEnrich <- function(test_samples, cross_probes) {

  # keeps track of progress
  num_completed <- 0

  # removes probes from ss data and target source if use xenograft data
  if (!is.null(cross_probes)) {
    hyb <- read.csv(cross_probes, stringsAsFactors = FALSE, header = FALSE)
    hyb <- as.list(hyb[, 1])
    data <- data[!(rownames(data) %in% hyb), ]
    data_av_std <- data_av_std[!(rownames(data_av_std) %in% hyb), ]
    av <- data_av_std[, 1]
    std <- data_av_std[, 2]
    all_targets <- UpdateTargets(data)
  }

  # averages replicates if there are multiple
  test_sample <- data.frame(Mean = rowMeans(test_samples, na.rm = TRUE))

  # compute zscores for all probes against the ss prior distributions
  zscores <- (test_sample - av) / std

  # creates df to store enrichment p-values
  final_p_vals <- data.frame(matrix(ncol = length(NRs), nrow = 1))
  colnames(final_p_vals) <- NRs

  NR_plots <- list() # holds z-scores for target genes of each NR

  for (d in 1:length(NRs)) { # for each nuclear receptor...

    # gets a list of target probes for the specified NR
    if (NRs[d] %in% names(all_targets)) {
      target_probes <- all_targets[[NRs[d]]]
    } else {
      next # if NR had less than 15 target genes go to next NR
    }

    # select just the probes that are targets of the current NR being evaluated
    targ_z_probe <- zscores[target_probes, , drop = FALSE]

    # select all other probes that are not listed as targets of the NR
    rest_z_probe <- zscores[!row.names(zscores) %in% target_probes, ,
      drop = FALSE
    ]

    # gets differentially expressed genes from probe z scores
    targ_z_gene <- DiffProbes2Genes(targ_z_probe)
    targ_z_gene_mappings <- targ_z_gene[[2]]
    targ_z_gene <- targ_z_gene[[1]]

    rest_z_gene <- DiffProbes2Genes(rest_z_probe)
    rest_z_gene_mappings <- rest_z_gene[[2]]
    rest_z_gene <- rest_z_gene[[1]]

    # generate boxplots of the data
    out_plot <- MakeBoxPlots(
      test_sample, targ_z_gene, targ_z_gene_mappings,
      NRs[d]
    )

    # figure out how many target and non-targets are dysregulated
    targ_above_thresh <- GetNumAboveThreshold(targ_z_gene)
    rest_above_thresh <- GetNumAboveThreshold(rest_z_gene)

    # design contingency table
    above_thresh <- c(targ_above_thresh, rest_above_thresh)
    below_thresh <- c(
      length(targ_z_gene) - targ_above_thresh,
      length(rest_z_gene) - rest_above_thresh
    )
    c_table <- cbind(above_thresh, below_thresh)

    # computes results
    result <- fisher.test(c_table, alternative = "greater")
    print(NRs[d])
    print(result[["p.value"]])
    NR_plots[[NRs[d]]] <- list(out_plot, length(targ_z_gene))
    final_p_vals[1, NRs[d]] <- result[["p.value"]]

    # calculates percentage complete
    num_completed <- num_completed + 1
    percent <- round(100 * num_completed / length(all_targets))

    incProgress(
      amount = 1 / length(all_targets),
      detail = paste(percent, "% complete", "")
    )
  }

  # removes NA values
  final_p_vals <- final_p_vals[, !apply(is.na(final_p_vals), 2, any)]

  # appends adjusted p-values
  padj <- p.adjust(final_p_vals[1, ], method = "fdr")
  final_p_vals <- rbind(final_p_vals, padj)

  # transposes and sorts dataframe
  final_p_vals <- t(final_p_vals)
  final_p_vals <- final_p_vals[order(final_p_vals[, 2]), ]
  NR_plots <- as.list(NR_plots[order(final_p_vals[, 2])]) # orders the list of plots
  # makes rownames into first column
  final_p_vals <- cbind(Row.Names = rownames(final_p_vals), final_p_vals)
  colnames(final_p_vals) <- c("NR", "raw p-value", "adj p-value")

  # adds gene symbol annotations to differentially expressed probes list
  syms <- as.character(glist[rownames(zscores)])
  syms <- replace(syms, syms == "NULL", "")
  zscores <- cbind(zscores, syms)
  # sorts differentially expressed probes list
  zscores <- zscores[order(abs(unlist(zscores[, 1])), decreasing = TRUE), ]

  # makes probe ID first column and names columns for both downloadable tables
  zscores <- cbind(Row.Names = rownames(zscores), zscores)
  colnames(zscores) <- c("Probeset ID", "Z-Score", "Gene Symbol")
  test_samples <- cbind(Row.Names = rownames(test_samples), test_samples)
  colnames(test_samples)[1] <- c("Probeset ID")

  return(list(NR_plots, final_p_vals, zscores, test_samples))
}
