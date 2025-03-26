# Functional Enrichment Analysis of Conserved Overlaps
library(clusterProfiler)
library(org.Hv.eg.db) # Create custom annotation DB (see below)
library(tidyverse)
library(enrichplot)

# 1. Prepare Custom Annotation Database ----------------------------------------
# Create from your eggNOG results (run once per strain)
create_custom_db <- function(emapper_file) {
  annotations <- read_tsv(emapper_file, comment = "#", 
                         col_names = c("query","seed_ortholog","evalue","score",
                                      "eggNOG_OGs","max_annot_lvl","COG_category",
                                      "Description","Preferred_name","GOs","EC",
                                      "KEGG_ko","KEGG_Pathway","KEGG_Module",
                                      "KEGG_Reaction","KEGG_rclass","BRITE",
                                      "KEGG_TC","CAZy","BiGG_Reaction","PFAMs")) %>%
    mutate(
      GOs = str_replace_all(GOs, ",", ";"),
      KEGG_Pathway = str_replace_all(KEGG_Pathway, ",", ";")
    )
  
  # Create OrgDb-like object
  custom_db <- data.frame(
    GID = annotations$query,
    GO = annotations$GOs,
    KEGG = annotations$KEGG_Pathway,
    COG = annotations$COG_category,
    stringsAsFactors = FALSE
  ) %>%
    separate_rows(GO, sep = ";") %>%
    separate_rows(KEGG, sep = ";") %>%
    filter(!is.na(GO) | !is.na(KEGG))
  
  return(custom_db)
}

# 2. Load Overlap Genes --------------------------------------------------------
# Assuming you have a list of overlapping gene IDs per strain
overlap_genes <- list(
  strain1 = c("geneA", "geneB", ...),
  strain2 = c("geneC", "geneD", ...),
  ...
)

# 3. Enrichment Analysis Pipeline ---------------------------------------------
run_enrichment <- function(gene_list, custom_db, ontology="GO", p_cutoff=0.05) {
  
  # Create enrichResult object
  if(ontology == "GO") {
    enrich_res <- enricher(
      gene = gene_list,
      pvalueCutoff = p_cutoff,
      pAdjustMethod = "BH",
      universe = custom_db$GID,
      TERM2GENE = dplyr::select(custom_db, term = GO, gene = GID)
    )
  } else if(ontology == "KEGG") {
    enrich_res <- enricher(
      gene = gene_list,
      pvalueCutoff = p_cutoff,
      pAdjustMethod = "BH", 
      universe = custom_db$GID,
      TERM2GENE = dplyr::select(custom_db, term = KEGG, gene = GID)
    )
  }
  
  # Simplify redundant terms
  if(ontology == "GO") {
    enrich_res <- simplify(
      enrich_res, 
      cutoff = 0.7, 
      by = "p.adjust", 
      select_fun = min
    )
  }
  
  return(enrich_res)
}

# 4. Cross-Strain Meta-Analysis ------------------------------------------------
# Combine genes from all strains
all_overlap_genes <- unique(unlist(overlap_genes))

# Create combined background from all strains
combined_db <- map_dfr(list.files("eggnog_files", full.names = TRUE), create_custom_db)

# Run enrichment
go_res <- run_enrichment(all_overlap_genes, combined_db, "GO")
kegg_res <- run_enrichment(all_overlap_genes, combined_db, "KEGG")

# 5. Visualization -------------------------------------------------------------
# Dotplot of enriched terms
dotplot(go_res, showCategory=20, font.size=8) + 
  ggtitle("GO Enrichment - Conserved Overlap Genes")

# Network plot showing term relationships
cnetplot(go_res, foldChange=1, circular=FALSE, colorEdge=TRUE)

# KEGG Module visualization
emapplot(kegg_res) + 
  scale_color_gradient(low="blue", high="red") +
  theme_minimal()

# 6. Biological Interpretation -------------------------------------------------
# Extract top stress-related terms
stress_terms <- go_res@result %>%
  filter(grepl("response to stress|heat shock|oxidat", Description, ignore.case=TRUE)) %>%
  arrange(p.adjust)

# Print key findings
cat("Key Enriched Processes:\n")
cat(paste0("- ", stress_terms$Description[1:5], " (FDR=", 
          signif(stress_terms$p.adjust[1:5],2), ")\n"))
