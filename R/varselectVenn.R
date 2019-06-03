#' Using Venn Diagrams to Compare High-Importance Variables Across Cancer Types
#'
#' This function generates a Venn diagram using RNA-seq data from the The Cancer Genome Atlas (TCGA) database.  Users specify which cancer types to include (varselectVenn currently supports 2- and 3-set Venn diagrams), as well as the target variable to predict.  The user-specified data is processed by a random forest classifier, and the variables (genes) are ranked by their influence on the model’s predictive power.  Users specify how many of the high-importance genes to retain, and a Venn diagram is generated that shows which genes are of high-importance among the different cancer types.  The function also returns a list object that specifies which genes were retained for each cancer type, as well as which genes were at the intersection of all specified cancers types. 
#'
#' @param types A vector of TCGA-supported acronyms that designate the type of cancer.  varselectVenn currently supports ACC, BLCA, KIRC, KIRP, LIHC, and THCA.  varselectVenn currently supports 2- and 3-set diagrams
#' 
#' @param num_var a numeric that specifies how many of the high-importance variables to retain.  The random forest classifier will rank each of the 20501 genes features in the order of their impact on the model.  Setting num_var to 100 for example, will retain the 100 most important variables for each cancer type to include in the Venn diagram.
#'
#' @param target the variable to be predicted.  varselectVenn currently supports tumor “patholigicstage” (which attempts to distinguish stage I tumors from stage II, III, and IV tumors), the patient’s “vitalststus” (a binary for whether the patient is alive or not), and the patient’s “gender”.
#'
#'
#' @details Cancer type acronyms:
#' ACC	  Adrenocortical Carcinoma,
#' BLCA	  Bladder Urothelial Carcinoma,
#' KIRC	  Kidney Renal Clear Cell Carcinoma,
#' KIRP	  Kidney Renal Papillary Cell Carcinoma,
#' LIHC	  Liver Heptocellular Carcinoma,
#' THCA	  Thyroid Carcinoma
#'
#'
#' @author Jacob Blamer, \email{jwilliamblamer@gmail.com}
#'
#'
#' @examples
#' varseleVenn(c("KIRP", "KIRC", "LIHC"), 100, "vitalstatus")
#' varseleVenn(c("BLCA", "THCA"), 75, "pathologicstage")
#'
#'
#' @import 
#' randomForest
#' e1071
#' TCGA2STAT
#' VennDiagram
#'
#' @export
varselectVenn <- function(types, num_var, target) {
  feat_selectOUT <- c()
  #### Import Data ####
  for (i in 1:length(types)) {
    varname <- paste("var_import", types[i], sep="")
    df1 <- getTCGA(types[i], "RNASeq2", "RPKM", clinical=T, cvars=target)
    df1.md <- df1$merged.dat
    
    if (target == "pathologicstage") {
      stage_i <- which(df1.md[,2] == "stage i")
      stage_i0 <- which(df1.md[,2] != "stage i")
      df1.md[stage_i,2] <- "si"
      df1.md[stage_i0,2] <- "si0"
      df1.md <- df1.md[-which(is.na(df1.md[,2])),]
    }
    
    print(paste("Random Forest Variable Selection for ", types[i]))
   
    ## Variable Selection w/RF
    rf <- randomForest(df1.md[3:20503],
                       as.factor(df1.md[,2]), 
                       ntree = 200)
    
    rf_nonzero_idxs <- which(rf$importance != 0)
    rf_var_sort <- head(sort(rf$importance[rf_nonzero_idxs,], decreasing = T), num_var)
    assign(varname, names(rf_var_sort))
    feat_selectOUT <- cbind(feat_selectOUT, names(rf_var_sort))
  }
  
  fin_df <- as.data.frame(feat_selectOUT)
   
  ## Venn Diagrams 
  if (length(types) == 2) {
    intsct <- intersect(fin_df$V1, fin_df$V2)
    venn_list <- list(fin_df$V1, fin_df$V2)
    names(venn_list) <- c(types[1], types[2])
    venn.plot <- venn.diagram(
      x = venn_list,
      euler.d = FALSE,
      scaled = FALSE,
      col = "transparent", 
      fill = c("cornflowerblue", "pink"),
      filename = "out_two.tiff",
      cat.pos = c(-20, 20),
      cat.dist = c(0.05, 0.05),
      cex = 2.5,
      cat.cex = 2.5,
    );
  }
    
    if (length(types) == 3) {
      intsct <- intersect(intersect(fin_df$V1, fin_df$V2), fin_df$V3)
      venn_list <- list(fin_df$V1, fin_df$V2, fin_df$V3)
      names(venn_list) <- c(types[1], types[2], types[3])
      venn.plot <- venn.diagram(
        x = venn_list,
        euler.d = FALSE,
        scaled = FALSE,
        col = "transparent", 
        fill = c("cornflowerblue", "yellow", "pink"),
        filename = "out_three.tiff",
        cat.pos = c(-20, 20, 180),
        cat.dist = c(0.05, 0.05, 0.05),
        cex = 2.5,
        cat.cex = 2.5,
      );
  }

  colnames(fin_df) <- types
  fin_list <- as.list(fin_df)
  fin_list$Intersect <- intsct
  
  return(fin_list)
  
}
  


  
