list.of.packages <- c("ggplot2", "SHAPforxgboost", "xgboost", "fastshap", "farver")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(SHAPforxgboost)
library(xgboost) # for xgboost
library(ggplot2)
library(fastshap)


generate_SHAP_plot <- function(shap_df, descriptors_df){
  print("shap_df")
  print(shap_df)
  print("descriptors_df")
  print(descriptors_df)
  top.6.names <- names(sort(abs(shap_df[y, trimmed]), decreasing = T)) #get the names of the absolute 5 highest SHAP-values
  print(top.6.names[1:6])
  top.6.df <- as.data.frame(shap_df[y, top.6.names[1:6]])
  top.6.df.t <- as.data.frame(t(top.6.df)) #transpose the dataframe, as to make it easier to plot
  top.6.df.t$V2 <- t(descriptors_df[i, colnames(top.6.df)]) #add corresponding descriptor values
  sum <- rowSums(subset(shap_df[y,], select = trimmed))
  sum
  sum_label <- paste("SHAP value \n (Predicted ", ifelse(sum >= 0, "druggable)", "less druggable)"))
  colnames(top.6.df.t) <- c("SHAP values", "Descriptor values") #This is needed for y in ggplot
  top.6.barplot <- ggplot(data = top.6.df.t, aes(y = `SHAP values`, x = reorder(rownames(top.6.df.t), abs(`SHAP values`)), fill = abs(`SHAP values`))) + 
    geom_bar(stat = "identity") + 
    ylab(sum_label) +
    ggtitle(rownames(descriptors_df[i,])) +
    theme_classic() + 
    theme(legend.position = "none", axis.text.x = element_text(size= 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 17), axis.title.y = element_blank(), title = element_text(size = 17)) +
    geom_text(aes(y = 1.8, label = ifelse((top.6.df.t$Descriptor > 99), signif(top.6.df.t$Descriptor, 3), ifelse((top.6.df.t$Descriptor < 99), signif(top.6.df.t$Descriptor, 2), "")))) +
    coord_flip(clip = "off") +
    annotate("text", label = "Descriptor \n value", y = 1.7, x = 7.0) +
    ylim(-2.0, 2.0) +
    scale_color_hue(direction = 1)
  png(paste0("../ind_shap_plots/",rownames(descriptors_df[i,]),"_SHAP.png"), width = 1500, height = 1500, unit = "px", res = 300)
  print(top.6.barplot)
  dev.off()
}

##### Load & prepare
print("Load & Prepare")

model <- xgb.load('../scripts/DrugPred_RNA.model')

trimmed = c("psa_r", "fr_hpb_atoms", "fr_buried_sl_atoms", "hsa", "exp_sl_sa", "no_bs_atoms", "InertialShapeFactor", "no_sl_atoms", "sl_bs_r", "PMI3", "sa_vol_r", "SpherocityIndex")

descriptors <- read.csv("../descriptor_values.csv", row.names = 1, stringsAsFactors=FALSE) #do not convert strings to factor, because that how the code below is written

print(descriptors)

descriptors_trimmed <- descriptors[,trimmed]

##### Predict
print("Predict")

pred <- predict(model,as.matrix(descriptors_trimmed))
pred



##### Create out-file
print("Out-file")
scores <- cbind(rownames(descriptors_trimmed), pred) #combine RNA labels with their scores, so that it can be uploaded to MYSQL
write.csv(scores, '../Predictions.csv', row.names = F, quote = F)

##### Create SHAP-value plot
print("SHAP-values ")
all.pred.rows <- which(descriptors$csa>0)

shap_values <- shap.values(xgb_model = model, X_train = as.matrix(descriptors_trimmed))
shap_values_df <- as.data.frame(shap_values$shap_score)
shap_values_df


y <- 0
for (i in all.pred.rows){
  y <- y +1
  print(y)
  shap_df <- shap_values_df
  print("SHAP")
  print(shap_df)
  descriptors_df <- descriptors_trimmed
  print("DESCRIPTORS")
  print(descriptors_df)
  print("#")
  top.6.names <- names(sort(abs(shap_df[y, trimmed]), decreasing = T)) #get the names of the absolute 5 highest SHAP-values
  print("##")
  print(top.6.names[1:6])
  top.6.df <- as.data.frame(shap_df[y, top.6.names[1:6]])
  top.6.df.t <- as.data.frame(t(top.6.df)) #transpose the dataframe, as to make it easier to plot
  top.6.df.t$V2 <- t(descriptors_df[i, colnames(top.6.df)]) #add corresponding descriptor values
  #top.6.df.t.s <- as.data.frame(top.6.df.t[sort(abs(top.6.df.t$V1), decreasing = T, index.return =T)[[2]], ]) # sort the data frame columns, while retaining the sign
  sum <- rowSums(subset(shap_df[y,], select = trimmed))
  sum
  sum_label <- paste("SHAP value \n (Predicted ", ifelse(sum >= 0, "druggable)", "less druggable)"))
  colnames(top.6.df.t) <- c("SHAP values", "Descriptor values") #This is needed for y in ggplot
  top.6.barplot <- ggplot(data = top.6.df.t, aes(y = `SHAP values`, x = reorder(rownames(top.6.df.t), abs(`SHAP values`)), fill = abs(`SHAP values`))) + 
    geom_bar(stat = "identity") + 
    ylab(sum_label) +
    ggtitle(rownames(descriptors_df[i,])) +
    theme_classic() + 
    theme(legend.position = "none", axis.text.x = element_text(size= 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 17), axis.title.y = element_blank(), title = element_text(size = 17)) +
    geom_text(aes(y = 1.8, label = ifelse((top.6.df.t$Descriptor > 99), signif(top.6.df.t$Descriptor, 3), ifelse((top.6.df.t$Descriptor < 99), signif(top.6.df.t$Descriptor, 2), "")))) +
    coord_flip(clip = "off") +
    annotate("text", label = "Descriptor \n value", y = 1.7, x = 7.0) +
    ylim(-2.0, 2.0) +
    scale_color_hue(direction = 1)
  png(paste0("../ind_shap_plots/",rownames(descriptors_df[i,]),"_SHAP.png"), width = 1500, height = 1500, unit = "px", res = 300)
  print(top.6.barplot)
  dev.off()
}

