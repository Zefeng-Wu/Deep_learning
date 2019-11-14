## 1. keras mnist data test
library(keras)
### plot confuse matrix
Plot_confuse_matrix<-function(model,x_test,y_test){
  require(tidyverse)
  perf = model %>% evaluate(x_test, y_test)
  acc     = perf$acc %>% round(3)*100
  y_pred  = model %>% predict_classes(x_test)
  y_real  = y_test %>% apply(1,function(x){return(which(x==1) - 1) })
  results = tibble(y_real = y_real , y_pred = y_pred ,
                   Correct = ifelse(y_real == y_pred,"yes","no"))
  title = 'Performance on 40% unseen data - Feed Forward Neural Network'
  xlab  = 'True class'
  ylab  = 'Predicted class'
  results %>%
    ggplot(aes(x = factor(y_real), y = factor(y_pred), colour = Correct)) +
    geom_point() +
    ggtitle(label = title, subtitle = paste0("Accuracy = ", acc,"%")) +
    xlab(xlab) +
    ylab(ylab) +
    scale_color_manual(labels = c('No', 'Yes'),
                       values = c('tomato','cornflowerblue')) +
    geom_jitter() +
    theme_bw()
}
Plot_confuse_matrix(model = network,x_test = test_images,y_test = test_labels)

##### 2. sequence classification based on keras
#### sequence classification

# extract TSS sequences
GetPromoterSeqs<-function(gtf_file="~/MyResearch/genome.db/TAIR/gtf/Arabidopsis_thaliana.TAIR10.31.protein_coding.gtf",
                          genome.fa_file="~/MyResearch/genome.db/TAIR/dna/Arabidopsis_thaliana.TAIR10.31.dna.toplevel.fa",
                          up_stream_length = 400,
                          down_stream_length = 200){
  library(Biostrings)
  library(GenomicFeatures)
  library(BSgenome)
  tx<-makeTxDbFromGFF(gtf_file)
  
  ara_dna_seqs<-readDNAStringSet(genome.fa_file)
  genes<-genes(tx)
  seqlengths(genes)<-width(ara_dna_seqs)[match(names(seqlengths(genes)),names(ara_dna_seqs))] 
  
  ## promoter ranges
  AG_promoters<-promoters(genes,upstream = up_stream_length,downstream = down_stream_length)
  AG_promoters<-trim(AG_promoters)
  AG_promoters<-restrict(AG_promoters,end=seqlengths(AG_promoters))
  
  ## promtoer sequences
  AG_promoter_seq<- getSeq(ara_dna_seqs, AG_promoters)
  names(AG_promoter_seq)<-AG_promoters$gene_id
  sequence_length = up_stream_length + down_stream_length 
  AG_promoter_seq<-AG_promoter_seq[width(AG_promoter_seq)==sequence_length]
  return(AG_promoter_seq)
}
AG_promoter_seq<-GetPromoterSeqs()
seq_length = width(AG_promoter_seq)[1]

## build tensor arrays
Seqs2matrix = function(seqs){
  return(do.call(what = rbind, args = strsplit(x = seqs, split = '')))
}
Tensor_array_build<-function(promoters_seqs){
  sequence_width<-width(promoters_seqs[1])
  sequence_num <-length(promoters_seqs)
  o_tensors <- array(rep(0,sequence_num*sequence_width*4), dim = c(sequence_num,sequence_width,4))
  rownames(o_tensors)<-names(AG_promoter_seq) ## add rownames
  
  a<-as.data.frame(promoters_seqs)
  b<-Seqs2matrix(seqs = a$x)
  rownames(b)<-rownames(a)
  
  for (i in 1:sequence_num){
    dna_i_sequence <- b[i,]
    message(i)
    temp_matrix<-matrix(0,sequence_width,4)
    colnames(temp_matrix)<-c("A","T","C","G")
    temp_matrix[,"A"]<-"A"
    temp_matrix[,"T"]<-"T"
    temp_matrix[,"C"]<-"C"
    temp_matrix[,"G"]<-"G"
    for (n in colnames(temp_matrix)){
      temp_matrix[,n]<-ifelse(temp_matrix[,n]==dna_i_sequence,1,0)
    }
    temp_matrix<-t(apply(temp_matrix,1,as.numeric))
    o_tensors[i,,] <-temp_matrix 
  }
  return(o_tensors)
}
o_tensors<-Tensor_array_build(promoters_seqs = AG_promoter_seq)

###### import imprinted genes ot other calsses of genes
df<-read.table("TPM_14seedling.txt",header = TRUE,row.names = 1)
vv<-apply(df,1,var)
library(OneR)
g<-as.character(bin(vv,nbins = 5,labels = c("0","1","2","3","4"),method = "content"))
names(g)<-rownames(df)
AG_promoter_seq<-AG_promoter_seq[names(AG_promoter_seq)%in%rownames(df)]


#### define training and testing data
train_test_class <- function(genes_list,ratio=0.9){
  train_sliced <-sample(seq(1:length(genes_list)),
                       ratio*length(genes_list),
                       replace = FALSE)
  train_genes<-genes_list[train_sliced]
  test_genes<-genes_list[-train_sliced]
  return(list(train_genes=train_genes,test_genes=test_genes))
}  
train_genes <-train_test_class(names(AG_promoter_seq),ratio = 0.6)$train_genes
test_genes <-train_test_class(names(AG_promoter_seq),ratio = 0.6)$test_genes

##### reshape array to long vector
x_train <- array_reshape(o_tensors[train_genes,,], c(length(train_genes), seq_length * 4))  ## extract arrary by rownames for test genes
x_test  <- array_reshape(o_tensors[test_genes,,],  c(length(test_genes), seq_length * 4)) 

##### label the train data and test data
y_train <- to_categorical(g[train_genes], num_classes = 5)  ## label pos and neg class for train data
y_test  <- to_categorical(g[test_genes],  num_classes = 5) # must input class variable from 0

#### build model
model <- keras_model_sequential() %>% 
  layer_dense(units  = seq_length * 4, activation = 'relu', input_shape = seq_length * 4) %>% 
  layer_dropout(rate = 0.4) %>% 
  layer_dense(units  = 100, activation  = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units  = 5, activation   = 'softmax')

## compile model
model %>% compile(
  loss      = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics   = c('accuracy')
)

## model training and testing
history = model %>% fit(
  x_train, y_train, 
  epochs = 150, 
  batch_size = 50, 
  validation_split = 0.2
)

### plot classifications results

Plot_confuse_matrix(model = model,x_test = x_test ,y_test = y_test )


### 2.CNN network
### CNN network

x_train <- array_reshape(o_tensors[train_genes,,], c(length(train_genes),seq_length,4,1))
x_test  <- array_reshape(o_tensors[test_genes,,],  c(length(test_genes),seq_length,4,1))

y_train <- to_categorical(g[train_genes], num_classes = 5)  ## label pos and neg class for train data
y_test  <- to_categorical(g[test_genes],  num_classes = 5) # must input class variable from 0

model <- keras_model_sequential() %>%
  layer_conv_2d(filters = 32, kernel_size = c(3,3), activation = 'relu',
                input_shape = c(seq_length,4,1)) %>%
  layer_dropout(rate = 0.25) %>% 
  layer_flatten() %>% 
  layer_dense(units  = seq_length*4, activation = 'relu') %>% 
  layer_dropout(rate = 0.4) %>% 
  layer_dense(units  = 100, activation  = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units  = 5, activation   = 'softmax')

model %>% compile(
  loss      = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics   = c('accuracy')
)

history = model %>% fit(
  x_train, y_train, 
  epochs = 150, 
  batch_size = 50, 
  validation_split = 0.2
)
Plot_confuse_matrix(model = model,x_test = x_test ,y_test = y_test)

### 3.regression analysis
### regression analysis

y_train <- vv[train_genes]  ## label pos and neg class for train data
y_test  <- vv[test_genes]

## model build
build_regression_model <- function() {
    model <- keras_model_sequential() %>%
    layer_dense(units = 64, activation = "relu",
                input_shape = dim(x_train)[2]) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 1)
  
    model %>% compile(
    loss = "mse",
    optimizer = optimizer_rmsprop(),
    metrics = list("mean_absolute_error")
  )
  model
}

model <- build_regression_model()
model %>% summary()

### training
history <- model %>% fit(
  x_train,
  y_train,
  epochs = 150,
  validation_split = 0.2,
  verbose = 0
)
### predicting
test_predictions <- model %>% predict(x_test)

