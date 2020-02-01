args = commandArgs(trailingOnly=TRUE)
library(DESeq2)

data<-read.delim(args[1], header=T, row.names=1)
samples=ncol(data)/3-1
for( i in 1:samples){
  #if( i == samples) next
  for( j in i:samples+1){
      if( i == j) next
      s1 <- (i-1)*3+1
      e1 <- (i-1)*3+3
      n1  <- strsplit(colnames(data)[s1], ".", fixed=T)[[1]][2]
      s2 <- (j-1)*3+1
      e2 <- (j-1)*3+3
      n2<- strsplit(colnames(data)[s2], ".", fixed=T)[[1]][2]
      out_f1 <- paste(paste("DESeq2", n1, n2, sep="-"), "out", sep=".")

      y<-data[,c(s1:e1, s2:e2)]

      print(paste( s1, e1, n1, s2,e2, n2, sep= "_"))
      print("Starting DESeq")
      
      #前処理(DESeqDataSetオブジェクトの作成)
      data.cl <- c(rep(1, 3), rep(2, 3))#G1群を1、G2群を2としたベクトルdata.clを作成
      colData <- data.frame(condition=as.factor(data.cl))#condition列にクラスラベル情報を格納したcolDataオブジェクトを作成
      d <- DESeqDataSetFromMatrix(countData=y, colData=colData, design=~condition)#DESeqDataSetオブジェクトdの作成

      #本番(DEG検出)
      #d <- estimateSizeFactors(d)           #正規化を実行した結果をdに格納
      #d <- estimateDispersions(d)           #モデル構築(ばらつきの程度を見積もっている)
      #d <- nbinomLRT(d, full= ~condition, reduced= ~1)#検定
      d <- DESeq(d)                          #DESeq2を実行
      tmp <- results(d)                      #実行結果を抽出
      p.value <- tmp$pvalue                  ##p-valueをp.valueに格納
      p.value[is.na(p.value)] <- 1           #NAを1に置換している
      q.value <- tmp$padj                    #adjusted p-valueをq.valueに格納
      q.value[is.na(q.value)] <- 1           #NAを1に置換している
      ranking <- rank(p.value)               #p.valueでランキングした結果をrankingに格納

      #ファイルに保存(テキストファイル)
      tmp <- cbind( n1, n2, rownames(data), y, p.value, q.value, ranking)#入力データの右側にDEG検出結果を結合したものをtmpに格納
      write.table(tmp, out_f1, sep="\t", append=F, quote=F, row.names=F)#tmpの中身を指定したファイル名で保存
  }
}
