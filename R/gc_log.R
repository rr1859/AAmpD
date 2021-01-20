gc_log=function(x,n){
  row_i=as.numeric(as.character(x[7]))
  if(row_i < 100){
    rows_k=row_i:(row_i+99)
  } else if(row_i > (nrow(counts_gc_2)-100)){
    rows_k=(row_i-99):row_i
  } else {
    rows_k=(row_i-50):(row_i+49)
  }
  return(as.numeric(as.character(x[n]))/mean(counts_gc_2[rows_k,n]))
}
