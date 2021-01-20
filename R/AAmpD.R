AAmpD = function(path_to_files,t_cluster, tumor_prefix, normal_prefix, normal_clusters,
      mapp_file ,blacklist, genome_file ,bin_size,sd, cut_off) {
  print(paste0("Processing cluster ", t_cluster))
  mappability=read.table(mapp_file, stringsAsFactors = FALSE)
  black_list=read.table(blacklist, stringsAsFactors = FALSE)
  black_list.gr = GRanges(
    black_list[,1], 
    IRanges(black_list[,2], black_list[,3])
  );
  genome=read.table(genome_file, stringsAsFactors = FALSE)
  chr_size=genome[,2]
  names(chr_size)=as.character(genome[,1])
  bins.gr=tileGenome(seqlengths=chr_size, tilewidth=bin_size, cut.last.tile.in.chrom = TRUE)
  idy = queryHits(findOverlaps(bins.gr, black_list.gr))
  
  chrom=read.table(genome_file, stringsAsFactors=FALSE)
  chrom.gr=GRanges(chrom[,1], 
                   IRanges(rep(1,nrow(chrom)), chrom[,2]))
  idy1 = queryHits(findOverlaps(bins.gr, chrom.gr, type = "within"));
  idy_final=idy1[-which(idy1 %in% idy)]
  gc_features=GCcontent(Hsapiens,bins.gr[idy_final,])
  
  t_peaks_size=read.table(paste0(path_to_files, "/", tumor_prefix,"_peak_size.txt"), stringsAsFactors = FALSE)
  t_counts=read.table(paste0(path_to_files, "/",tumor_prefix, "_cluster", t_cluster , "_nopeaks_counts.txt"), stringsAsFactors = FALSE)
  t_seq_depth=countLines(paste0(path_to_files, "/",tumor_prefix, "_cluster", t_cluster , ".bed"))
  t_norm_counts = ((t_counts[idy_final,4]/(bin_size-t_peaks_size[idy_final,4]))/t_seq_depth)*10^10
  if(length(normal_clusters) > 1){
    n_norm_counts_all=NULL
    for(j in normal_clusters){
      print(paste0("normal", j))
      n_peaks_size=read.table(paste0(path_to_files, "/",normal_prefix, "_peak_size.txt"), stringsAsFactors = FALSE)
      n_counts=read.table(paste0(path_to_files, "/",normal_prefix, "_cluster", j , "_nopeaks_counts.txt"), stringsAsFactors = FALSE)
      n_seq_depth=countLines(paste0(path_to_files, "/",normal_prefix, "_cluster", j , ".bed"))
      n_norm_counts_all=cbind(n_norm_counts_all, ((n_counts[idy_final,4]/(bin_size-n_peaks_size[idy_final,4]))/n_seq_depth)*10^10)
    }
    n_norm_counts=rowMeans(n_norm_counts_all)
  } else {
    n_peaks_size=read.table(paste0(path_to_files, "/",normal_prefix, "_peak_size.txt"), stringsAsFactors = FALSE)
    n_counts=read.table(paste0(path_to_files, "/",normal_prefix, "_cluster", normal_clusters , "_nopeaks_counts.txt"), stringsAsFactors = FALSE)
    n_seq_depth=countLines(paste0(path_to_files, "/",normal_prefix, "_cluster", normal_clusters , ".bed"))
    n_norm_counts = ((n_counts[idy_final,4]/(bin_size-n_peaks_size[idy_final,4]))/n_seq_depth)*10^10
  }

  counts_gc=data.frame(chr=t_counts[idy_final,1], start=t_counts[idy_final,2], end=t_counts[idy_final,3],tumor=t_norm_counts, control=n_norm_counts,
                       gc=gc_features[,1], mapp=mappability[idy_final,1])
  counts_gc_2=counts_gc[which(counts_gc$gc > 0),]
  counts_gc_2=counts_gc_2[order(counts_gc_2$gc),]
  counts_gc_2$row_id=1:nrow(counts_gc_2)
  counts_gc_2$chr=factor(counts_gc_2$chr, levels=c(paste("chr", 1:22, sep = ""), "chrX", "chrY"))
  counts_gc_2=counts_gc_2[order(counts_gc_2$chr, counts_gc_2$start),]
  norm_counts_gc=apply(counts_gc_2,1, gc_log, n=4)
  control_counts_gc=apply(counts_gc_2,1, gc_log, n=5)
  norm_counts_normal_log=log(norm_counts_gc/control_counts_gc,2)
  norm_counts_normal_log[which(is.na(norm_counts_normal_log) | is.infinite(norm_counts_normal_log))]=0
  CNA.object <- CNA(cbind(norm_counts_normal_log),
                    counts_gc_2$chr,counts_gc_2$start,
                    data.type="logratio",sampleid=tumor_prefix)
  smoothed.CNA.object <- smooth.CNA(CNA.object)
  segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1,undo.splits="sdundo", undo.SD=sd)
  #segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
  
  des=describe(segment.smoothed.CNA.object$output$seg.mean)
  pdf(paste0(tumor_prefix,"_", normal_prefix, "_DNAcopy_cluster",t_cluster, ".pdf"), width=21, height=21)
  plot(segment.smoothed.CNA.object, plot.type="w", pt.cols = c("darkseagreen3", "darkslategray1"), zeroline=FALSE, main=samplename)
  abline(h=des$mean+(cut_off*des$sd), col="blue")
  dev.off()
  
  final_amplicons=subset(segment.smoothed.CNA.object$output, seg.mean > des$mean+(cut_off*des$sd))
  if(nrow(final_amplicons) > 0){
    final_amplicons$counts=0
    final_amplicons$mapp=0
    counts_gc.gr=GRanges(
      counts_gc_2[,1], 
      IRanges(counts_gc_2[,2], counts_gc_2[,3]))
    for(i in 1:nrow(final_amplicons)){
      final_amplicons.gr=GRanges(
        final_amplicons[i,2], 
        IRanges(as.numeric(final_amplicons[i,3]), as.numeric(final_amplicons[i,4])))
      idy = queryHits(findOverlaps(counts_gc.gr, final_amplicons.gr))
      final_amplicons[i,]$counts = mean(as.numeric(counts_gc_2[idy,]$tumor))
      final_amplicons[i,]$mapp = mean(as.numeric(counts_gc_2[idy,]$mapp))
    }
    final_amplicons_2=subset(final_amplicons, counts > quantile(counts_gc_2$tumor,0.99) & mapp > 0.8) 
  }
  return(data.frame(final_amplicons_2))
}

