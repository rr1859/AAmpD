getBackgroundReads = function(snap_obj, clusters, snap_file, peak_file, bin_size, genome_file,sample_prefix, output_folder, path_to_bgreads){
  genome=read.table(genome_file, stringsAsFactors = FALSE)
  chr_size=genome[,2]
  names(chr_size)=as.character(genome[,1])
  tiles=tileGenome(seqlengths=chr_size, tilewidth=bin_size, cut.last.tile.in.chrom = TRUE)
  write.table(data.frame(tiles)[,1:3], paste0(genome_file, ".", bin_size, ".bed"), sep = "\t", quote = FALSE, col.names=FALSE, row.names=FALSE)
  dir.create(output_folder, showWarnings=FALSE)
  for(i in clusters){
    print(paste("Processing cluster", i))
    write.table(snap_obj@barcode[snap_obj@cluster == i],
                paste(sample_prefix, "_cluster",i, "_bc.txt",sep = ""),
                quote = FALSE, col.names=FALSE, row.names=FALSE)
    system2(command="/mnt/silencer2/home/yangli/apps/anaconda2/bin/snaptools",
            args=c("dump-fragment",
                   "--snap-file", snap_file,
                   "--output-file", paste0(output_folder, "/", sample_prefix, "_cluster",i, ".bed"),
                   "--barcode-file", paste0(sample_prefix, "_cluster",i, "_bc.txt")))
    
    system2(command=paste0("bash"),
            args=c("bg_reads_50kb.sh", paste0(sample_prefix, "_cluster",i), gsub(".bed", "", peak_file), paste0(genome_file, ".", bin_size, ".bed"), output_folder))
  }
  system2(command="mv", args= c(paste0(output_folder, '/', gsub(".bed", "", peak_file), "_bins.bed"), paste0(output_folder, '/', sample_prefix, "_peak_size.txt")))
}
