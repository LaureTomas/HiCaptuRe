#' Digest genome for a specific restriction enzyme
#'
#' This function takes a genome installed and generates its digest for a given restriction enzyme
#'
#' @param genome character with the name of the genome version
#' @param RE_name name of the restriction enzyme
#' @param motif recognition motif of the restriction enzyme
#' @param cut_position cut position of the restriction enzyme inside the motif
#' @param select_chr a character vector containing the specific chromosomes to used from this genome, if NULL all chromosomes will be used
#' @param PAR_mask a logical value where the Y chromosome should exclude the Pseudoautosomical Regions (PAR) or not
#' @param PAR_file a full path to a file containing the coordinates of Y chromosome PAR with at least 3 colums with header: seqnames, start, end
#' @param ... extra arguments for read.table
#'
#' @return list object with 2 elements
#'
#' @note The package provides for a PAR coordinates file only for Homo sapiens for the genome version 38
#'
#'
#' @importFrom Biostrings replaceAt DNAStringSet matchPattern
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqnames
#'
#'
#' @export
digest_genome <- function(genome="GRCh38",RE_name="hindIII",motif="AAGCTT",cut_position=1,select_chr=c(1:22,"X","Y"),PAR_mask=T,PAR_file=NULL,...)
{
  genome_name <- genome
  ## Load reference genome from installed ones
  genome <-  tryCatch({
    load_genome(genome = genome)
  }, error = function(e) {
    stop(e)
  })
  ## If Pseudoautosomal Regions (PAR) masked, read provided file for human
  if (PAR_mask)
  {
    if (is.null(PAR_file))
    {
      PAR_file <- system.file("extdata", paste0("PAR_",gsub(" ","_",metadata(genome)$organism),"_coordinates.txt"), package="HiCaptuRe")
      if (file.exists(PAR_file))
      {
        PAR <- read.delim(PAR_file,header = T,...)
      } else
      {
        stop(paste("There is no PAR file provided for",metadata(genome)$organism,"\n It must be a headed file with seqnames,start,end columns"))
      }
    } else
    {
      if (file.exists(PAR_file))
      {
        PAR <- read.delim(PAR_file,header = T,...)
      } else
      {
        stop("The PAR file provided doesn't exist")
      }
    }
  } ## PAR mask

  ## Select primary chromosomes
  chrs <- GenomeInfoDb::seqnames(genome)

  if (!is.null(select_chr))
  {
    chrs <- chrs[chrs %in% select_chr]
  }
  if (length(chrs) == 0)
  {
    stop("Chromosomes selected in 'select_chr' are not present in this genome. Please try changing 'select_chr' or setting it to 'NULL'")
  }

  p <- progressr::progressor(steps = length(chrs))


  ## Digest genome by chromosomes
  digest <- data.frame()
  for (chr in chrs)
  {
    p(sprintf(paste("Digesting",chr)))
    if (PAR_mask)
    {
      if (grepl(unique(PAR$seqnames),chr))
      {
        chr_seq <- Biostrings::replaceAt(genome[[chr]], IRanges::IRanges(PAR$start[1], PAR$end[1]),
                                         Biostrings::DNAStringSet(strrep("N",length(PAR$start[1]:PAR$end[1]))))
        chr_seq <- Biostrings::replaceAt(chr_seq, IRanges::IRanges(PAR$start[2], PAR$end[2]),
                                         Biostrings::DNAStringSet(strrep("N",length(PAR$start[2]:PAR$end[2]))))
      }
      chr_seq <- genome[[chr]]
    }
    else
    {
      chr_seq <- genome[[chr]]
    }

    m <- Biostrings::matchPattern(motif, chr_seq)

    correct <- start(m)-1+cut_position

    starts <- c(1,correct+1)
    ends <- c(correct,length(chr_seq))

    df <- data.frame(seqnames=chr,start=starts,end=ends)

    digest <- rbind(digest,df)
  }

  if (is.null(PAR_file))
  {
    PAR_file <- "NULL"
  } else
  {
    PAR_file <- normalizePath(PAR_file)
  }

  digest$fragment_ID <- 1:nrow(digest)
  output <- list(digest=digest,
                 parameters=c("Genome"=genome_name,
                              "Genome_Package"=genome@pkgname,
                              "Restriction_Enzyme"=RE_name,
                              "Motif"=motif,
                              "Cut_Position"=cut_position,
                              "Selected_Chromosomes"=paste(select_chr,collapse = ","),
                              "PAR_mask"=PAR_mask,
                              "PAR_file"=PAR_file),
                 seqinfo=seqinfo(genome)[chrs])
  return(output)
}


