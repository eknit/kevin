#' read_in_data

read_in_data <- function(runfile_id){
  library(plyr)
  library(lattice)
  #Add line to find whether subdirectory already exists
  readin <- read.csv(paste(runfile_id, "/", runfile_id, ".csv", sep=""), header=F, skip=7, stringsAsFactors=T)[, c(1:3, 5, 6,9:11)]
  
  
  split_point <- match("Drift Corrected", readin[,1]) #Find where spreadsheet divides between drift corrected and not drift corrected
  raw <- readin[1:(split_point-2),]
  dc <- readin[(split_point+3):nrow(readin), c(1, 4:8)]
  
  data <- merge(raw, dc, by="V1")
  names(data) <- c("Ps", "ID", "Wt", "NugR", "d15NR", "CugR", "d13CR", "d18OR", "Nugdc", "d15Ndc", "Cugdc", "d13Cdc", "d18Odc")
  
  #make numeric things numeric
  data[, c(1, 3:13)] <- sapply(data[, c(1, 3:13)], function(x) as.numeric(as.character(x))) 
  
  #Add pcC, pcN and CN ratio						#NB these are based on drift corrected values
  data$pcC <- data$Cugdc/data$Wt/10
  data$pcN <- data$Nugdc/data$Wt/10 
  data$CN <- data$Cugdc/data$Nugdc*14/12
  data$Runfile <- runfile_id
  data
}

#' make_archive()
make_archive <- function(){
  archive <- data.frame(samplename=data$ID, sampleweight = data$Wt, RS="S", Owner=Owner, PreparedBy =PreparedBy, Funding=Funding, Sampletype=0, Taxa="", site="", age.context=0, ugN=data$Nugdc, d15N=data$d15Ndc, ugC=data$Cugdc, d13C=data$d13Cdc, CN=data$CN, labno = paste(runfile_id, data$Ps, sep="-"))
  archive <- archive[order(data$Ps),]
  list.of.salanines <- data$Ps[archive$samplename=="SALANINE"]
  first.s <- list.of.salanines[1]
  first.r <- first.s-8
  list.of.ralanines <- c(2, first.r, list.of.salanines)
  archive$RS <- "S"
  archive$RS[list.of.ralanines] <- "R"
  
  Sampletype <- rep(SampleType, nrow(archive))
  standards <- which(archive$samplename %in% c("ALANINE", "SALANINE", RM1.name, RM2.name))
  Sampletype[standards] <- "Standard"
  archive$Sampletype <- Sampletype
  
  header <- matrix(nrow=3, ncol=16)
  header[1, 1:2] <- c("Run number", paste(runfile_id))
  header[2, 1:2] <- c("User", "")
  header[1:2, 3:16] <- ""
  header[3,] <-  c("sample name", "sample weight", "R/S", "Owner", "Prepared by", "Funding", "Sample type", "Taxa", "Site", "age/context", "ug N", "d15N AIR", "ug C", "d13C VPDB", "C/N molar", "Lab number")
  final_archive <- rbind(header, as.matrix(archive))
  colnames(final_archive) <- NULL
  write.csv(final_archive, paste(runfile_id, "/final.archive.", runfile_id, ".csv", sep=""), row.names=F)
}
#' make standards_true
make_standards_true <- function(){
  standards_true <- data.frame(rbind(
    c(RM1.name, RM1T.C, RM1Tsd.C, RM1T.N, RM1Tsd.N), 
    c("SALANINE", -26.91, 0, -1.63, 0),
    c(RM2.name, RM2T.C, RM2Tsd.C, RM2T.N, RM2Tsd.N)))
  standards_true[,2:5] <- sapply(standards_true[,2:5], function(x) as.numeric(as.character(x)))
  standards_true
}
#' make_normalized

make_normalized <- function(data){
  if (!is.null(delete)) {data <- subset(data, !(Ps %in% delete))} #Remove lines to be deleted
  standards <- which(data$ID %in% c("ALANINE", "SALANINE", RM1.name, RM2.name))
  standards_subset <- data[standards, ]
  standards_subset <- standards_subset[standards_subset$ID != "ALANINE" & standards_subset$Ps > 8,]
  standards_meas <- ddply(standards_subset, "ID", function (x) 
    c(meas_d13C=mean(x$d13Cdc), meas_sd=sd(x$d13Cdc), meas_d15N=mean(x$d15Ndc), meas_sd=sd(x$d15Ndc)))
  #standards_true <- data.frame(rbind(
   # c(RM1.name, RM1T.C, RM1Tsd.C, RM1T.N, RM1Tsd.N), 
   # c("SALANINE", -26.91, 0, -1.63, 0),
   # c(RM2.name, RM2T.C, RM2Tsd.C, RM2T.N, RM2Tsd.N)))
 # standards_true[,2:5] <- sapply(standards_true[,2:5], function(x) as.numeric(as.character(x)))
  standards_true <- make_standards_true()
 
  get_run_values <- function(isotope){
    delta_dc <- data[, names(data)==paste(isotope, "dc", sep="")]
    RM1M <- mean(delta_dc[data$ID==RM1.name])
    RM2M <- mean(delta_dc[data$ID==RM2.name]) 
    RM1Msd <- sd(delta_dc[data$ID==RM1.name])
    RM2Msd <- sd(delta_dc[data$ID==RM2.name])
    alaninesd <- sd(delta_dc[data$ID=="SALANINE" & data$Ps > Ps_cutoff])
    mean.alanine <- mean(delta_dc[data$ID=="SALANINE" & data$Ps > Ps_cutoff])
    return(rbind(RM1M, RM2M, RM1Msd, RM2Msd, alaninesd, mean.alanine))
  }
  
  
  normalize_x <- function(x, isotope){
    if (isotope=="d13C"){
      index=2
    }
    else if (isotope=="d15N"){
      index=4
    }
    RM1Tsd <- standards_true[1, index+1]
    RM2Tsd <- standards_true[3, index+1]
    RM1T <- standards_true[1, index]
    RM2T <- standards_true[3, index]
    RM1M <- get_run_values(isotope)[1]
    RM2M <- get_run_values(isotope)[2]
    RM1Msd <- get_run_values(isotope)[3]
    RM2Msd <- get_run_values(isotope)[4]
    alaninesd <- get_run_values(isotope)[5]
    measured_vals <- c(RM1T, RM1M, RM2T, RM2M, x)
    matrix1 <- diag(c(RM1Tsd, RM1Msd, RM2Tsd, RM2Msd, alaninesd))
    matrix2 <- matrix(rep(measured_vals, 5), ncol=5)
    matrix3 <- matrix1 + matrix2
    raw2true <- function(cv) {cv[1] + (cv[5]-cv[2])*((cv[1]-cv[3])/	(cv[2]-cv[4]))}
    normalized <- raw2true(measured_vals)
    finalerror <- sqrt(sum((normalized - apply(matrix3, 2, raw2true))^2))
    return(c(normalized, finalerror))
  }
  
  
  data$normd13C <- sapply(data$d13Cdc, function(x) normalize_x(x, "d13C"))[1,]
  data$d13Csd <- sapply(data$d13Cdc, function(x) normalize_x(x, "d13C"))[2,]
  data$normd15N <- sapply(data$d15Ndc, function(x) normalize_x(x, "d15N"))[1,]
  data$d15Nsd <- sapply(data$d15Ndc, function(x) normalize_x(x, "d15N"))[2,]
  data <- data[order(data$Ps),]
  write.csv(data, paste(runfile_id, "/Alldata.csv", sep=""))
  data
}
#' d13C
d13C <- expression(paste(delta^{13},"C (\u2030) drift corrected"))

#' d15N
d15N <- expression(paste(delta^{15},"N (\u2030) drift corrected"))

#' Make standards subset
#' 
make_standards_subset <- function(data){
  standards <- which(data$ID %in% c("ALANINE", "SALANINE", RM1.name, RM2.name))
  standards_subset <- data[standards, ]
  standards_subset <- standards_subset[standards_subset$ID != "ALANINE" & standards_subset$Ps > 8,]
  standards_subset$ID <- factor(standards_subset$ID, levels=c(RM1.name, "SALANINE", RM2.name))
  standards_subset
}

make_samples <- function(data){
  standards <- which(data$ID %in% c("ALANINE", "SALANINE", RM1.name, RM2.name))
  samples <- data[-standards,]
  samples
}

#' make_carbon_standard_figure
make_carbon_standards_figure <- function(){
  RM1 <- data[data$ID==paste(RM1.name),]
  RM2 <- data[data$ID==paste(RM2.name),]
  alanine <- data[data$ID=="SALANINE" & data$Ps > 10,]
  standards_subset <- make_standards_subset(data)
  samples <- make_samples(data)
  #standards <- which(data$ID %in% c("ALANINE", "SALANINE", RM1.name, RM2.name))
  #standards_subset <- data[standards, ]
  #standards_subset <- standards_subset[standards_subset$ID != "ALANINE" & standards_subset$Ps > 8,]
  #standards_subset$ID <- factor(standards_subset$ID, levels=c(RM1.name, "SALANINE", RM2.name))
  #samples <- data[-standards,]
  standards_true <- data.frame(rbind(
    c(RM1.name, RM1T.C, RM1Tsd.C, RM1T.N, RM1Tsd.N), 
    c("SALANINE", -26.91, 0, -1.63, 0),
    c(RM2.name, RM2T.C, RM2Tsd.C, RM2T.N, RM2Tsd.N)))
  standards_true[,2:5] <- sapply(standards_true[,2:5], function(x) as.numeric(as.character(x)))
  
  ###### For Carbon 
  par(mfrow=c(1,3))
  par(mar=c(4,5,10,1))
  
  for (i in 1:3){
    batch <- standards_subset[standards_subset$ID==levels(standards_subset$ID)[i],]
    ymax <- max(c(batch$d13Cdc, (standards_true[i,2]-standards_true[i,3]), (standards_true[i,2]+standards_true[i,3])))
    ymin <- min(c(batch$d13Cdc, (standards_true[i,2]-standards_true[i,3]), (standards_true[i,2]+standards_true[i,3])))
    boxplot(batch$d13Cdc, ylab=d13C, ylim=c(ymin, ymax))
    par(new=T)
    plot(rep(1, nrow(batch)), batch$d13Cdc, axes=F, xlab="", ylab="", ylim=c(ymin, ymax))
    abline(h=standards_true[i, 2], col="red", lwd=2)
    abline(h=standards_true[i, 2]+standards_true[i, 3], col="red", lty=2)
    abline(h=standards_true[i, 2]-standards_true[i, 3], col="red", lty=2)
    text(0.8, standards_true[i, 2], "Expected value", col="red", pos=3)
    text(1.1, batch$d13Cdc, batch$Ps)
    mtext(batch$ID[1])
  }
  par(new=T)
  par(mfrow=c(1,1))
  plot(1:10, 1:10, axes=F, xlab="", ylab="", type="n")
  mtext("Drift corrected d13C for Two Reference Materials and S-Alanine", line=7)
  mtext("Look at the spread of the measured values of References Materials and Alanines. What is the spread in permil?", cex=0.7, line=6, adj=0)
  mtext("Are the measured deltas of the alanines centred (randomly distributed) around their expected mean, in red?", cex=0.7, line=5.5, adj=0)
  mtext("Or is there a bias in the measurement? Are any of the alanines wrong? Position in run printed for identification", cex=0.7, line=5, adj=0)
  dev.copy2pdf(file=paste(runfile_id, "/plot1.pdf", sep=""), encoding="WinAnsi")
}

#' make_nitrogen_standards_figure

make_nitrogen_standards_figure <- function(){
  RM1 <- data[data$ID==paste(RM1.name),]
  RM2 <- data[data$ID==paste(RM2.name),]
  alanine <- data[data$ID=="SALANINE" & data$Ps > 10,]
  standards_subset <- make_standards_subset(data)
  samples <- make_samples(data)
  #standards <- which(data$ID %in% c("ALANINE", "SALANINE", RM1.name, RM2.name))
  #standards_subset <- data[standards, ]
  #standards_subset <- standards_subset[standards_subset$ID != "ALANINE" & standards_subset$Ps > 8,]
  #standards_subset$ID <- factor(standards_subset$ID, levels=c(RM1.name, "SALANINE", RM2.name))
  #samples <- data[-standards,]
  standards_true <- make_standards_true()
  #standards_true <- data.frame(rbind(
  #  c(RM1.name, RM1T.C, RM1Tsd.C, RM1T.N, RM1Tsd.N), 
  #  c("SALANINE", -26.91, 0, -1.63, 0),
  #  c(RM2.name, RM2T.C, RM2Tsd.C, RM2T.N, RM2Tsd.N)))
  #standards_true[,2:5] <- sapply(standards_true[,2:5], function(x) as.numeric(as.character(x)))
  
  par(mfrow=c(1,3))
  par(mar=c(4,5,10,1))
  
  for (i in 1:3){
    batch <- standards_subset[standards_subset$ID==levels(standards_subset$ID)[i],]
    ymax <- max(c(batch$d15Ndc, (standards_true[i,4]-standards_true[i,5]), (standards_true[i,4]+standards_true[i,5])))
    ymin <- min(c(batch$d15Ndc, (standards_true[i,4]-standards_true[i,5]), (standards_true[i,4]+standards_true[i,5])))
    boxplot(batch$d15Ndc, ylab=d15N, ylim=c(ymin, ymax))
    par(new=T)
    plot(rep(1, nrow(batch)), batch$d15Ndc, axes=F, xlab="", ylab="", ylim=c(ymin, ymax))
    abline(h=standards_true[i, 4], col="red", lwd=2)
    abline(h=standards_true[i, 4]+standards_true[i, 5], col="red", lty=2)
    abline(h=standards_true[i, 4]-standards_true[i, 5], col="red", lty=2)
    text(0.8, standards_true[i, 4], "Expected value", col="red", pos=3)
    text(1.1, batch$d15Ndc, batch$Ps)
    mtext(batch$ID[1])
  }
  par(new=T)
  par(mfrow=c(1,1))
  plot(1:10, 1:10, axes=F, xlab="", ylab="", type="n")
  mtext("Drift corrected d15N for Two Reference Materials and S-Alanine", line=7)
  mtext("Look at the spread of the measured values of References Materials and Alanines. What is the spread in permil?", cex=0.7, line=6, adj=0)
  mtext("Are the measured deltas of the alanines centred (randomly distributed) around their expected mean, in red?", cex=0.7, line=5.5, adj=0)
  mtext("Or is there a bias in the measurement? Are any of the alanines wrong? Position in run printed for identification", cex=0.7, line=5, adj=0)
  dev.copy2pdf(file=paste(runfile_id,"/plot2.pdf", sep=""), encoding="WinAnsi")
}

#' make_drift_correction_figure
#' 
make_drift_correction_figure <- function(){
  par(mfrow=c(2,2))
  par(mar=c(4,5,2,1))
  xtext <- "Position in Run"
  plot(data$Ps, data$normd13C, xlab=xtext, ylab=d13C, main="Change in d13C through run?")
  plot(data$Ps, data$normd15N, xlab=xtext, ylab=d15N, main="Change in d15N through run?" )
  D15N <- expression(paste("Raw - drift-corrected ", Delta^{15},"N (\u2030)"))
  D13C <- expression(paste("Raw - drift-corrected ", Delta^{13},"C (\u2030)"))
  plot(data$Ps, data$d13Cdc-data$d13CR, ylab=D13C, xlab=xtext, main="Effect of drift correction vs position")
  plot(data$Ps, data$d15Ndc-data$d15NR, ylab=D15N, xlab=xtext, main="Effect of drift correction vs position")
  dev.copy2pdf(file=paste(runfile_id,"/plot3.pdf", sep=""), encoding="WinAnsi")
}

#' make_carbon_normalization_figure

make_carbon_normalization_figure <- function(){
  standards_subset <- make_standards_subset(data)
  samples <- make_samples(data)
  standards_true <- make_standards_true()
  standards_meas <- ddply(standards_subset, "ID", function (x) 
    c(meas_d13C=mean(x$d13Cdc), meas_sd=sd(x$d13Cdc), meas_d15N=mean(x$d15Ndc), meas_sd=sd(x$d15Ndc)))
  par(mfrow=c(1,1))
  par(mar=c(5,5,5,5))
  ylims=c(min(standards_true[,2]-1), max(standards_true[,2]+1)) 
  xlims=c(min(standards_meas[,2]-1), max(standards_meas[,2]+1))
  plot(standards_meas[,2], standards_true[,2], type="n", ylab=expression(paste("True ", delta^{13},"C (\u2030)")), 
       xlab=expression(paste("Measured ", delta^{13},"C (\u2030)")), ylim=ylims, xlim=xlims)
  arrows(standards_meas[c(1,3),2]-standards_meas[c(1,3),3], standards_true[c(1,3),2], standards_meas[c(1,3),2]+standards_meas[c(1,3),3], standards_true[c(1,3),2], angle=90, code=3, length=0.01)
  arrows(standards_meas[c(1,3),2], standards_true[c(1,3),2]-standards_true[c(1,3),3], standards_meas[c(1,3),2], standards_true[c(1,3),2]+standards_true[c(1,3),3], angle=90, code=3, length=0.01)
  text(standards_meas[c(1,3),2], standards_true[c(1,3),2], c(paste(RM1.name), paste(RM2.name)), cex=0.6, adj=c(0,1))
  slope <- abs(standards_true[1,2]-standards_true[3,2])/abs(standards_meas[1,2]-standards_meas[3,2])
  intercept <- standards_true[1,2]-(standards_meas[1,2]*slope)
  abline(a=intercept, b=slope, col="red")
  legend("topleft", paste("slope =", round(slope,4), ", intercept =", round(intercept,4)), col="red", lty=1)
  par(new=T)
  plot(samples$d13Cdc, samples$normd13C, axes=F, ylab="", xlab="", ylim=ylims, xlim=xlims)
  mtext("Measured versus True (normalized) d13C", line=3)
  mtext("How much effect does the normalization regression have? How different is the slope from 1,", line=2, cex=0.7)
  mtext("how different is the intercept from zero? Are the error bars for the two reference materials as small", line=1, cex=0.7)
  mtext("as they should be? Do your samples lie on the regression line between the two reference materials?", line=0, cex=0.7)
  dev.copy2pdf(file=paste(runfile_id,"/plot4.pdf", sep=""), encoding="WinAnsi")
}

#' make_nitrogen_normalization_figure
make_nitrogen_normalization_figure <- function(){
  standards_subset <- make_standards_subset(data)
  samples <- make_samples(data)
  standards_true <- make_standards_true()
  standards_meas <- ddply(standards_subset, "ID", function (x) 
    c(meas_d13C=mean(x$d13Cdc), meas_sd=sd(x$d13Cdc), meas_d15N=mean(x$d15Ndc), meas_sd=sd(x$d15Ndc)))
  
  par(mfrow=c(1,1))
  par(mar=c(5,5,5,5))
  ylims=c(min(standards_true[,4]-1), max(standards_true[,4]+1)) 
  xlims=c(min(standards_meas[,4]-1), max(standards_meas[,4]+1))
  plot(standards_meas[,4], standards_true[,4], type="n", ylab=expression(paste("True ", delta^{15},"N (\u2030)")), 
       xlab=expression(paste("Measured ", delta^{15},"N (\u2030)")), ylim=ylims, xlim=xlims)
  arrows(standards_meas[c(1,3),4]-standards_meas[c(1,3),5], standards_true[c(1,3),4], standards_meas[c(1,3),4]+standards_meas[c(1,3),5], standards_true[c(1,3),4], angle=90, code=3, length=0.01)
  arrows(standards_meas[c(1,3),4], standards_true[c(1,3),4]-standards_true[c(1,3),5], standards_meas[c(1,3),4], standards_true[c(1,3),4]+standards_true[c(1,3),5], angle=90, code=3, length=0.01)
  text(standards_meas[c(1,3),4], standards_true[c(1,3),4], c(paste(RM1.name), paste(RM2.name)), cex=0.6, adj=c(0,1))
  slope <- abs(standards_true[1,4]-standards_true[3,4])/abs(standards_meas[1,4]-standards_meas[3,4])
  intercept <- standards_true[1,4]-(standards_meas[1,4]*slope)
  abline(a=intercept, b=slope, col="red")
  ylims=c(min(standards_true[,4]-1), max(standards_true[,4]+1)) 
  xlims=c(min(standards_meas[,4]-1), max(standards_meas[,4]+1))
  legend("topleft", paste("slope =", round(slope,4), ", intercept =", round(intercept,4)), col="red", lty=1)
  par(new=T)
  plot(samples$d15Ndc, samples$normd15N, axes=F, ylab="", xlab="", ylim=ylims, 
       xlim=xlims)
  mtext("Measured versus True (normalized) d15N", line=3)
  mtext("How much effect does the normalization regression have? How different is the slope from 1,", line=2, cex=0.7)
  mtext("how different is the intercept from zero? Are the error bars for the two reference materials as small", line=1, cex=0.7)
  mtext("as they should be? Do your samples lie on the regression line between the two reference materials?", line=0, cex=0.7)
  dev.copy2pdf(file=paste(runfile_id,"/plot4.pdf", sep=""), encoding="WinAnsi")
}


#' make_CN_figure

make_CN_figure <- function(){
  standards_subset <- make_standards_subset(data)
  samples <- make_samples(data)
  par(mfrow=c(2,2))
  par(mar=c(4,5,3,1))
  CNmin <- if (min(samples$CN) < 2.9) {min(samples$CN)} else { 2.8}
  CNmax <- if (max(samples$CN) > 3.6) {max(samples$CN)} else { 3.7}
  plot(samples$CN, samples$normd13C, xlim=c(CNmin, CNmax), xlab="C/N ratio", ylab=d13C)
  abline(v=2.9, col="red")
  abline(v=3.6, col="red")
  plot(samples$CN, samples$normd15N, xlim=c(CNmin, CNmax), xlab="C/N ratio", ylab=d15N)
  abline(v=2.9, col="red")
  abline(v=3.6, col="red")
  xtext <- expression(paste("Sample weight (", mu,"g)"))
  plot(samples$CugR, samples$normd13C, ylab=d13C, xlab=xtext)
  plot(samples$NugR, samples$normd15N, ylab=d15N, xlab=xtext)
  par(mfrow=c(1,1))
  mtext("Do the C/N ratios fall within a good range?", cex=0.7, line=2)
  mtext("Within the good range, is there any possibility of a trend between C/N and d13C or d15N?", cex=0.7, line=1)
  par(mfrow=c(2,1))
  mtext("Similarly, was the weight of C and N measured high enough?", cex=0.7, line=2)
  mtext("Is there any biasing effect of weight on isotopic ratio, even for the OK samples?", cex=0.7, line=1)
  dev.copy2pdf(file=paste(runfile_id,"/plot4.pdf", sep=""), encoding="WinAnsi")
}

#' final_compile
final_compile <- function(){
  data <- read_in_data(runfile_id)
  data <- make_normalized(data)
  make_carbon_standards_figure()
  make_nitrogen_standards_figure()
  make_drift_correction_figure()
  make_carbon_normalization_figure()
  make_nitrogen_normalization_figure()
  make_CN_figure()
  samples <- make_samples(data)
  samples <- samples[order(samples$Ps),]
  write.csv(samples, paste(runfile_id, "/Samples.", runfile_id, ".csv", sep=""))
}

#' make_merged()

make_merged <- function(directory){
  setwd(directory)
  files <- dir("~/Dropbox/dropbox AGRICURB/Runfiles") 
  read_data  <- function(files, x){
    this.file <- paste(files[x], "/Samples.", files[x], ".csv", sep="")
    if (file.exists(this.file)){
      read.csv(this.file,  stringsAsFactors=F)
    }
  }
  result <- suppressWarnings(lapply(files, read_data))
  cond <- sapply(result, function(x) length(x)>1)
  not_null <- result[cond]
  merged_all <- Reduce(function(x, y) merge(x, y, all=TRUE), not_null)
  merged_all
}

#' get_codes
#' Gets ID codes (first three letters is default, but can change)
get_codes <- function(df, start=1, stop=3){
  codes <- substring(df$ID, start, stop)
}

#' extract 
#' Extracts specific sites codes - can take a vector of sites codes
extract <- function(df, codes){
  df$Site <- get_codes(df)
  sub_df <- df[df$Site %in% codes,]
  sub_df
}

#' merge_plants
#' Checks for C-only and N-only runs, merges them together keeping unique column ids

merge_plants <- function(data){
  Ndata <- data[data$normd13C==0 | is.na(data$normd13C),]
  for (i in 1:nrow(Ndata)){
    if (is.na(Ndata$RunfileN)[i]) Ndata$RunfileN[i] <- Ndata$Runfile[i]
  }
  names(Ndata) <- paste(names(Ndata), "N", sep="")
  names(Ndata)[names(Ndata)=="RunfileN"] <- "Runfile"
  names(Ndata)[names(Ndata)=="RunfileNN"] <- "RunfileN"
  names(Ndata)[names(Ndata)=="IDN"] <- "ID"
  names(Ndata)[names(Ndata)=="pcNN"] <- "pcN"
  names(Ndata)[names(Ndata)=="normd15NN"] <- "normd15N"
  names(Ndata)[names(Ndata)=="d15NsdN"] <- "d15Nsd"
  include <- c("PsN", "ID", "WtN", "NugRN", "d15NRN", "NugdcN", "d15NdcN", "pcN", "normd15N", "d15Nsd", "RunfileN")
  Nkeep <- Ndata[include]
  
  Cdata <- data[data$normd15N==0,]
  include <- c("Ps", "ID", "Wt", "CugR", "d13CR", "Cugdc", "d13Cdc", "pcC", "normd13C", "d13Csd", "Runfile")
  Ckeep <- Cdata[include]
  merged <- merge(Ckeep, Nkeep, by="ID", all=T)
  merged$CalcCN <- (merged$pcC*14)/(merged$pcN*12)
  merged
}

