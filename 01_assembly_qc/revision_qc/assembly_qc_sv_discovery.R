#!/usr/bin/env Rscript
# Minimal, reviewer-targeted assembly metric versus SV discovery workflow.
# Usage:
# Rscript assembly_qc_sv_discovery_minimal.R haplotype_sv.tsv assembly_qc_minimal.tsv output_directory

suppressPackageStartupMessages({library(data.table); library(ggplot2); library(glmmTMB)})
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 3L) stop("Usage: Rscript assembly_qc_sv_discovery_minimal.R haplotype_sv.tsv assembly_qc_minimal.tsv output_directory")
sv <- fread(a[[1]]); qc <- fread(a[[2]]); out <- a[[3]]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
need <- function(x, nms, label) { z <- setdiff(nms, names(x)); if (length(z)) stop(label, " missing: ", paste(z, collapse=", ")) }
need(sv, c("Haplotype","Sample","Gender","Population","INS","DEL","INStelo","DELtelo","INScenSat","DELcenSat"), "SV table")
need(qc, c("Sample","AssemblySpanBp","ContigN50Bp","ConsensusQV","Sex","Population","BroadRegion"), "QC table")
if (anyDuplicated(qc$Sample)) stop("QC table has duplicate Sample values")
for (x in c("INS","DEL","INStelo","DELtelo","INScenSat","DELcenSat")) { sv[[x]] <- as.numeric(sv[[x]]); if(anyNA(sv[[x]]) || any(sv[[x]]<0)) stop("Invalid ",x) }
h <- sv[, .(n_hap=uniqueN(Haplotype), n_sex=uniqueN(Gender), n_pop=uniqueN(Population)), by=Sample]
if (any(h$n_hap != 2L | h$n_sex != 1L | h$n_pop != 1L)) stop("Every sample must have exactly two haplotypes and consistent metadata")
d <- sv[, .(Sex_from_sv=first(Gender), Population_from_sv=first(Population),
             SV_INS=sum(INS), SV_DEL=sum(DEL), INS_TELO=sum(INStelo), DEL_TELO=sum(DELtelo),
             INS_CENSAT=sum(INScenSat), DEL_CENSAT=sum(DELcenSat)), by=Sample]
d[, `:=`(SV_INS_NONDIFFICULT=SV_INS-INS_TELO-INS_CENSAT, SV_DEL_NONDIFFICULT=SV_DEL-DEL_TELO-DEL_CENSAT)]
if (any(d$SV_INS_NONDIFFICULT<0 | d$SV_DEL_NONDIFFICULT<0)) stop("Regional SV classes overlap; do not subtract them until their definitions are reconciled")
d[, `:=`(SV_TOTAL=SV_INS+SV_DEL, SV_TOTAL_NONDIFFICULT=SV_INS_NONDIFFICULT+SV_DEL_NONDIFFICULT)]
d <- merge(d, qc, by="Sample", all.x=TRUE, sort=FALSE)
if (anyNA(d$Sex)) stop("QC rows missing for: ", paste(d[is.na(Sex),Sample], collapse=", "))
if (any(d$Sex != d$Sex_from_sv | d$Population != d$Population_from_sv)) stop("Sex or Population disagrees between inputs")
for (x in c("AssemblySpanBp","ContigN50Bp","ConsensusQV")) { d[[x]] <- as.numeric(d[[x]]); if(anyNA(d[[x]]) || any(d[[x]]<=0)) stop("Invalid QC metric: ",x) }
d[, `:=`(Sex=factor(Sex), Population=factor(Population), BroadRegion=factor(BroadRegion),
          z_log_span=as.numeric(scale(log10(AssemblySpanBp))), z_log_n50=as.numeric(scale(log10(ContigN50Bp))),
          z_qv=as.numeric(scale(ConsensusQV)))]
fwrite(d,file.path(out,"01_sample_level_dataset.tsv"),sep="\t")
fwrite(data.table(item=c("N_samples","N_populations","N_broad_regions"),value=c(nrow(d),uniqueN(d$Population),uniqueN(d$BroadRegion))),file.path(out,"00_audit.tsv"),sep="\t")

# One QC metric per model: these metrics are correlated and are not interpreted as independent causal effects.
rhs <- function(metric, x) { p <- c(metric,"Sex"); if(nlevels(x$BroadRegion)>1L) p<-c(p,"BroadRegion"); if(nlevels(x$Population)>1L) p<-c(p,"(1|Population)"); paste(p,collapse=" + ") }
fit_one <- function(x,outcome,metric,subset) {
  x <- droplevels(as.data.frame(x)); f <- as.formula(paste(outcome,"~",rhs(metric,x)))
  model <- tryCatch(glmmTMB(f,data=x,family=nbinom2()),error=function(e)e)
  if(inherits(model,"error")) return(data.table(Subset=subset,Outcome=outcome,Metric=metric,N=nrow(x),Estimate=NA_real_,SE=NA_real_,CountRatio=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,Convergence=NA_integer_,PositiveDefiniteHessian=NA,Status=conditionMessage(model)))
  co <- summary(model)$coefficients$cond; b<-co[metric,"Estimate"]; s<-co[metric,"Std. Error"]
  data.table(Subset=subset,Outcome=outcome,Metric=metric,N=nrow(x),Estimate=b,SE=s,CountRatio=exp(b),CI_low=exp(b-1.96*s),CI_high=exp(b+1.96*s),P=co[metric,"Pr(>|z|)"],Convergence=model$fit$convergence,PositiveDefiniteHessian=model$sdr$pdHess,Status="OK")
}
outcomes <- c("SV_TOTAL","SV_INS","SV_DEL","SV_TOTAL_NONDIFFICULT","SV_INS_NONDIFFICULT","SV_DEL_NONDIFFICULT")
metrics <- c("z_log_span","z_log_n50","z_qv")
r <- rbindlist(lapply(outcomes, function(o) {
  rbindlist(lapply(metrics, function(m) fit_one(d, o, m, "All")))
}), fill=TRUE)
r[,FDR:=p.adjust(P,"BH")]
fwrite(r,file.path(out,"03_negative_binomial_associations.tsv"),sep="\t")

labels <- c(z_log_span="Diploid assembled span",z_log_n50="Contig N50",z_qv="Consensus QV")
main <- r[Subset=="All" & Outcome %in% c("SV_TOTAL","SV_INS","SV_DEL") & Status=="OK"]
main[,MetricLabel:=factor(labels[Metric],levels=unname(labels))]
main[,Outcome:=factor(Outcome,levels=c("SV_TOTAL","SV_INS","SV_DEL"),labels=c("Total SV","Insertion","Deletion"))]
p <- ggplot(main,aes(MetricLabel,CountRatio,colour=Outcome))+geom_hline(yintercept=1,linetype=2,colour="grey50")+geom_errorbar(aes(ymin=CI_low,ymax=CI_high),position=position_dodge(.65),width=.2)+geom_point(position=position_dodge(.65),size=2.4)+scale_y_log10()+theme_classic(base_size=12)+theme(axis.text.x=element_text(angle=30,hjust=1))+labs(x=NULL,y="SV discovery count ratio per 1-SD metric increase (95% CI)",colour=NULL,subtitle="Adjusted for sex, broad region, and population random intercept")
ggsave(file.path(out,"04_main_association_forest.pdf"),p,width=10,height=5)
long <- rbindlist(list(data.table(Value=d$AssemblySpanBp,SV_TOTAL=d$SV_TOTAL,BroadRegion=d$BroadRegion,Metric="Diploid assembled span"),data.table(Value=d$ContigN50Bp,SV_TOTAL=d$SV_TOTAL,BroadRegion=d$BroadRegion,Metric="Contig N50")))
p <- ggplot(long,aes(Value,SV_TOTAL,colour=BroadRegion))+geom_point(alpha=.75,size=2)+geom_smooth(method="lm",formula=y~x,colour="black",se=TRUE)+facet_wrap(~Metric,scales="free_x")+scale_x_log10(labels=scales::label_number(scale_cut=scales::cut_si("bp")))+scale_y_continuous(labels=scales::label_number(big.mark=","))+theme_classic(base_size=12)+labs(x=NULL,y="Per-sample SV count",colour="Broad region",subtitle="Descriptive plots; adjusted models are reported separately")
ggsave(file.path(out,"05_descriptive_span_n50_scatter.pdf"),p,width=10,height=4.6)
message("Done: ",normalizePath(out))

