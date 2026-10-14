## =====================================================================








## =====================================================================
source(file.path(Sys.getenv("EHR_AUDIT_REPO_ROOT", unset = getwd()),
                 "R", "common", "00_paths.R"))
suppressPackageStartupMessages({
  for(p in c("data.table","mice","cluster","diptest","MASS","glmnet","mclust"))
    if(!requireNamespace(p,quietly=TRUE)) install.packages(p,repos="https://cloud.r-project.org")
  library(data.table); library(mice); library(cluster); library(diptest); library(MASS)
  library(glmnet); library(mclust) })
set.seed(20240601)
OUT_DIR<-file.path(DIR_OUTPUT,"qc"); dir.create(OUT_DIR,showWarnings=FALSE,recursive=TRUE)

MICE_M<-5L; MICE_MAXIT<-10L; PRIMARY_IMP<-1L; BOOT_B<-500L
REUSE_EXISTING_MICE<-TRUE
FEATURES33<-c("bilirubin_total_max","alt_max","pao2fio2ratio_min","abs_lymphocytes_min","lactate_max",
 "ph_min","pco2_max","calcium_min","calcium_max","ptt_max","inr_max","temperature_min","temperature_max",
 "urine_output_24h_ml","glucose_max","aniongap_max","potassium_min","potassium_max","hemoglobin_min",
 "sodium_min","sodium_max","wbc_max","platelets_min","bicarbonate_min","chloride_min","chloride_max",
 "bun_max","creatinine_max","resp_rate_max","gcs_min","spo2_min","heart_rate_max","mbp_min")
BASE_NUM<-c("age","sex","sofa_score","aki_stage_0_24h")
winsor<-function(v,p=c(.01,.99)){q<-quantile(v,p,na.rm=TRUE);v[v<q[1]]<-q[1];v[v>q[2]]<-q[2];v}
KM<-function(m,k,ns=25)kmeans(m,k,nstart=ns,iter.max=100,algorithm="Lloyd")
ci2<-function(m){km<-KM(m,2);km$tot.withinss/sum(sweep(m,2,colMeans(m))^2)}
dip_axis_p<-function(m){km<-KM(m,2);diptest::dip.test(as.numeric(m%*%(km$centers[1,]-km$centers[2,])))$p.value}
gap_rise<-function(m,Kmax=12,B=50,nsub=3000){mm<-if(nrow(m)>nsub)m[sample.int(nrow(m),nsub),]else m
  cg<-clusGap(mm,function(x,k)KM(x,k,10),K.max=Kmax,B=B,d.power=2,spaceH0="scaledPCA");as.logical(which.max(cg$Tab[,"gap"])==Kmax)}
sig_eff<-function(m,B=100){ci<-ci2(m);Sig<-cov(m);cn<-replicate(B,ci2(MASS::mvrnorm(nrow(m),colMeans(m),Sig)));100*(mean(cn)-ci)/mean(cn)}
fauc<-function(y,p){ok<-!is.na(y)&!is.na(p);y<-y[ok];p<-p[ok];n1<-sum(y==1);n0<-sum(y==0);if(n1==0||n0==0)return(NA);r<-rank(p);(sum(r[y==1])-n1*(n1+1)/2)/(n1*n0)}

## ---------- 1. load + winsor + MICE ----------
d<-fread(PATH_FINAL_FULL); mid<-intersect(c("stay_id","hadm_id"),names(d))[1]
bc<-fread(file.path(DIR_DATA,"baseline_covars.csv"))
core<-d[,c(mid,"age","gender","sofa_score","mortality_30d","make30",FEATURES33),with=FALSE]
core[,sex:=as.integer(toupper(substr(as.character(gender),1,1))=="M")]
core<-merge(core,bc[,c(mid,"aki_stage_0_24h"),with=FALSE],by=mid,all.x=TRUE)
core<-core[complete.cases(core[,c("age","sex","sofa_score","aki_stage_0_24h","mortality_30d","make30"),with=FALSE])]
RAW<-copy(core); for(j in FEATURES33) RAW[[j]]<-winsor(RAW[[j]])          # winsor 1/99 before MICE
mice_path<-file.path(DIR_MODEL,"mice_primary.rds")
if (isTRUE(REUSE_EXISTING_MICE) && file.exists(mice_path)) {
  cat(sprintf("n=%d ; loading existing MICE object: %s\n", nrow(RAW), mice_path))
  imp<-readRDS(mice_path)
  if (!inherits(imp, "mids")) stop("Existing mice_primary.rds is not a mids object.")
  if (imp$m < MICE_M) stop("Existing MICE object has fewer imputations than MICE_M. Delete mice_primary.rds or lower MICE_M.")
} else {
  cat(sprintf("n=%d ; running MICE (PMM, m=%d, maxit=%d) on %d features ...\n",nrow(RAW),MICE_M,MICE_MAXIT,length(FEATURES33)))
  imp<-mice(RAW[,..FEATURES33],m=MICE_M,maxit=MICE_MAXIT,method="pmm",printFlag=FALSE,seed=2024)
  saveRDS(imp,mice_path)
}

## ---------- 2. per-imputation: standardize, cluster, align to mortality ----------
y30<-RAW$mortality_30d
lab_list<-vector("list",MICE_M); Xstd_list<-vector("list",MICE_M)
for(mm in seq_len(MICE_M)){
  Xc<-scale(as.matrix(complete(imp,mm)[,FEATURES33])); Xstd_list[[mm]]<-Xc
  km<-KM(Xc,2); cl<-km$cluster
  hi<-which.max(tapply(y30,cl,mean)); lab<-ifelse(cl==hi,1L,2L)   # 1 = higher 30d mortality (C1)
  lab_list[[mm]]<-lab
}
## cross-imputation stability (ARI)
aris<-c(); for(a in 1:(MICE_M-1)) for(b in (a+1):MICE_M) aris<-c(aris,mclust::adjustedRandIndex(lab_list[[a]],lab_list[[b]]))
cat(sprintf("\n[Stability] Paired across-imputation ARI: mean=%.3f (range %.3f–%.3f)  <- higher=stable assignment (stability is not discreteness)\n",mean(aris),min(aris),max(aris)))


disc<-rbindlist(lapply(seq_len(MICE_M),function(mm){Xc<-Xstd_list[[mm]]
  data.table(imputation=mm, dip_axis_p=signif(dip_axis_p(Xc),3),
             gap_rises_to_max=gap_rise(Xc), sig_effect_pct=round(sig_eff(Xc),3))}))
cat("\n===== Discreteness: separate results for 5 imputations (historical expectations: high dip P / gap TRUE / effect~0.2-0.3%)=====\n"); print(disc)

fwrite(disc,file.path(OUT_DIR,"mice_discreteness_byimp.csv"))
fwrite(data.table(pairwise_ARI=aris),file.path(OUT_DIR,"mice_ari.csv"))

## ---------- 4. B2 decisive contrast per imputation + Rubin pooling ----------

b2_rubin<-function(outcome){
  th<-Vp<-numeric(MICE_M); auc_b<-auc_bp<-auc_fp<-auc_fpl<-numeric(MICE_M)
  for(mm in seq_len(MICE_M)){
    Xc<-Xstd_list[[mm]]; y<-RAW[[outcome]]
    BASE<-as.matrix(RAW[,..BASE_NUM]); lab<-factor(lab_list[[mm]])
    fold<-integer(length(y)); for(cl in c(0,1)){ix<-which(y==cl);fold[ix]<-sample(rep_len(1:10,length(ix)))}
    cvp<-function(dz){p<-numeric(length(y));for(f in 1:10){tr<-fold!=f;te<-fold==f
      m<-cv.glmnet(dz[tr,],y[tr],family="binomial",alpha=.5,nfolds=5);p[te]<-as.numeric(predict(m,dz[te,],s="lambda.min",type="response"))};p}
    Dbase<-model.matrix(~.,data.frame(BASE))[,-1,drop=FALSE]
    Dbp  <-model.matrix(~.,data.frame(BASE,lab))[,-1,drop=FALSE]
    Dfp  <-cbind(BASE,Xc)
    Dfpl <-cbind(BASE,Xc,model.matrix(~lab)[,-1,drop=FALSE])
    p_b<-cvp(Dbase);p_bp<-cvp(Dbp);p_fp<-cvp(Dfp);p_fpl<-cvp(Dfpl)
    auc_b[mm]<-fauc(y,p_b);auc_bp[mm]<-fauc(y,p_bp);auc_fp[mm]<-fauc(y,p_fp);auc_fpl[mm]<-fauc(y,p_fpl)
    th[mm]<-fauc(y,p_fpl)-fauc(y,p_fp)
    d<-replicate(BOOT_B,{i<-sample.int(length(y),replace=TRUE);fauc(y[i],p_fpl[i])-fauc(y[i],p_fp[i])});Vp[mm]<-var(d)
  }
  thbar<-mean(th);W<-mean(Vp);Bv<-var(th);Tt<-W+(1+1/MICE_M)*Bv;se<-sqrt(Tt)
  data.table(outcome=outcome,
    auc_base=round(mean(auc_b),4),auc_base_phen=round(mean(auc_bp),4),
    auc_feat=round(mean(auc_fp),4),auc_feat_lab=round(mean(auc_fpl),4),
    incr_base_to_basephen=round(mean(auc_bp)-mean(auc_b),4),
    decisive_dAUC_feat_plus_label=round(thbar,4),
    dAUC_lo=round(thbar-1.96*se,4),dAUC_hi=round(thbar+1.96*se,4))
}
cat("\n===== B2 (Rubin-pooled across ",MICE_M," imputations) =====\n",sep="")
b2<-rbind(b2_rubin("mortality_30d"),b2_rubin("make30")); print(b2)

## ---------- 5. save primary MICE matrix + labels (imp#PRIMARY) with alignment check ----------
Xp<-Xstd_list[[PRIMARY_IMP]]; sid<-RAW[[mid]]
rc<-cor(rank(Xp[,"lactate_max"]),rank(RAW[["lactate_max"]]),use="complete.obs")
if(rc<0.9) stop("Alignment check failed (rank correlation < 0.9)")
Xsave<-as.data.table(as.data.frame(Xp)); Xsave[,stay_id:=as.integer(sid)]; setcolorder(Xsave,c("stay_id",FEATURES33))
saveRDS(Xsave,file.path(DIR_MODEL,"X_primary33_std_mice.rds"))
saveRDS(data.table(stay_id=as.integer(sid),cluster_k2=lab_list[[PRIMARY_IMP]]),file.path(DIR_MODEL,"labels_primary_mice.rds"))

fwrite(disc,file.path(OUT_DIR,"mice_discreteness_byimp.csv"))
fwrite(b2,file.path(OUT_DIR,"mice_b2_rubin.csv"))
fwrite(data.table(pairwise_ARI=aris),file.path(OUT_DIR,"mice_ari.csv"))
cat(sprintf("\n[Saved] X_primary33_std_mice.rds / labels_primary_mice.rds (imp#%d, alignment rankcor=%.3f)\n",PRIMARY_IMP,rc))
cat("Downstream input roles: PREPROC_MATRIX_RDS = X_primary33_std_mice.rds; labels = labels_primary_mice.rds.\n")
cat("Historical interpretation rule: a pooled feat_lab - feat delta AUC interval containing 0 does not establish added label information.\n")
