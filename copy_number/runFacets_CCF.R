computeCCF=function(vaf, tt, minor, purity, multiplicity=NULL){

  bs=vaf
  cg=tt
  major=tt-minor
  vlen=length(major)

  out=diff1=rep(list(NA),vlen)
  for(i in 1:vlen){
    #multiplicity. e.g, major-minor:3-1 multiplicity include 1, 2, 3
    majori=major[i]
    if(!is.na(majori)){
      mm=1:major[i]
      out[[i]]=bs[i]*(2*(1-purity[i])+cg[i]*purity[i])/mm/purity[i]
      diff1[[i]]=abs(out[[i]]-1) #distance to 1
    }
  }


  if(is.null(multiplicity)){
    column.list=lapply(diff1,which.min)
    #catch integer(0) to maintain list length, NA allowed when unlist
    column.list[sapply(column.list, function(x)length(x)==0)] = NA
    column.list=unlist(column.list)
  }else{
    column.list = multiplicity
    }

  #select the multiplicity that gives ccf closest to 1
  #ccf=out[cbind(1:nrow(out), column.list)]

  rawccf=unlist(lapply(1:length(out),function(x)out[[x]][column.list[x]]))

  ccf=rawccf
  ccf[major==0]=NA
  ccf[ccf>1]=1
  ccf[is.na(minor)]=NA

  return(list(ccf=ccf, rawccf=rawccf, out=out, multiplicity=column.list))
}

library(Hmisc)
confCCF=function(alt, ref, tt, minor, purity, multiplicity, alpha=0.05){

  depth=alt+ref
  ci=binconf(alt, depth, alpha, method=c("exact"))
  lowerci=ci[,2]
  upperci=ci[,3]

  lower=computeCCF(vaf=lowerci, tt, minor, purity, multiplicity)$ccf
  upper=computeCCF(vaf=upperci, tt, minor, purity, multiplicity)$ccf

  return(list(lower=lower, upper=upper))

}
