# Missing Data Analysis

TODO(Gerard): Missing exploration, in processs

\label{sec:Missing_data_analysis}

In this section, the patterns of missing data will be discussed and analyzed, in order to assess its influence in the analyses performed in previous sections.

Figure \ref{fig:nabars} and Table \ref{taula:nas} show the missing data distribution over time and for each dose. Clearly, the fact that an observation is missing is associated with the dose and the time (fisher tests p-values= 0.005 and 0.007, respectively). Actually, there are more missings in low and medium dose than in high dose, and also, there are more missings as time increases. Missings in the PCV are not produced completely at random, then.

```{r nabars, fig.height=3, fig.width=5, fig.cap="Percentage of missings for each time and dose with respect to the whole sample. \\label{fig:nabars}", fig.pos="H"}

par(cex=0.8)
a <- barplot(table(is.na(cows$pcv.f), paste(cows$dose, cows$time, sep="_"))[2,c(4:9, 1:3)]/nrow(cows)*100, beside=T, space=c(0,0,0,1,0,0,1,0,0), xaxt="n", ylab="%")
abline(h=0)
axis(1, at=a[c(2,5,8),], paste("Dose", levels(cows$dose)), tick=F, line=1.5)
axis(1, at=a, rep(paste("T", 1:3, sep=""),3), tick=F, line=-0.5)
```

```{r taula:nas, results='asis'}
tau <- xtable(aggregate(list(Missing=is.na(cows$pcv.f), Available=!is.na(cows$pcv.f)), list(Time=cows$time, Dose=cows$dose), sum)[,c(2,1,3:4)],
              caption="Number of missings and available data in the outcome.",
              label="taula:nas")
print(tau, comment=F, include.rownames=F)
```

For all this, it is pretty clear that there is some systematic mechanism that produces the missing data. Probably, given that 8 of the 9 cows that have some missings have the first missing after being unhealthy, it would be logic to understand that the missings are produced because of the death or that the cow is dropped out from the study due to having too many health issues. This makes sense also when looking at Figure ????: the cow with ID 5, for example, would have died (or left the study) before taking the time 3 low dose measure, and if one assumes that its sequence of doses was H-L-M, this explains the missings for the medium dose. Following this rationale, cows number 1, 2 and 4 would have died at the time 3 of dose L (their sequence would end with L, since there are no more missings); cows number 3, 8 and 10 before time 3 of the medium dose (ending with M dose) and cows with ID 6 and 9 before time 2 of the L dose. In this sense, our sample could be easily biased because we only have data in the low and medium dose for the cows with a less severe infection (due to the death of the cows that have a worse prognostic and that they are more likely to die with a dose other than high).
Therefore, the pattern for missings in this data is, probably, missings not at random (MNAR). 

Since all the analyses performed previously were done with the complete cases (missing PCV's were omited), the results obtained could be far from the reality (or not). The complete cases analysis can only deal with MCAR type of missing data. Therefore, several approaches will be used to analyse the data taking this into account.
                                                                               
                                                                               It is natural to think that, if the animals that give up the study are due to deaths they are unhealthy for the non-observed times. This could justify (in some way) subtituting the missings for the last value observed (last observation carried forward). Also a two pattern mixed model will also be adjusted in order to try to deal better with the MNCAR.
                                                                               
                                                                               ```{r}
                                                                               # ultima <- sapply(1:10, function(x) cows[is.na(cows$pcv) & cows$id==x, c("time", "dose")], simplify=F)
                                                                               #   if(nrow(x)<3){
                                                                               #     cows[cows$time==(min(x$time)-1) & cows$dose==unique(x$dose),"pcv"]
                                                                               #   } else {
                                                                               #     dosi <- names(table(x$dose))[table(x$dose)!=3]
                                                                               #     cows[cows$time==(min(x$time[x$dose==dosi])-1) & cows$dose==dosi,"pcv"]                                                                                                                                                                        
                                                                               #   }
                                                                               # })
                                                                               # cows.locf <- 
                                                                               ```
                                                                               