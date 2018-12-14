# How to re-built RNA velocity package

## reference: http://kbroman.org/pkg_primer/pages/build.html

## git clone the rna velocyto repository to local directory:
## git clone https://github.com/velocyto-team/velocyto.R.git /Users/yanwengong/Documents/fall_2018/nie_lab/projects/RNA_velocity/scripts

## change the name of velocity.R to velocyto.YG.R to avoid overwriting in the following places:
## 1. directory name; 2. DESCRIPTION: Package:*; 3. NAMESPACE: useDynLib(*) 4. man/velocyto.YG.R-package.Rd 

## terminal&at where the folder is: R CMD build velocyto.YG.R
## then, there will be a tar: velocyto.YG.R_0.6.tar.gz

## terminal: R CMD INSTALL velocyto.YG.R_0.6.tar.gz
## the intalled library will be under: /Library/Frameworks/R.framework/Versions/3.5/Resources/library/ 

library(velocyto.YG.R)
