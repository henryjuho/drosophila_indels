# Reviewer Comments
 
- [x] 1. Separate INDEL analysis into insertions and deletions
- [x] 2. Consider DFE for different lengths
    - [ ] run model per length bin, maybe 1+2bp, 3bp, 4+5bp, 6bp? using the a class model
- [x] 3. run with high classes c = 10
    - [x] see if power analysis in Kousathanas and Keightley 2013 or Keightley and Eyre-Walker 2010
    - [x] Keightley and Eyre-Walker 2010 table 5 
        - as sample size increases better able to distinguish models with different numbers of classes based on likelihood
        - best fitting model was 4 classes of sites, yet for n=16 and n=64, only 1 and 2 class models would converge
        - n=256 3 class would converge
        - n=1024 4 class would converge
    - [ ] if not do power analysis, simulate 2 class data with increasing pop size, see how often 1 class model is rejected at each sample size
- [x] 4. table 4 statistics on nonsense mutations would be a good comparison
    - nonsense pi: = 5.83e-6, theta = 9.13e-6, Tajima's D = -1.5
    - del cds  pi: = 2.74e-5, theta = 4.17e-5, Tajima's D = -1.48
    - ins cds  pi: = 1.47e-5, theta = 1.88e-5, Tajima's D = -0.94
- [x] 5. Small contrast between frame shift and in frame, yet frame shift more disruptive
    - [x] frame shift occur at higher frequency and are more deleterious, in frame occur at lower frequency but persist longer, result is that estimates are closer together
- [x] 6. Allele frequency distribution for INDELs
    - [x] add AR, see [plot](regional_indel_sfs.pdf)
- [ ] 7. Mutation rate influence pol error rate, yet neutral ref doesn't all pol error to vary
- [ ] 8. Discuss that insertions look neutral but dels don't and compare to Leushkin et al.
- [x] 9. Synonymous as neutral ref
    - [x] could use fourfold or ARs
    - [ ] check length distribution of AR INDELs
- [x] 10. Page 2 column 2 inherit should be inherent
- [ ] 11. Page 2 column 2 start of paragraph 'An additional challenge...' is poorly written
 
