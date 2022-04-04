
# STAR alignment
```
        ==========     _____ _    _ ____  _____  ______          _____  
        =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
	  v2.0.2

//========================== featureCounts setting ===========================\\
||                                                                            ||
||             Input files : 24 BAM files                                     ||
||                                                                            ||
||                           GRCm39_104_25T_S1_trimmedAligned.sortedByCoo ... ||
||                           GRCm39_104_26DHT_S11_trimmedAligned.sortedBy ... ||
||                           GRCm39_104_27E2_S13_trimmedAligned.sortedByC ... ||
||                           GRCm39_104_DHT-18_S16_trimmedAligned.sortedB ... ||
||                           GRCm39_104_DHT-8_S14_trimmedAligned.sortedBy ... ||
||                           GRCm39_104_DHT-9_S15_trimmedAligned.sortedBy ... ||
||                           GRCm39_104_E2-10_S9_trimmedAligned.sortedByC ... ||
||                           GRCm39_104_E2-11_S10_trimmedAligned.sortedBy ... ||
||                           GRCm39_104_E2-19_S12_trimmedAligned.sortedBy ... ||
||                           GRCm39_104_KO-16_S17_trimmedAligned.sortedBy ... ||
||                           GRCm39_104_KO-21_S18_trimmedAligned.sortedBy ... ||
||                           GRCm39_104_KO-23_S19_trimmedAligned.sortedBy ... ||
||                           GRCm39_104_KO-24_S20_trimmedAligned.sortedBy ... ||
||                           GRCm39_104_Testosteron-17_S4_trimmedAligned. ... ||
||                           GRCm39_104_Testosteron-5_S2_trimmedAligned.s ... ||
||                           GRCm39_104_Testosteron-6_S3_trimmedAligned.s ... ||
||                           GRCm39_104_Ukapszula-22_S8_trimmedAligned.so ... ||
||                           GRCm39_104_Ukapszula-2_S5_trimmedAligned.sor ... ||
||                           GRCm39_104_Ukapszula-3_S6_trimmedAligned.sor ... ||
||                           GRCm39_104_Ukapszula-4_S7_trimmedAligned.sor ... ||
||                           GRCm39_104_intakt-13_S21_trimmedAligned.sort ... ||
||                           GRCm39_104_intakt-14_S22_trimmedAligned.sort ... ||
||                           GRCm39_104_intakt-15_S23_trimmedAligned.sort ... ||
||                           GRCm39_104_intakt-20_S24_trimmedAligned.sort ... ||
||                                                                            ||
||             Output file : featureCounts_GRCm39.104_O.txt                   ||
||                 Summary : featureCounts_GRCm39.104_O.txt.summary           ||
||              Paired-end : no                                               ||
||        Count read pairs : no                                               ||
||              Annotation : Mus_musculus.GRCm39.104.gtf (GTF)                ||
||      Dir for temp files : ./                                               ||
||                                                                            ||
||                 Threads : 1                                                ||
||                   Level : meta-feature level                               ||
||      Multimapping reads : not counted                                      ||
|| Multi-overlapping reads : counted                                          ||
||   Min overlapping bases : 1                                                ||
||                                                                            ||
\\============================================================================//

//================================= Running ==================================\\
||                                                                            ||
|| Load annotation file Mus_musculus.GRCm39.104.gtf ...                       ||
||    Features : 842117                                                       ||
||    Meta-features : 55416                                                   ||
||    Chromosomes/contigs : 39                                                ||
||                                                                            ||
|| Process BAM file GRCm39_104_25T_S1_trimmedAligned.sortedByCoord.out.ba ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 22144066                                             ||
||    Successfully assigned alignments : 10198011 (46.1%)                     ||
||    Running time : 0.36 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_26DHT_S11_trimmedAligned.sortedByCoord.out ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 2556540                                              ||
||    Successfully assigned alignments : 964715 (37.7%)                       ||
||    Running time : 0.04 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_27E2_S13_trimmedAligned.sortedByCoord.out. ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 80782635                                             ||
||    Successfully assigned alignments : 48314947 (59.8%)                     ||
||    Running time : 1.29 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_DHT-18_S16_trimmedAligned.sortedByCoord.ou ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 5664966                                              ||
||    Successfully assigned alignments : 2326641 (41.1%)                      ||
||    Running time : 0.09 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_DHT-8_S14_trimmedAligned.sortedByCoord.out ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 2809763                                              ||
||    Successfully assigned alignments : 1261576 (44.9%)                      ||
||    Running time : 0.05 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_DHT-9_S15_trimmedAligned.sortedByCoord.out ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 3631095                                              ||
||    Successfully assigned alignments : 1488043 (41.0%)                      ||
||    Running time : 0.06 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_E2-10_S9_trimmedAligned.sortedByCoord.out. ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 4229242                                              ||
||    Successfully assigned alignments : 1771850 (41.9%)                      ||
||    Running time : 0.07 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_E2-11_S10_trimmedAligned.sortedByCoord.out ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 4111532                                              ||
||    Successfully assigned alignments : 1650072 (40.1%)                      ||
||    Running time : 0.07 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_E2-19_S12_trimmedAligned.sortedByCoord.out ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 4955260                                              ||
||    Successfully assigned alignments : 2062911 (41.6%)                      ||
||    Running time : 0.08 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_KO-16_S17_trimmedAligned.sortedByCoord.out ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 5253954                                              ||
||    Successfully assigned alignments : 2338536 (44.5%)                      ||
||    Running time : 0.09 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_KO-21_S18_trimmedAligned.sortedByCoord.out ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 1671477                                              ||
||    Successfully assigned alignments : 1058783 (63.3%)                      ||
||    Running time : 0.03 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_KO-23_S19_trimmedAligned.sortedByCoord.out ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 5635547                                              ||
||    Successfully assigned alignments : 3239875 (57.5%)                      ||
||    Running time : 0.10 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_KO-24_S20_trimmedAligned.sortedByCoord.out ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 3390658                                              ||
||    Successfully assigned alignments : 1359871 (40.1%)                      ||
||    Running time : 0.05 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_Testosteron-17_S4_trimmedAligned.sortedByC ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 1909175                                              ||
||    Successfully assigned alignments : 768262 (40.2%)                       ||
||    Running time : 0.03 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_Testosteron-5_S2_trimmedAligned.sortedByCo ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 3625419                                              ||
||    Successfully assigned alignments : 1576848 (43.5%)                      ||
||    Running time : 0.06 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_Testosteron-6_S3_trimmedAligned.sortedByCo ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 820246                                               ||
||    Successfully assigned alignments : 342358 (41.7%)                       ||
||    Running time : 0.02 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_Ukapszula-22_S8_trimmedAligned.sortedByCoo ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 3939648                                              ||
||    Successfully assigned alignments : 1654796 (42.0%)                      ||
||    Running time : 0.06 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_Ukapszula-2_S5_trimmedAligned.sortedByCoor ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 2871650                                              ||
||    Successfully assigned alignments : 1195361 (41.6%)                      ||
||    Running time : 0.05 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_Ukapszula-3_S6_trimmedAligned.sortedByCoor ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 2876855                                              ||
||    Successfully assigned alignments : 1184307 (41.2%)                      ||
||    Running time : 0.05 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_Ukapszula-4_S7_trimmedAligned.sortedByCoor ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 1413231                                              ||
||    Successfully assigned alignments : 587388 (41.6%)                       ||
||    Running time : 0.02 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_intakt-13_S21_trimmedAligned.sortedByCoord ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 3244407                                              ||
||    Successfully assigned alignments : 1320932 (40.7%)                      ||
||    Running time : 0.05 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_intakt-14_S22_trimmedAligned.sortedByCoord ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 2241221                                              ||
||    Successfully assigned alignments : 902681 (40.3%)                       ||
||    Running time : 0.04 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_intakt-15_S23_trimmedAligned.sortedByCoord ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 1522951                                              ||
||    Successfully assigned alignments : 617807 (40.6%)                       ||
||    Running time : 0.03 minutes                                             ||
||                                                                            ||
|| Process BAM file GRCm39_104_intakt-20_S24_trimmedAligned.sortedByCoord ... ||
||    Single-end reads are included.                                          ||
||    Total alignments : 2873608                                              ||
||    Successfully assigned alignments : 1184396 (41.2%)                      ||
||    Running time : 0.05 minutes                                             ||
||                                                                            ||
|| Write the final count table.                                               ||
|| Write the read assignment summary.                                         ||
||                                                                            ||
|| Summary of counting results can be found in file "featureCounts_GRCm39.10  ||
|| 4_O.txt.summary"                                                           ||
||                                                                            ||
\\============================================================================//

```

