## README
Scripts and workflows to produce the analyses in

Ganz, Luquette, Bizzotto et al. Contrasting somatic mutation patterns in aging human neurons and oligodendrocytes. Cell, 2024.

## Output directories and figure correspondence
Figure numbering in these scripts does not correspond to the published numbering.
Panel numberings do not correspond to the final arrangements, but descriptive names
can be used to match.

Main figure mappings:
* fig1 -> Figure 1
    * panel C: fig1_corrected
* fig2_sigprofilerextractor -> Figure 2
    * 
* figX -> Figure 3
* fig3_sigprofilerextractor -> Figure 4
* fig4_corrected -> Figure 5
    * panels A and C: fig4
* fig5_sigprofilerextractor -> Figure 6
* fig6 -> Figure 6
    * panel A: fig6_corrected

Supplemental figure mappings:
Figure S1 was created manually.
* suppfigX -> Figure S2
* suppfig2 -> Figure S3, Figure S4A-C,F
* suppfig2_corrected -> Figure S4D-E
* fig2_corrected (panel_c_all_A.pdf) -> Figure S4G
* suppfig3_corrected -> Figure S5
* suppfig4 -> Figure S6 
* suppfig5 -> Figure S7

**_corrected** directories indicate correction for location-specific sensitivity.
Analyses that were corrected for sensitivity were performed with and without
correction for comparison.

**_sigprofilerextractor** directories indicate fitting of mutational signatures
by SigProfilerExtractor. This replaced an earlier ad-hoc method.
