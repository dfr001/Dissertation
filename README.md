# Dissertation
Code used when writing my dissertation, "[Wilcoxon Rank Sum Tests to Detect One-Sided Mixture Alternatives in Group Sequential Clinical Trials](https://escholarship.org/uc/item/3f28b3mc)."

The `grp.seq.mix` function is used to determine arm size needed for a study based on a given design alternative. See `r_code_doc.pdf` for full documentation for `grp.seq.mix`.

`cluster_est.R` contains code used for simulating estimation of parameters in the one-sided mixture alterative group sequential setting as well as automation of the creation of tables able to be copied and pasted into LaTeX. First, a function is defined to simulate a number of estimates for a particular design alternative scenario. Then a number of scenarios are run and stored. Finally, the stored results are read and put into a form to be easily transferred to LaTeX.
