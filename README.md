# Pre-processing

### Count tables
- We make sure that column and row names make sense
    - To be easy to use with the functionalities included here we require *unique* column names in the form ```<Temperature>_<Treatment>_<Time>_<Replicate>```.
    Note that this particular way of structuring the metadata stems from the biological questions at hand!

    Some comments regarding the **non-statistical analysis** normalization (this is not relevant for DEA): 
    
    Control treatment (Mock) must be named "control". There is some ambiguity here. Either we work with mock treatments as control per temperature, in which case temperature effect could not be seen, or we choose one temperature as reference. 

    In what follows I assume that 21 degrees is the "control" temperature, as I think it makes the most sense. So that we should call it "21C".

    Hence, the references w.r.t control will be taken w.r.t the 21 degrees mock treatment. In that sense, mock treatment at elevated temperature will get its own set of samples!

    Since we do not care about the time effects (in themselves), each time point gets its own reference.

- The row names must be coherent with the ID in the annotation .gff file. We will work with TAIR identifiers (e.g., AT1G22750).

- Counts must be integers, so *round* the data provided by Salmon. (For now we are ignoring the tximport issue, to discuss)

### DNA & annotation
- Apparently we do not have chromosome names but "contigs". We will work with that
- We will format the raw TAIR12 (from keygene) into simple IDs.

From the raw data, we can run, to construct our formatted annotation. We will be using this extensively. Code for this in ...

# Analysis (no formal hypohtesis testing!)
## RNA-seq processing

We assume that rows and column names are formatted as expected, this is fundamental as every function here follows the same structure. As an example
we will be following the older UU data, and by temperature all samples will have 21 degrees (there is nov variation).

For the following "non-statistical" analysis we will work on log2FC w.r.t control (to discuss). This is an easy to interpret metric that corrects for gene length and irrelevant biological processes present in control and treatment.

First, we need to make the samples comparable, correcting for sequence depth. To do so we follow DESeq2 Median-of-Ratios normalization. It is relevant to stop
here for a second to consider this, as it will be *key* in understanding what DESeq2 is doing when we talk about statistics.

### Median of ratios normalization

Step 1. Compute geometric mean per gene $i$ across $m$ samples $j$: 

$$K_i^{R} = \left( \prod_{j=1}^{m} K_{ij} \right)^{1/m}.$$

Step 2. Compute the ratios of counts (per gene $i$ and sample $j$) to geometric means (replace nan for 0):
$$\text{ratio} = \frac{K_{ij}}{K_i^{R}}$$
Step 3. Compute the median ratio per sample (over the genes whose geometric mean is not 0), the size factors:

$$s_j = \text{median}_{i:K_i^{R}\neq 0}\left(\frac{K_{ij}}{K_i^{R}}\right)$$

Step 4. Normalize counts

$$K^{\text{norm}}_{ij} = \frac{K_{ij}}{s_j}$$

We will not discuss today why this is better than something else. As a general justification, doing this is consistent with the normalization done by DESeq2 later on the statistics. 

Note, however, this **does not correct for gene length**. Hence, this is still not a valid metric to compare across genes.

To do a PCA we need the genes to be comparable across samples, but also, we need the genes to be in a similar "scale". In any case, we are going to avoid a lot of headches by working on log2FC.

### $\log_2 \text{FC}$ construction

By focusing on the change of a gene with respect to itself, we remove gene length bias. Furthermore, we remove some other irrelevant sources of variation.

Sample names follow 
```<Temperature>_<Treatment>_<Time>_<Replicate>```. Let $K^{\text{norm}}_{i,(T,\tau,t,r)}$ be the median-of-ratios normalized count of gene $i$ at temperature $T$, treatment $\tau$, time $t$ and replicate $r$.

**Common reference**. Both temperature and treatment matter, so every sample is compared against the same reference condition, defined by a control temperature $T_0$ and a control treatment $\tau_0$ (the control_temperature and control_treatment arguments). Replicate numbers are arbitrary labels, so we do not pair replicates; we average the reference over its $n_0$ replicates. **Time is the one factor where we do not need a common reference, so the reference is taken at the same time point $t$ as the sample being compared**:

$$\bar{K}_{i,t}^{\text{ref}} = \frac{1}{n_0} \sum_{r'=1}^{n_0} K^{\text{norm}}_{i,(T_0,\,\tau_0,\,t,\,r')} .$$

Log fold change. For every sample that is not part of the reference, with a pseudocount of 1 so that zero counts are well defined:

$$\log_2 \text{FC}_{i,(T,\tau,t,r)} = \log_2 \left( \frac{K^{\text{norm}}_{i,(T,\tau,t,r)} + 1}{\bar{K}_{i,t}^{\text{ref}} + 1} \right).$$

Values are then clipped to $[-c, c]$ (default $c = 4$) so that a few extreme genes do not dominate the PCA:

$$\log_2 \text{FC}^{\text{clip}} = \max\!\big(-c,\ \min(c,\ \log_2 \text{FC})\big).$$

Code to do this from raw counts, assuming proper previous formating in code/Counts_Processing/utils.py.

#### NOTE! Alternatively, TPM normalization might be used before PCA. But never plotting without correction for sequencing depth and gene length!!
There is already a TPM file in the (cloud) folders!

## PCA

With the normalized data, this is absolutely strightfoward. I provide code for some nice classic PCA plot on the log2FC data (or any other normalization for that matter) in code/Analysis/PCA.py.

#### NOTE! This is very different from the PCA "recommended" in the DESeq2 pipeline.

## Sparse RDA (*Advanced*)
(Discuss if we want to get into this now!)

This is much more advanced than PCA. It gives similar insight but it also explains the regulatory motifs underpining the variation in a biplot. Strictly speaking is not a PCA plot (and mathematically it gets somewhat tricky). Overall this is a prediction model, an optimization process constrained by some hyperparameters 
that are choosen using CV.

We will not get into the math of this. For today, one can look at this as an extension of PCA that uses promoter regions and takes more time to compute because we are tuning some hyperparameters.

This requires installation of the MotifCompass package (which is not available if I do not send you invitation to collaborate). I provide the code to get the biplot in code/Analysis/sRDA.py. However, note, that much more can be done with this!

### Further requirements
Here we need genome information. 

We will construct 2. One based on 6-mers, which does require extra steps to map it to TFs, another based on the JASPAR's PWM, which already indicates 
the TF. We do not use DAP-seq because that mapping was done in TAIR10 and it does not seem easy to use it on TAIR12.

# Statistical analysis
Now we switch gears and programming language!

## Motif Enrichment Analysis





## Further issues to discuss

- What to do with the untreated samples? In both experiments