#!/usr/bin/env python3
"""Model PRSS copy number with polygenic scores.

Fit the model in the reference population and compare predicted and observed
values in Bajau. Inputs include PLINK genotypes, phenotypes and covariates.
Outputs include model validation summaries and population comparisons.
"""

import pandas as pd
import numpy as np
from pandas_plink import read_plink
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
from scipy import stats
import matplotlib.pyplot as plt
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def main():
    """ Main analysis workflow """
    try:
        # Stage 1: Data preparation and quality control
        logger.info("=== Stage 1: Data preparation and QC ===")
        # 1.1 Read reference population data (Southeast Asian agricultural populations)
        logger.info("Load reference PLINK data ...")
        ref_geno, ref_fam, ref_snps = load_plink_data('../malaylike_plink/33Malay_like.GRCh38.dv161.chr1_22XY.norm.QC.SNP')
        # 1.2 Read target population data (Bajau)
        logger.info("Load Bajau PLINK data ...")
        bajau_geno, bajau_fam, bajau_snps = load_plink_data('../bajau_plink/32Bajau_like.GRCh38.dv161.chr1_22XY.norm.QC.SNP')
        # 1.3 Read phenotype data (Gene copy number)
        logger.info("Load phenotype data ...")
        ref_pheno = load_phenotype_data('../Malay.PRSS_CNV.persample.txt', ref_fam.values)
        bajau_pheno_observed = load_phenotype_data('../Bajau.PRSS_CNV.persample.txt', bajau_fam.values)
        # 1.4 Read covariate data
        logger.info("Load covariate data ...")
        ref_covariates = load_covariates('../65merge.cov.txt', ref_fam.values)
        # 1.5 Find shared SNPs (Use the same SNPs in the reference and target populations)
        common_snps = list(set(ref_geno.columns) & set(bajau_geno.columns))
        # Keep autosomal SNPs only (common_snps format: chr1_41981_A_G)
        common_snps = [snp for snp in common_snps if snp.split('_')[0] in [f'chr{i}' for i in range(1, 23)]]
        logger.info(f"Variants shared by the reference and target populations:  {len(common_snps)}  SNPs")
        # 2025-08-25 20:58:22,927 - INFO - Variants shared by the reference and target populations: 11872363 SNPs
        ref_geno_common = ref_geno.loc[:, ref_geno.columns.intersection(common_snps)]
        bajau_geno_common = bajau_geno.loc[:, bajau_geno.columns.intersection(common_snps)]
        # Stage 2: Fit the PGS model in the reference population
        logger.info("=== Stage 2: Fit the PGS model in the reference population  ===")
        # 2.1 Run GWAS(Adjust for covariates)
        logger.info("Run GWAS ( Adjust for covariates )...")
        gwas_results = perform_gwas_with_covariates(
            genotype_df=ref_geno_common,
            phenotype_s=ref_pheno,
            covariates_df=ref_covariates,
            min_n=15
        )
        # 2.2 Select significant SNPs(Typically use P value threshold)
        significant_snps = select_significant_snps(gwas_results, pval_threshold=5e-8)
        logger.info(f"SNPs with P < 5e-8:  {len(significant_snps)}  significant SNP")
        # Extract significant SNP weights and genotypes
        snp_weights = gwas_results.loc[significant_snps, 'beta']
        ref_geno_sig = ref_geno_common[significant_snps]
        bajau_geno_sig = bajau_geno_common[significant_snps]
        # 2.3 Calculate polygenic scores
        logger.info("Calculate polygenic scores ...")
        ref_pgs = calculate_pgs(ref_geno_sig, snp_weights)
        bajau_pgs_predicted = calculate_pgs(bajau_geno_sig, snp_weights)
        # 2.4 Validate PGS model performance in the reference population
        logger.info("Validate PGS model performance ...")
        validate_pgs_model(ref_pheno, ref_pgs, trait_name="PRSS")
        # Stage 3: Compare predicted and observed values in Bajau
        logger.info("=== Stage 3: Bajau population analysis  ===")
        # 3.1 Statistical comparison: Paired t-test and linear regression
        compare_predicted_observed(
            predicted_pgs=bajau_pgs_predicted,
            observed_pheno=bajau_pheno_observed,
            trait_name="PRSS"
        )
        logger.info("Analysis complete !")
    except Exception as e:
        logger.error(f"Analysis failed : {str(e)}")
        raise

def load_plink_data(prefix):
    """
     Read PLINK binary files
    Parameters:
    -----------
    prefix : str
        PLINK file prefix ( Without the extension )
    Returns:
    --------
    geno_df : DataFrame
         Genotype matrix ( Samples  x SNP)
    bim : DataFrame
        SNP Information
    fam : DataFrame
         Sample metadata
    """
    (bim, fam, G) = read_plink(prefix)
    geno_array = G.compute()  # shape (num_snps, num_samples)
    geno_df = pd.DataFrame(geno_array.T, columns=bim['snp'].values, index=fam['iid'].values)
    # Create variant map
    snps = pd.DataFrame({
        'SNP': bim['snp'].values,
        'CHR': bim['chrom'].values,
        'POS': bim['pos'].values
    })
    return geno_df, fam['iid'], snps

def load_phenotype_data(file_path, sample_ids):
    """
     Load phenotype data
    Parameters:
    -----------
    file_path : str
         Phenotype file path
    sample_ids : array-like
         Samples ID List ,  Used to align records
    Returns:
    --------
    phenotype_s : Series
         Phenotype data ,  Indexed by sample ID
    """
    pheno_df = pd.read_csv(file_path, sep='\t', index_col=0)
    # Align phenotype and genotype sample order
    phenotype_s = pheno_df.loc[sample_ids].iloc[:, 0]
    phenotype_s.name = 'copy_number'
    return phenotype_s

def load_covariates(file_path, sample_ids):
    """
     Load and standardize covariates
    Parameters:
    -----------
    file_path : str
         Covariate file path
    sample_ids : array-like
         Samples ID List ,  Used to align records
    Returns:
    --------
    covariates_df : DataFrame
         Standardized covariate matrix
    """
    covar_df = pd.read_csv(file_path, sep='\t', index_col=0)
    covariates_df = covar_df.loc[sample_ids]
    # Remove age Column
    if 'age' in covariates_df.columns:
        covariates_df = covariates_df.drop(columns=['age'])
    # Standardize continuous covariates (Age, PCs)
    continuous_vars = ['pc1', 'pc2', 'pc3', 'pc4', 'pc5']
    for var in continuous_vars:
        if var in covariates_df.columns:
            scaler = StandardScaler()
            covariates_df[var] = scaler.fit_transform(covariates_df[[var]].values)
    return covariates_df

def perform_gwas_with_covariates(genotype_df, phenotype_s, covariates_df, min_n=15):
    """
     Vectorized  GWAS,  Adjust for covariates
    Parameters
    ----------
    genotype_df : DataFrame (samples x SNPs)
    phenotype_s : Series
    covariates_df : DataFrame
    min_n : int,  Minimum sample size
    Returns
    -------
    gwas_results : DataFrame
        index = SNP
        columns = ['beta','se','pval']
    """
    # Align samples
    common_samples = genotype_df.index.intersection(phenotype_s.index).intersection(covariates_df.index)
    G = genotype_df.loc[common_samples].to_numpy(dtype=float)   # n x m
    y = phenotype_s.loc[common_samples].to_numpy(dtype=float)   # n
    C = covariates_df.loc[common_samples].to_numpy(dtype=float) # n x k
    snps = genotype_df.columns
    n, m = G.shape
    k = C.shape[1]
    if n < min_n:
        return pd.DataFrame(columns=['beta','se','pval'])
    # Build the projection matrix M = I - Q Q^T, Remove covariate effects
    # First QR Decompose the covariate matrix + Intercept
    Q, R = np.linalg.qr(np.column_stack([np.ones(n), C]))  # n x (k+1)
    M = np.eye(n) - Q @ Q.T                                 # n x n
    y_resid = M @ y  #  Phenotype residuals
    G_resid = M @ G  # SNP Matrix residuals , n x m
    # Exclude SNPs with zero standard deviation
    sd = np.std(G_resid, axis=0)
    valid = sd > 1e-8
    G_resid = G_resid[:, valid]
    snps_valid = snps[valid]
    # beta
    beta = (G_resid.T @ y_resid) / np.sum(G_resid * G_resid, axis=0)
    # Residuals
    resid_matrix = y_resid[:, None] - G_resid * beta[None, :]
    sigma2 = np.sum(resid_matrix ** 2, axis=0) / (n - k - 2)
    se = np.sqrt(sigma2 / np.sum(G_resid ** 2, axis=0))
    tval = beta / se
    pval = 2 * (1 - stats.t.cdf(np.abs(tval), df=n - k - 2))
    gwas_results = pd.DataFrame({
        'beta': beta,
        'se': se,
        'pval': pval
    }, index=snps_valid)
    return gwas_results

def select_significant_snps(gwas_df, pval_threshold=5e-8):
    """
     Select SNPs using a P-value threshold
    Parameters:
    -----------
    gwas_df : DataFrame
        GWAS Results
    pval_threshold : float
         Significance threshold
    Returns:
    --------
    significant_snps : list
         Significant SNP Name list
    """
    significant_snps = gwas_df[gwas_df['pval'] < pval_threshold].index.tolist()
    return significant_snps

def calculate_pgs(genotype_df, weights_s):
    """
     Calculate polygenic scores
    Parameters:
    -----------
    genotype_df : DataFrame
         Genotype matrix
    weights_s : Series
        SNP weights ( Effect sizes )
    Returns:
    --------
    pgs_s : Series
         Polygenic score results
    """
    # Ensure SNP Use the same order
    aligned_geno = genotype_df[weights_s.index]
    # PGS = Σ(SNP Genotypes * SNP weights)
    pgs = aligned_geno.dot(weights_s)
    return pgs

def validate_pgs_model(observed_pheno, predicted_pgs, trait_name):
    """
     Validate PGS Predictive performance in the reference population
    Parameters:
    -----------
    observed_pheno : Series
         Observed phenotype values
    predicted_pgs : Series
         Predicted PGS Value
    trait_name : str
         Trait name ( Used for plots and output )
    """
    # Merge data and drop missing values
    df = pd.DataFrame({'observed': observed_pheno, 'predicted': predicted_pgs})
    df = df.dropna()
    if len(df) < 5:
        logger.warning("Insufficient sample size ,  Insufficient data for model validation")
        return
    # Calculate the correlation coefficient
    corr, p_value = stats.pearsonr(df['observed'], df['predicted'])
    logger.info(f"{trait_name} PGS Model validation  -  Correlation coefficient : {corr:.3f} (P = {p_value:.2e})")
    # Plot validation results
    plt.figure(figsize=(8, 6))
    plt.scatter(df['predicted'], df['observed'], alpha=0.6, s=30)
    # Add a trend line
    z = np.polyfit(df['predicted'], df['observed'], 1)
    p = np.poly1d(z)
    plt.plot(df['predicted'], p(df['predicted']), "r--", alpha=0.8)
    plt.xlabel('Predicted PGS')
    plt.ylabel('Observed Copy Number')
    plt.title(f'PGS Model Validation for {trait_name}\n(r = {corr:.3f}, p = {p_value:.2e})')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'pgs_validation_{trait_name}.png', dpi=300, bbox_inches='tight')
    plt.close()

def compare_predicted_observed(predicted_pgs, observed_pheno, trait_name):
    """
     Compare predicted and observed values
    Parameters:
    -----------
    predicted_pgs : Series
         Predicted PGS Value
    observed_pheno : Series
         Observed phenotype values
    trait_name : str
         Trait name
    """
    # Align data
    df = pd.DataFrame({'predicted': predicted_pgs, 'observed': observed_pheno})
    df = df.dropna()
    if len(df) < 5:
        logger.warning("Insufficient sample size ,  Insufficient data for comparison")
        return
    logger.info(f"=== {trait_name}  Analysis results  ===")
    # 1. Descriptive statistics
    logger.info(f"Mean predicted value : {df['predicted'].mean():.3f} ± {df['predicted'].std():.3f}")
    logger.info(f"Mean observed value : {df['observed'].mean():.3f} ± {df['observed'].std():.3f}")
    # 2. Paired T Test
    t_stat, p_val_ttest = stats.ttest_rel(df['observed'], df['predicted'])
    logger.info(f"Paired T Test : t = {t_stat:.3f}, p = {p_val_ttest:.3e}")
    # 3. Linear regression analysis
    X = sm.add_constant(df['predicted'])
    model = sm.OLS(df['observed'], X).fit()
    intercept = model.params['const']
    slope = model.params['predicted']
    intercept_pval = model.pvalues['const']
    slope_pval = model.pvalues['predicted']
    logger.info(f"Regression model : Observed = {intercept:.8f} + {slope:.8f} * Predicted")
    logger.info(f"Intercept significance : p = {intercept_pval:.3e}")
    logger.info(f"Slope significance : p = {slope_pval:.3e}")
    # 4. Interpretation
    if np.isclose(intercept, 0, atol=0.1) and np.isclose(slope, 1, atol=0.1):
        logger.info("Conclusion :  No significant difference between observed and predicted values ,  Consistent with neutral evolution ( Ancestry-related inheritance )")
    elif intercept > 0.1 and np.isclose(slope, 1, atol=0.1):
        logger.info("Conclusion :  Significant positive deviation ,  Suggests positive selection")
    elif intercept < -0.1 and np.isclose(slope, 1, atol=0.1):
        logger.info("Conclusion :  Significant negative deviation ,  Suggests negative selection")
    else:
        logger.info("Conclusion :  Complex pattern of deviation ,  Requires further investigation")
    # 5. Plot the comparison
    plt.figure(figsize=(10, 8))
    # Scatter plot
    plt.scatter(df['predicted'], df['observed'], alpha=0.7, s=50, label='Bajau Individuals')
    # Regression line
    x_range = np.linspace(df['predicted'].min(), df['predicted'].max(), 100)
    y_pred = model.predict(sm.add_constant(x_range))
    plt.plot(x_range, y_pred, 'r-', linewidth=2, label='Regression Line')
    # y=x Line (Neutral expectation)
    plt.plot(x_range, x_range, 'k--', linewidth=2, label='y = x (Neutral Expectation)')
    plt.xlabel('Predicted Copy Number (from PGS Model)', fontsize=12)
    plt.ylabel('Observed Copy Number', fontsize=12)
    plt.title(f'{trait_name}: Observed vs. Predicted in Bajau Population\n'
              f'Intercept = {intercept:.2f} (p = {intercept_pval:.2e}), '
              f'Slope = {slope:.2f} (p = {slope_pval:.2e})', fontsize=14)
    plt.legend(fontsize=11)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'bajau_comparison_{trait_name}.png', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main()
