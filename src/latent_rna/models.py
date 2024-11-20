import numpy as np
import pandas as pd
from skfda.representation.basis import BSplineBasis
from skfda.representation.grid import FDataGrid
from skfda.preprocessing.dim_reduction import FPCA
from sklearn.decomposition import PCA
from tqdm import tqdm
from pathlib import Path
import pickle
from typing import Iterator
from .coverage import CoverageData

class Model:
    def __init__(self, fpca: bool, fpca_x_values: str, fpca_basis: str):
        """Initialize a model object

        Args:
            fpca: Whether to fit a functional PCA model instead of regular PCA.
            fpca_x_values: Whether to use bin numbers or genomic positions as
              x-values for functional PCA. Options are 'bin' and 'pos'.
            fpca_basis: Basis function to use for functional PCA. Options are
              'discrete' for discretized FPCA directly on the data, and 'spline'
              for a 4th-order B-spline basis.
        """
        self.model = None
        self.features = None
        self.fpca = fpca
        self.fpca_x_values = fpca_x_values
        self.fpca_basis = fpca_basis

    def fit(self, df: pd.DataFrame, var_expl: float, n_pcs_max: int):
        """Fit functional PCA model to normalized coverage for one gene
        
        Saves enough PCs to explain var_expl variance or n_pcs_max, whichever
        is fewer.

        Args:
            df: The input DataFrame containing normalized coverage per bin per
              sample.
            var_expl: Maximum variance explained by the PCs kept per gene. Pass
              0 or 1 for no variance explained cutoff.
            n_pcs_max: Maximum number of PCs to keep per gene. Pass 0 for no
              cutoff.
        """
        if self.fpca:
            n_bins_with_var = (df.std(axis=1) > 0).sum()
            x = df.values.T
            if self.fpca_x_values == 'bin':
                fd = FDataGrid(x, grid_points=np.arange(x.shape[1]))
            else:
                fd = FDataGrid(x, grid_points=df.index.get_level_values('pos'))
            # FPCA currently defaults to n_components=3, so we need to set it explicitly:
            # (There's no option for variance explained cutoff, so we'll subset the PCs later.)
            if n_pcs_max == 0:
                n_pcs_max = x.shape[0]
            n_comp = min(n_pcs_max, x.shape[0], n_bins_with_var)
            if n_comp == 0:
                return
            if self.fpca_basis == 'spline':
                self.SPLINE_ORDER = 4
                if n_comp < self.SPLINE_ORDER:
                    return
                fd = fd.to_basis(BSplineBasis(n_basis=n_comp, order=self.SPLINE_ORDER))
            # Must pass weights to avoid 'Matrix is not positive definite' error:
            model = FPCA(n_components=n_comp, _weights=np.ones(x.shape[1]))
            try:
                model.fit(fd)
            except Exception as e:
                print(f"Error occurred during model fitting: {e}", flush=True)
                return
            if var_expl not in {0, 1}:
                subset_pcs_fpca(model, var_expl, self.fpca_basis)
        else:
            df = df.loc[df.std(axis=1) > 0, :]
            if df.shape[0] == 0:
                return
            x = df.values.T
            # Center and scale the data:
            xmean = x.mean(axis=0)
            xstd = x.std(axis=0)
            x = (x - xmean) / xstd
            model = PCA(n_components=var_expl) if var_expl not in {0, 1} else PCA()
            try:
                model.fit(x)
            except np.linalg.LinAlgError:
                print("SVD did not converge, retrying with added noise", flush=True)
                x += np.random.normal(0, 1e-6, size=x.shape)
                model.fit(x)
            subset_pcs_pca(model, n_pcs_max)
            self.features = pd.DataFrame(index=df.index) # Save feature stats for transforming new data
            self.features['mean'] = xmean
            self.features['std'] = xstd
        self.model = model

    def transform(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply functional PCA model to normalized coverage for one gene
        
        Args:
            df: The input DataFrame containing normalized coverage per bin per
              sample.

        Returns:
            The DataFrame of transformed data, with gene ID and PC number as
            index levels and sample IDs as columns.
        """
        # Filter df to include same features as the model:
        if self.fpca:
            x = df.values.T
            if self.fpca_x_values == 'bin':
                fd = FDataGrid(x, grid_points=np.arange(x.shape[1]))
            else:
                fd = FDataGrid(x, grid_points=df.index.get_level_values('pos'))
            if self.fpca_basis == 'spline':
                n_basis = self.model.components_.n_basis
                fd = fd.to_basis(BSplineBasis(n_basis=n_basis, order=self.SPLINE_ORDER))
            mat = self.model.transform(fd).T
        else:
            df = df.loc[self.features.index, :]
            x = df.values.T
            # Center and scale the data using stats from the model-fitting data:
            x = (x - self.features['mean'].values) / self.features['std'].values
            mat = self.model.transform(x).T
        pc_names = [f'PC{i + 1}' for i in range(mat.shape[0])]
        out = pd.DataFrame(mat, index=pc_names, columns=df.columns)
        out.index.set_names('PC', inplace=True)
        # Add gene_id to index:
        out.reset_index(inplace=True)
        out.insert(0, 'gene_id', df.index.get_level_values('gene_id')[0])
        out.set_index(['gene_id', 'PC'], inplace=True)
        return out

def subset_pcs_fpca(fpca: FPCA, var_expl: float, fpca_basis: str):
    """Subset the FPCA model using variance explained cutoff
    
    Args:
        fpca: The input FPCA model.
        var_expl: Maximum variance explained by the PCs kept per gene.
        fpca_basis: Basis function to use for functional PCA. Options are
          'discrete' for discretized FPCA directly on the data, and 'spline'
          for a 4th-order B-spline basis.
    """
    # Find the number of PCs needed to explain var_expl variance:
    n_pcs_max = np.argmax(np.cumsum(fpca.explained_variance_ratio_) >= var_expl) + 1
    if fpca_basis == 'discrete':
        if fpca.components_.shape[0] <= n_pcs_max:
            return
        fpca.components_.data_matrix = fpca.components_.data_matrix[:n_pcs_max, :]
    else:
        if fpca.components_.coefficients.shape[0] <= n_pcs_max:
            return
        fpca.components_.coefficients = fpca.components_.coefficients[:n_pcs_max, :]
    fpca.explained_variance_ = fpca.explained_variance_[:n_pcs_max]
    fpca.explained_variance_ratio_ = fpca.explained_variance_ratio_[:n_pcs_max]
    fpca.singular_values_ = fpca.singular_values_[:n_pcs_max]

def subset_pcs_pca(pca: PCA, n_pcs_max: int):
    """Subset the PCA model to n_pcs_max PCs to reduce file size
    
    Args:
        pca: The input PCA model.
        n_pcs_max: Maximum number of PCs to keep per gene. Pass 0 for no cutoff.
    """
    if n_pcs_max != 0 and n_pcs_max < pca.n_components_:
        pca.components_ = pca.components_[:n_pcs_max, :]
        pca.n_components_ = n_pcs_max
        pca.explained_variance_ = pca.explained_variance_[:n_pcs_max]
        pca.explained_variance_ratio_ = pca.explained_variance_ratio_[:n_pcs_max]
        pca.singular_values_ = pca.singular_values_[:n_pcs_max]
        pca.noise_variance_ = None # Could compute if necessary, but omitting otherwise

def fit_batch(norm_covg_dirs: list, batch: int, var_expl_max: float, n_pcs_max: int, fpca: bool = False, fpca_x_values: str = 'bin', fpca_basis: str = 'discrete') -> dict:
    """Fit functional PCA models for all genes in a batch
    
    Args:
        norm_covg_dirs: List of directories containing normalized coverage data.
        batch: ID of the batch to fit.
        var_expl_max: Maximum variance explained by the PCs kept per gene. Pass
          0 or 1 for no variance explained cutoff.
        n_pcs_max: Maximum number of PCs to keep per gene. Pass 0 for no cutoff.
        fpca: Whether to fit a functional PCA model instead of regular PCA.
        fpca_x_values: Whether to use bin numbers or genomic positions as
          x-values for functional PCA. Options are 'bin' and 'pos'.
        fpca_basis: Basis function to use for functional PCA. Options are
          'discrete' for discretized FPCA directly on the data, and 'spline' for
          a 4th-order B-spline basis.

    Returns:
        Dictionary of fitted models, with keys 'var_expl_max', 'n_pcs_max', and
        'models'. The 'models' key contains a dictionary with gene IDs as keys
        and Model objects as values.
    """
    covg = CoverageData(norm_covg_dirs, batch)
    models = {'var_expl_max': var_expl_max, 'n_pcs_max': n_pcs_max, 'models': {}}
    for gene_id, x in tqdm(covg.by_gene(), total=len(covg.genes), desc="Fitting models"):
        model = Model(fpca, fpca_x_values, fpca_basis)
        model.fit(x, var_expl_max, n_pcs_max)
        if model.model is not None:
            models['models'][gene_id] = model
    return models

def fit(norm_covg_dirs: list, batches: list, var_expl_max: float, n_pcs_max: int, output_dir: Path, fpca: bool = False, fpca_x_values: str = 'bin', fpca_basis: str = 'discrete'):
    """Fit functional PCA models for all genes in one or more batches
    
    Args:
        norm_covg_dirs: List of directories containing normalized coverage data.
        batches: List of batch IDs to fit.
        var_expl_max: Maximum variance explained by the PCs kept per gene. Pass
          0 or 1 for no variance explained cutoff.
        n_pcs_max: Maximum number of PCs to keep per gene. Pass 0 for no cutoff.
        output_dir: Directory in which to save model pickle files.
        fpca: Whether to fit a functional PCA model instead of regular PCA.
        fpca_x_values: Whether to use bin numbers or genomic positions as
          x-values for functional PCA. Options are 'bin' and 'pos'.
        fpca_basis: Basis function to use for functional PCA. Options are
          'discrete' for discretized FPCA directly on the data, and 'spline' for
          a 4th-order B-spline basis.
    """
    output_dir.mkdir(exist_ok=True)
    for batch in batches:
        print(f'=== Fitting models for batch {batch} ===', flush=True)
        models = fit_batch(norm_covg_dirs, batch, var_expl_max, n_pcs_max, fpca, fpca_x_values, fpca_basis)
        outfile = output_dir / f'models_batch_{batch}.pickle'
        with open(outfile, 'wb') as f:
            pickle.dump(models, f)
        print(f'Models saved to {outfile}', flush=True)

def transform_batch(norm_covg_dir: Path, batch: int, models_dir: dict) -> Iterator[pd.DataFrame]:
    """Apply functional PCA transformation to all genes in a batch
    
    Args:
        norm_covg_dir: Directory of per-batch numpy binary files with normalized
          coverage.
        batch: ID of the batch to transform.
        models_dir: Directory of saved models to load and use for transformation.

    Yields:
        DataFrame of transformed data for one gene, with gene ID and PC number
        as index levels and sample IDs as columns.
    """
    covg = CoverageData([norm_covg_dir], batch)
    print('  Loading models...', flush=True)
    models_file = models_dir / f'models_batch_{batch}.pickle'
    with open(models_file, 'rb') as f:
        models = pickle.load(f)
    for gene_id, x in tqdm(covg.by_gene(), total=len(covg.genes), desc="  Generating phenotypes per gene"):
        if gene_id in models['models']:
            yield models['models'][gene_id].transform(x)

def transform(norm_covg_dir: Path, models_dir: dict, n_batches: int) -> pd.DataFrame:
    """Apply functional PCA transformation to all genes
    
    Args:
        norm_covg_dir: Directory of per-batch numpy binary files with normalized
          coverage.
        models_dir: Directory of saved models to load and use for transformation.
        n_batches: Number of batches in the data. Latent phenotypes from all
          batches will be computed and concatenated.
    """
    out = []
    for batch in range(n_batches):
        print(f'Transforming batch {batch}', flush=True)
        out.extend(transform_batch(norm_covg_dir, batch, models_dir))
    out = pd.concat(out)
    out = out.reset_index()
    return out

