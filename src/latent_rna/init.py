from pathlib import Path
import shutil
import pkg_resources
import yaml

def init_project(project_dir: Path):
    """Initialize a new latent-rna project directory
    
    Args:
        project_dir: Path to create new project in
    """
    # Create directory structure
    project_dir.mkdir(exist_ok=False)
    (project_dir / 'covg_norm').mkdir(exist_ok=True)
    (project_dir / 'gene_bins').mkdir(exist_ok=True) 
    (project_dir / 'models').mkdir(exist_ok=True)
    (project_dir / 'phenotypes').mkdir(exist_ok=True)

    # Copy default config
    config_template = pkg_resources.resource_string(
        'latent_rna', 'config.default.yaml'
    ).decode('utf-8')
    
    config_path = project_dir / 'config.yaml'
    with open(config_path, 'w') as f:
        f.write(config_template)
        
    print(f"""
Project initialized at {project_dir}

Directory structure created:
  {project_dir}/
    config.yaml  # Edit this file to configure your project
    covg_norm/   # Normalized coverage data
    gene_bins/   # Gene bins
    models/      # Fitted models
    phenotypes/  # Phenotype data
    """)
