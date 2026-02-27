from pathlib import Path
import shutil
from importlib.resources import files
import yaml

def copy_resource(resource_name: str, outfile: Path):
    """Copy a resource file to a given output file"""
    content = files('laddr.resources').joinpath(resource_name).read_text(encoding='utf-8')
    with open(outfile, 'w', encoding='utf-8') as f:
        f.write(content)

def init_project(project_dir: Path, config_type: str = 'default', template: str = 'both'):
    """Initialize a new LaDDR project directory
    
    Args:
        project_dir: Path to create new project in.
        config_type: Type of config file to initialize project with. "default"
          includes parameters for the recommended binning and model fitting
          methods. "example" includes a config and coverage manifest file set up
          to run the example dataset. "extended" includes all possible config
          parameters.
        template: Specify "snakemake" to include a snakefile, "shell" to include
          a shell script with the basic commands for running the project, or
          "both" to include both.
    """
    # Create directory structure
    if project_dir.exists():
        raise FileExistsError(f"Directory {project_dir} already exists")
    project_dir.mkdir()

    copy_resource(f'config.{config_type}.yaml', project_dir / 'config.yaml')
    message = f"""Project initialized{" for example data" if config_type == "example" else ""} at {project_dir}

Files created:
  {project_dir}/
    config.yaml: """
    if config_type in ['default', 'extended']:
        message += "Edit this file to configure your project\n"
    else:
        message += "Parameters and input data paths for the example dataset\n"

    if config_type == 'example':
        copy_resource('coverage_manifest.tsv', project_dir / 'coverage_manifest.tsv')
        message += f"    coverage_manifest.tsv: Table specifying example datasets, samples, and coverage files\n"

    if template in ['both', 'shell']:
        copy_resource('run.sh', project_dir / 'run.sh')
        message += f"    run.sh: Shell script for running the project\n"
    if template in ['both', 'snakemake']:
        copy_resource('Snakefile', project_dir / 'Snakefile')
        message += f"    Snakefile: Snakemake file for running the project\n"

    print(message, end='')
