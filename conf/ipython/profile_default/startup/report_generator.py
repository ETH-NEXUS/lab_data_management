import argparse
import sh
import papermill as pm
from ipylab import JupyterFrontEnd
import os

import os


def generate_report(notebook_path, experiment, label):
    quarto = sh.Command("quarto").bake("render")

    input_dir = os.path.dirname(os.path.dirname(notebook_path))
    output_dir = os.path.join(input_dir, 'output', experiment)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    base_filename = os.path.basename(notebook_path).replace('.ipynb', '')
    output_filename = f"{base_filename}_{experiment}_output.ipynb"
    output_path = os.path.join(output_dir, output_filename)

    parameters = {"experiment": experiment, "label": label}
    pm.execute_notebook(notebook_path, output_path, parameters=parameters or {})

    _args = ["--to", "pdf", "-M", "echo:false", "--no-execute"]
    quarto(output_path, *_args)

    pdf_path = output_path.replace('.ipynb', '.pdf')
    app = JupyterFrontEnd()
    app.commands.execute("docmanager:open", {"path": pdf_path})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate PDF report.')
    parser.add_argument('--notebook_path', type=str, help='Path to the notebook')
    parser.add_argument('--experiment', type=str, help='Experiment name')
    parser.add_argument('--label', type=str, help='Label name')

    args = parser.parse_args()
    generate_report(args.notebook_path, args.experiment, args.label)
