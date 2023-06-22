from ipylab import JupyterFrontEnd
from ipylab.commands import CommandPalette
from os.path import split, join, sep


def generateReport(*args, **kwargs):
    import ipynbname
    import sh

    quarto = sh.Command("quarto").bake("render")
    nb = str(ipynbname.path())
    _, nb_path = split(nb)
    nb_path = join(sep, nb_path).replace("ipynb", "pdf")

    _args = ["--to", "pdf", "-M", "echo:false", "--no-execute"]
    quarto(nb, *_args)

    app = JupyterFrontEnd()
    app.commands.execute("docmanager:open", {"path": nb_path})


app = JupyterFrontEnd()
app.commands.add_command(
    "generate_report", execute=generateReport, label="Generate PDF Report"
)
palette = CommandPalette()
palette.add_item("generate_report", "_LDM Commands")
