# Install jupyter
pip install jupyterlab

# Add horizon theme
jupyter labextension install @mohirio/jupyterlab-horizon-theme

# Add Black formatting (use CTRL + B to format)
pip install black
jupyter labextension install @ryantam626/jupyterlab_code_formatter
pip install jupyterlab_code_formatter
jupyter serverextension enable --py jupyterlab_code_formatter
