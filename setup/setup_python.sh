# # Install python
# pyenv install 3.8.3
# pyenv global 3.8.3

# # Install requirements
# pip install -r python_base_requirements.txt

# Add code formatting to jupyterlab
jupyter labextension install @ryantam626/jupyterlab_code_formatter
pip install jupyterlab_code_formatter
jupyter serverextension enable --py jupyterlab_code_formatter

# Add theme to jupyterlab
jupyter labextension install @mohirio/jupyterlab-horizon-theme