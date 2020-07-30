FORTUNE=$(fortune -n 62 -s computers science wisdom fortunes  | tr '\r\n\t' ' ')
echo "
\|/          (__)
      \------(oo)  ${FORTUNE}
       ||    (__)
       ||w--||     \|/
   \|/
"

# ZSH configuration
export ZSH="/Users/gregorygydush/.oh-my-zsh"
ZSH_DISABLE_COMPFIX=true

plugins=(
    git
    osx
    poetry
)

source $ZSH/oh-my-zsh.sh
source ~/.zshprompt

# # Initialize pyenv
if command -v pyenv 1>/dev/null 2>&1; then
  eval "$(pyenv init -)"
  eval "$(pyenv virtualenv-init -)"
fi

# Fixing PATH
export PYENV_ROOT="$HOME/.pyenv"
export PATH="$PYENV_ROOT/bin:$PATH"
export PATH="$PATH:$HOME/bin"
export PATH="$HOME/.local/bin:$PATH"