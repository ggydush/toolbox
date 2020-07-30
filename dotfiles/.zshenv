# Colorize man pages
export LESS_TERMCAP_mb=$(printf "\e[1;31m")
export LESS_TERMCAP_md=$(printf "\e[1;31m")
export LESS_TERMCAP_me=$(printf "\e[0m")
export LESS_TERMCAP_se=$(printf "\e[0m")
export LESS_TERMCAP_so=$(printf "\e[1;44;33m")
export LESS_TERMCAP_ue=$(printf "\e[0m")
export LESS_TERMCAP_us=$(printf "\e[1;32m")

# pyenv disable warning
export PYENV_VIRTUALENV_DISABLE_PROMPT=1

# Aliases
alias update="source ~/.zshrc"
alias cdd='cd ~/dev/'

# Text file viewer
vtxt () {
  tabview $1
}

# Change virtual envs automatically
function cd() {
  builtin cd "$@"
  if [[ -z "$VIRTUAL_ENV" ]] ; then
      if [[ -d ./.venv ]] ; then
        source ./.venv/bin/activate
      fi
  else
      parentdir="$(dirname "$VIRTUAL_ENV")"
      if [[ "$PWD"/ != "$parentdir"/* ]] ; then
        deactivate
      fi
  fi
}