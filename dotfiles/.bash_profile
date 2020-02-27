export PATH
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

# Aliases
alias ll="ls -lrh"
alias jup='jupyter notebook --ip=$(hostname) --port=8080'
alias jl='jupyter-lab --ip=$(hostname) --port=8080'
alias svg='use Graphviz; snakemake --dag | dot -Tsvg > workflow.svg'

# Prompt
parse_git_branch() {
  git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/ (\1)/'
}

green=$(tput setaf 30);
white=$(tput setaf 15);
yellow=$(tput setaf 222);

PS1="";
PS1+="\[${green}\]\u"; # user
PS1+="\[${green}\] [\@"]; #timestamp
PS1+="\[${white}\] in ";
PS1+="\[${yellow}\]\w"; # working directory
PS1+="\[${white}\]\$(parse_git_branch)"
PS1+="\[${yellow}\]\n➞  "
PS1+="\[${white}\]"; # reset color
export PS1;

# Add space between commands
PROMPT_COMMAND="echo;$PROMPT_COMMAND"