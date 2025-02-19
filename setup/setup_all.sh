# Install Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

# # Install GUI apps
brew cask install calendar-366
brew cask install evernote
brew cask install google-backup-and-sync
brew cask install google-chrome
brew cask install google-drive-file-stream
brew cask install iterm2
brew cask install kite
brew cask install hyperdock
brew cask install rocket
brew cask install slack
brew cask install spotify
brew cask install vanilla
brew cask install visual-studio-code
brew cask install zoomus

# Install command line apps
brew install fortune
brew install node
brew install pyenv
brew install zsh

# Install ohmyzsh
sh -c "$(curl -fsSL https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh)"

# Setup environment
SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
ln -s $(dirname $SCRIPTPATH)/dotfiles/.zsh* ~/

# Setup python
sh setup_python.sh

# Poetry completions
mkdir $ZSH/plugins/poetry
poetry completions zsh > $ZSH/plugins/poetry/_poetry

echo "
There are some things I cannot automate...

Sync Vscode settings
https://marketplace.visualstudio.com/items?itemName=Shan.code-settings-sync

Install latest Android messages desktop
https://github.com/chrisknepper/android-messages-desktop/releases/

Install Microsoft Todo
https://todo.microsoft.com/tasks/

Install Alfred and change settings to Google Drive location

Add iterm settings profile in ~/dev/toolbox/dotfiles/

Login to all other Chrome, Evernote, Todo, Messages, etc. to sync settings
"