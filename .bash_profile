# ~/.bash_profile

# This file is sourced by bash for login shells.  The following line
# runs your .bashrc and is recommended by the bash info pages.

[[ -f ~/.bashrc ]] && . ~/.bashrc

###dsource ~/.profile

#terminal customizations
export PS1="\[\033[36m\]\u\[\033[m\]@\[\033[37m\]\h:\[\033[33m\]\w\[\033[m\]$ "
export CLICOLOR=1
export LSCOLORS=fxFxBxDxCxegedabagacad
#

# some settings to prefer homebrew paths in case it exists:
if which -s brew ; then
    #PATH="/usr/local/bin:/usr/local/sbin:$PATH"
    PATH="/usr/local/bin:/usr/local/sbin:$PATH"
    #PYTHONPATH="/usr/local/lib/python2.7/site-packages:$PYTHONPATH"
fi

alias showFiles='defaults write com.apple.finder AppleShowAllFiles YES; killall Finder /System/Library/CoreServices/Finder.app'

alias hideFiles='defaults write com.apple.finder AppleShowAllFiles NO; killall Finder /System/Library/CoreServices/Finder.app'

alias topp='top -o cpu -O +rsize -s 5 -n 20'

alias ls="/bin/ls -aF"

alias desktop="cd ~/Desktop"

alias cbmcc="cd ~/Documents/Academia/CBM/CCclasses/CC/"

alias hy299="ssh asd223@hy299-server.icmb.cornell.edu"

alias bscb="ssh asd223@bscb-teaching.cb.bscb.cornell.edu"

alias rsfinal="rsync -arz --stats /Users/ashleysdoane/Documents/Academia/CBM/CCclasses/CC/BTRY6831/final/ asd223@bscb-teaching.cb.bscb.cornell.edu:/home/local/CORNELL/asd223/final/"


# added by Anaconda 2.0.1 installer
export PATH="/Users/ashleysdoane/anaconda/bin:$PATH"

export NETBOX_HOME="/Users/ashleysdoane/netbox"

export PATH=$PATH:$NETBOX_HOME/bin

# mysql

#alias mysql=/usr/local/mysql/bin/mysql
#alias mysqladmin=/usr/local/mysql/bin/mysqladmin
#export PATH=/usr/local/mysql/bin:$PATH
date

export PATH=/Users/ashleysdoane/anaconda/bin:/usr/local/bin:/usr/local/sbin:/Users/ashleysdoane/libraries/samtools-1.2:/Users/ashleysdoane/anaconda/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/usr/local/git/bin:/usr/texbin:/Users/ashleysdoane/netbox/bin:/usr/local/mysql/bin

export PATH=/Users/ashleysdoane/anaconda/bin:/usr/local/bin:/usr/local/sbin:/Users/ashleysdoane/libraries/samtools-1.2:/Users/ashleysdoane/anaconda/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin:/usr/local/git/bin:/usr/texbin:/Users/ashleysdoane/netbox/bin:/usr/local/mysql/bin:/usr/local/mysql/bin
