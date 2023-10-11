#!/bin/zsh

# Ensure that FL_PREFIX and PL_PREFIX are defined or exit if they're not
: "${FL_PREFIX:?Variable FL_PREFIX not set}"
: "${PL_PREFIX:?Variable PL_PREFIX not set}"

# Remove FLUXpipe installation and related PERL5LIB path
echo "Uninstalling fluxpipe..."
cd "$FL_PREFIX/fluxpipe"
perl -Mlocal::lib=local/ -e 'print "$$_\n" for @INC' | xargs -I {} rm -rf {}
# Create a backup before using sed in-place editing
cp ~/.zshrc ~/.zshrc.bak
sed -i.bak '/export PERL5LIB=$FL_PREFIX\/fluxpipe\/local\/lib\/perl5:$PERL5LIB/d' ~/.zshrc

# Remove local::lib setup from .zshrc
# Create another backup for this sed operation
cp ~/.zshrc ~/.zshrc.bak2
sed -i.bak2 '/eval "$$(perl -I$PL_PREFIX\/lib\/perl5\/lib\/perl5 -Mlocal::lib=$PL_PREFIX\/lib\/perl5)"/d' ~/.zshrc

# Deactivate and remove conda environment
conda deactivate
conda env remove -n fluxenv
conda deactivate

# Remove files and directories
rm -rf "$FL_PREFIX/fluxpipe/local"
rm -rf "$FL_PREFIX/fluxpipe/fluxpipe.egg-info"
rm -rf "$FL_PREFIX/fluxpipe/blib"
rm -rf "$FL_PREFIX/fluxpipe/_Inline"
rm -rf "$FL_PREFIX/fluxpipe/__pycache__"
rm -rf "$FL_PREFIX/fluxpipe/.pytest_cache"
rm -rf "$FL_PREFIX/fluxpipe/pm_to_blib"


# Uninstall Perl modules installed via cpanm
# cd "$PL_PREFIX"
# cpanm --uninstall File::ShareDir PDL::Graphics::Gnuplot Math::RungeKutta
# ... Add other dependencies to be removed ...

# Optional sections are kept commented.
# Uncomment them if you want to use them and ensure you understand the repercussions.

# Optional: Uninstall perlbrew and its installations (Warning: this will remove ALL perlbrew installations)
# perlbrew off
# rm -rf ~/perl5/perlbrew

# Optional: Uninstall Homebrew (Warning: this will remove ALL Homebrew installations and packages)
# Uncomment the following two lines ONLY if you're sure about removing Homebrew:
# /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall.sh)"

echo "\tFluxpipe uninstall complete!\n\n"
