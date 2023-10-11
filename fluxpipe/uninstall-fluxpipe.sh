#!/bin/zsh

# Remove FLUXpipe installation and related PERL5LIB path
echo "Uninstalling fluxpipe!"
cd $FL_PREFIX/fluxpipe
perl -Mlocal::lib=local/ -e 'print "$$_\n" for @INC' | xargs -I {} rm -rf {}
sed -i '/export PERL5LIB=$FL_PREFIX\/fluxpipe\/local\/lib\/perl5:$PERL5LIB/d' ~/.zshrc

# Remove local::lib setup from .zshrc
sed -i '/eval "$$(perl -I$PL_PREFIX\/lib\/perl5\/lib\/perl5 -Mlocal::lib=$PL_PREFIX\/lib\/perl5)"/d' ~/.zshrc

# Deactivate and remove conda environment
conda deactivate
conda env remove -n fluxenv

# Uninstall Perl modules installed via cpanm
#cd $(PL_PREFIX)
#cpanm --uninstall File::ShareDir PDL::Graphics::Gnuplot Math::RungeKutta
# ... Add other dependencies to be removed ...

# Optional: Uninstall perlbrew and its installations (Warning: this will remove ALL perlbrew installations)
# perlbrew off
# rm -rf ~/perl5/perlbrew

# Optional: Uninstall Homebrew (Warning: this will remove ALL Homebrew installations and packages)
# /bin/bash -c "$$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall.sh)"

echo "Uninstallation complete!"