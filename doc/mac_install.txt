# Macintosh Installation

Note that for some quasi-locked down Macintosh systems where the user does not have primary access to an admin user, a slightly different installation is required, outlined here.

Here's the path to installing FLUX:

- Admin to install homebrew perl, or local perl

(Admin) brew install perl

- Admin to install homebrew packages

(Admin) brew install gnuplot

(Admin) brew install fftw

- Perldlrc

    require('PDL/default.perldlrc');

    use PDL::AutoLoader;
    $PDL::AutoLoader::Rescan=1;

    1;

Make sure to push to @INC (find from current perldlrc file)

- User directory (~/lib/) CPAN install PDL

PERL_MM_OPT="INSTALL_BASE=$HOME/lib/perl5"
cpan local::lib
echo 'eval "$(perl -I$HOME/lib/perl5/lib/perl5 -Mlocal::lib=$HOME/lib/perl5)"' >> ~/.zshrc

cpan File::ShareDir
cpan PDL::Graphics::Gnuplot
cpan Math::RungeKutta

- User directory (~/lib/) FLUX compilation

export PL_PREFIX='/Users/clowder/perl5'
export FL_PREFIX='/Users/clowder/lib'

Make everything