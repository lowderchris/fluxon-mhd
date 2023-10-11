use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME      => 'fluxpipe',
    VERSION   => '0.01',
    PREREQ_PM => {},           # e.g., { 'Some::Module' => 0.01 },
    dist      => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean     => { FILES    => 'fluxpipe-*' },

    # EXE_FILES => ['uninstaller_script.sh'],

    # Add other attributes as needed
);

sub MY::postamble {
    return <<'MAKEFILE';

everything: install_fluxpipe

install_fluxpipe:
	./install-fluxpipe.sh;

install_fp:
	./install-fluxpipe.sh;

uninstall_fluxpipe:
	./uninstall-fluxpipe.sh

uninstall_fp:
	./uninstall-fluxpipe.sh
MAKEFILE
}
