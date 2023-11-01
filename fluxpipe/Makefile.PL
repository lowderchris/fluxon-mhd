use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME      => 'fluxpipe',
    VERSION   => '0.01',
    PREREQ_PM => {

        # 'Some::Module'    => 0.01,
        # 'Another::Module' => 1.42,
        # Add your dependencies here
    },
    ABSTRACT => 'Fluxpipe - a pipeline for flux analysis',
    AUTHOR   => 'Gilly <gilly@swri.org>',
    LICENSE  => 'perl',
    dist     => { COMPRESS => 'gzip -9f', SUFFIX => 'gz' },
    clean    => { FILES    => 'fluxpipe-*' },

  # EXE_FILES => ['scripts/my_script.pl'],    # Add your executable scripts here
);

sub MY::postamble {

    # Define custom targets without creating a circular dependency
    return <<'MAKEFILE';
.PHONY: install_fluxpipe uninstall_fluxpipe

all :: install_fluxpipe

install_fluxpipe:
	@echo "Installing fluxpipe..."
	@./install-fluxpipe.sh

uninstall_fluxpipe:
	@echo "Uninstalling fluxpipe..."
	@./uninstall-fluxpipe.sh

MAKEFILE
}
