use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME      => 'fluxpipe',
    VERSION   => '0.01',
    PREREQ_PM => {
        # 'File::ShareDir',
        # 'PDL::Graphics::Gnuplot',
        # 'Math::RungeKutta',

        # 'Some::Module'    => 0.01,
        # 'Another::Module' => 1.42,
        # Add your dependencies here
    },
    ABSTRACT => 'Fluxpipe - a pipeline for fluxon analysis',
    AUTHOR   => 'Gilly <gilly@swri.org>',
    LICENSE  => 'perl',
    INST_LIB => "$ENV{'FL_PREFIX'}",
    dist     => { COMPRESS => 'gzip -9f', SUFFIX => 'gz' },
    clean    => { FILES    => 'fluxpipe-*' },
);

sub MY::postamble {

    # Define custom targets without creating a circular dependency
    return <<'MAKEFILE';
.PHONY: install_fluxpipe uninstall_fluxpipe

# Define the target "all" as a dependency of other targets using ::
all :: PREFIX_PATHS.sh install_fluxpipe

# This target runs the PREFIX_PATHS.sh script to set FL_PREFIX
PREFIX_PATHS.sh ::
	./PREFIX_PATHS.sh

# Combine both Linux and macOS install into a single install_fluxpipe target
install_fluxpipe ::
	@if [ `uname -s` = "Linux" ]; then \
		./install-fluxpipe-ubuntu.sh; \
    elif [ `uname -s` = "Darwin" ]; then \
		./install-fluxpipe-macos.sh; \
    else \
		echo "Unsupported operating system"; \
	fi

uninstall_fluxpipe ::
	@./uninstall-fluxpipe.sh
# @echo "Uninstalling fluxpipe..."

MAKEFILE
}
