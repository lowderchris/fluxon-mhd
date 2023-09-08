#!/usr/bin/perl
# use strict;
use warnings;
use File::Find;

=head1 NAME

perl+paths - Allows Perl scripts to find the fluxon-mhd libraries.

=head1 SYNOPSIS

    use perl_paths.pm;

    if (defined $ENV{'FL_PREFIX'}) {
    my $envpath = "$ENV{'FL_PREFIX'}/flux_pipe/perl_paths.pm";
    print "\nAttempting to require: $envpath\n";

    # Check if the file exists and is readable
    if (-e $envpath && -r _) {
        require $envpath;
        print "\n\n";
    } else {
        warn "File does not exist or is not readable: $envpath\n\n";
    }
    } else {
        warn "Environment variable FL_PREFIX is not set.\n\n";
    }

=head1 DESCRIPTION

This file gives Perl scripts access to the fluxon-mhd libraries.

=cut

=head2 add_to_incs

    add_to_incs($dir, @exclusions);

Adds a directory and its subdirectories to the @INC, @PDLLIB, and @PDL_INC arrays.

=over 4

=item $dir

The directory to add.

=item @exclusions

An array of directory names to exclude.

=back

=cut

sub add_to_incs {
    my ($dir, @exclusions) = @_;
    find(sub {
        if (-d $File::Find::name) {
            # Skip excluded directories
            foreach my $exclusion (@exclusions) {
                return if $File::Find::name =~ /$exclusion/;
            }
            push @INC, $File::Find::name;
            push @PDLLIB, $File::Find::name;
            push @PDL_INC, $File::Find::name;
        }
    }, $dir);
    return 1;
}

=head2 make_perl_incs

    make_perl_incs();

Populates the INC arrays based on the FL_PREFIX environment variable.

=cut

sub make_perl_incs {
    my $root_directory = $ENV{'FL_PREFIX'};
    my @exclusions = ('/t', '.git', '_Inline', 'sandbox', 'blib', '_pycache__');

    # Check if the environment variable is set
    if (defined $root_directory) {
        add_to_incs($root_directory, @exclusions);
    } else {
        warn 'Environment variable FL_PREFIX is not set. Skipping.\n';
    }
    return 1;
}

=head2 fix_envs

    fix_envs();

Sets the Perl environment variables based on the PL_PREFIX environment variable.

=cut

sub fix_envs {
    my $perl_path = $ENV{'PL_PREFIX'};
    $ENV{PERL_LOCAL_LIB_ROOT} = $perl_path;
    $ENV{PERL5LIB} = "$perl_path/lib/perl5";
    return 1;
}

########################################################################
# Main Code
# ----------------------------------------------------------------------

fix_envs();

make_perl_incs();

print "\n\n->Added flux directories to INC arrays<-";

my $do_print = 0;
if ($do_print) {
    print "\n\nINC has:\n ";
    print map { " $_\n" } @INC;
    print "--------------------------\n";

    print "\nPDL_INC has:\n ";
    print map { " $_\n" } @PDL_INC;
    print "--------------------------\n";

    print "\nPDLLIB has:\n ";
    print map { " $_\n" } @PDLLIB;
    print "--------------------------\n\n";
}

1;
