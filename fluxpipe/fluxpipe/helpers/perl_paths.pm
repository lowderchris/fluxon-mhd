# Function to add a directory and its subdirectories to various INC arrays
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


# Function to populate INC arrays based on the FL_PREFIX environment variable
sub make_perl_incs {
    my $root_directory = $ENV{'FL_PREFIX'};
    my @exclusions = ("/t", ".git", "_Inline", "sandbox", "blib", "_pycache__");

    # Check if the environment variable is set
    if (defined $root_directory) {
        add_to_incs($root_directory, @exclusions);
    } else {
        warn "Environment variable FL_PREFIX is not set. Skipping.\n";
    }
    return 1;
}


# Function to set Perl environment variables
sub fix_envs {
    print("Fixing Envs!");
    my $perl_path = $ENV{'PL_PREFIX'};
    $ENV{PERL_LOCAL_LIB_ROOT} = $perl_path;
    $ENV{PERL5LIB} = $ENV{PERL5LIB}.": $perl_path/lib/perl5";
    return 1;
}

########################################################################
# Main Code
# ----------------------------------------------------------------------
fix_envs();
make_perl_incs();
print "->Added flux directories to INC arrays<-\n";

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
