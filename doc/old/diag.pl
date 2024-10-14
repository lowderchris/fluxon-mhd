#!/usr/bin/perl
use strict;
use warnings;
use File::Find;

# Function to search for the autoload directory
sub find_autoload_dir {
    my @search_paths = @_;
    my $autoload_dir;

    find(sub {
        if ($_ eq 'Flux' && $File::Find::dir =~ m!/auto/share/dist$!) {
            $autoload_dir = $File::Find::dir;
            $autoload_dir =~ s!/share/dist$!/share/dist/Flux!;
            $File::Find::prune = 1;  # Stop searching once found
        }
    }, @search_paths);

    return $autoload_dir;
}

# Define the search paths (e.g., @INC and PERL5LIB)
my @search_paths = (@INC, split(':', $ENV{'PERL5LIB'}));

# Search for the autoload directory
my $autoload_dir = find_autoload_dir(@search_paths);

# Print the @INC paths with each path on its own line
print "Checking Perl \@INC paths...\n";
foreach my $inc (@INC) {
    print "  @INC path: $inc\n";
}

# Check if the PERL5LIB environment variable is set and print it nicely
print "\nChecking PERL5LIB environment variable...\n";
if ($ENV{'PERL5LIB'}) {
    my @perl5lib_paths = split(':', $ENV{'PERL5LIB'});
    foreach my $path (@perl5lib_paths) {
        print "  PERL5LIB path: $path\n";
    }
} else {
    print "  PERL5LIB is not set.\n";
}

# Check if the PDLLIB environment variable is set and print it nicely
print "\nChecking PDLLIB environment variable...\n";
if ($ENV{'PDLLIB'}) {
    my @pdllib_paths = split(':', $ENV{'PDLLIB'});
    foreach my $path (@pdllib_paths) {
        print "  PDLLIB path: $path\n";
    }
} else {
    print "  PDLLIB is not set.\n";
}

# Check if Flux module can be loaded
print "\nAttempting to load Flux module...\n";
eval {
    require Flux;
    Flux->import();
    print "  Flux module loaded successfully!\n";
};
if ($@) {
    print "  Error loading Flux module: $@\n";
}

# Check if the autoload directory for Flux is accessible
print "\nChecking if autoload directory exists...\n";
if ($autoload_dir && -d $autoload_dir) {
    print "  Autoload directory exists: $autoload_dir\n";
} else {
    print "  Autoload directory does NOT exist or is inaccessible.\n";
}
