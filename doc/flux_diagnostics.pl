use strict;
use warnings;

# ANSI color codes
my $green = "\e[32m";
my $red = "\e[31m";
my $reset = "\e[0m";

# Print a header for better readability
print "\n" . "=" x 80 . "\n";
print "Perl Flux Diagnostics\n";
print "=" x 80 . "\n\n";

# Print environment variable values set during Flux installation
print "Environment Variables Set During Flux Installation:\n";
print "  PL_PREFIX          : " . ($ENV{'PL_PREFIX'} ? "$ENV{'PL_PREFIX'}\n" : "${red}PL_PREFIX is unset.$reset\n");
print "  FL_PREFIX          : " . ($ENV{'FL_PREFIX'} ? "$ENV{'FL_PREFIX'}\n" : "${red}FL_PREFIX is unset.$reset\n");
print "  FL_MHDLIB          : " . ($ENV{'FL_MHDLIB'} ? "$ENV{'FL_MHDLIB'}\n" : "${red}FL_MHDLIB is unset.$reset\n");
print "  DESIRED_PERL_VERSION: " . ($ENV{'DESIRED_PERL_VERSION'} ? "$ENV{'DESIRED_PERL_VERSION'}\n" : "${red}DESIRED_PERL_VERSION is unset.$reset\n");
print "  PDLLIB             : " . ($ENV{'PDLLIB'} ? "$ENV{'PDLLIB'}\n" : "${red}PDLLIB is unset.$reset\n");
print "\n" . "=" x 80 . "\n";

# Check the current working directory
my $cwd = `pwd`;
chomp($cwd);
print "Checking current working directory...\n";
print "  Current working directory: $cwd\n";
print "\n" . "=" x 80 . "\n";

# Check @INC paths
print "Checking Perl \@INC paths...\n";
foreach my $inc (@INC) {
    print "  INC path: $inc\n";
}
print "\n" . "=" x 80 . "\n";

# Check the PERL5LIB environment variable
print "Checking PERL5LIB environment variable...\n";
if ($ENV{'PERL5LIB'}) {
    my @perl5lib_paths = split(':', $ENV{'PERL5LIB'});
    foreach my $path (@perl5lib_paths) {
        print "  PERL5LIB path: $path\n";
    }
} else {
    print "  ${red}PERL5LIB is not set.$reset\n";
}
print "\n" . "=" x 80 . "\n";

# Check the PDLLIB environment variable
print "Checking PDLLIB environment variable...\n";
if ($ENV{'PDLLIB'}) {
    my @pdllib_paths = split(':', $ENV{'PDLLIB'});
    foreach my $path (@pdllib_paths) {
        print "  PDLLIB path: $path\n";
    }
} else {
    print "  ${red}PDLLIB is not set.$reset\n";
}
print "\n" . "=" x 80 . "\n";

# Check if the autoload directory exists
print "Checking if autoload directory exists...\n";

# Get the PDLLIB environment variable
my $pdl_lib = $ENV{'PDLLIB'};

# If PDLLIB is set, split it into individual directories and search for Flux autoload dir
if ($pdl_lib) {
    my @dirs = split /:/, $pdl_lib;
    my $autoload_dir = "";

    # Search for the correct Flux autoload directory
    foreach my $dir (@dirs) {
        if ($dir =~ /Flux/) {  # Adjust this condition if needed
            $autoload_dir = $dir;
            last;
        }
    }

    # Check if the autoload directory was found and exists
    if ($autoload_dir && -d $autoload_dir) {
        print "  Autoload directory found: $autoload_dir\n";
    } else {
        print "  ${red}Autoload directory does NOT exist or was not found in PDLLIB.$reset\n";
    }
} else {
    print "  ${red}PDLLIB environment variable is not set. Please ensure the autoload directory is properly configured.$reset\n";
}
print "\n" . "=" x 80 . "\n";

# Attempt to load Flux module
print "Attempting to load Flux module...\n";
eval {
    require Flux;
    Flux->import();
    print "  ${green}Flux module loaded successfully!$reset\n";
    print "  Flux module path: " . $INC{"Flux.pm"} . "\n";
};
if ($@) {
    print "  ${red}Error loading Flux module: $@$reset\n";
}

# Final success or failure message
if (!$@) {
    print "\n" . "=" x 80 . "\n";
    print "${green}SUCCESS: Flux Diagnostics Passed!$reset\n";
} else {
    print "\n" . "=" x 80 . "\n";
    print "${red}FAILURE: Flux Diagnostics Encountered an Error.$reset\n";
}
print "=" x 80 . "\n";