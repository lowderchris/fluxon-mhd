
=head1 NAME

pipe_helper - Utility Functions for File and Environment Management

=head1 SYNOPSIS

    use pipe_helper;
    shorten_path($string);
    find_highest_numbered_file($directory);
    set_env_variable($variable, $value);
    # ... and so on

=head1 DESCRIPTION

This Perl script provides utility functions for managing files and environment variables.
It includes functions for shortening file paths, finding the highest-numbered file in a directory,
setting and checking environment variables, and more.

=head1 FUNCTIONS

=head2 configurations

    configurations($debug, $config_name, $config_filename);

Reads and processes a configuration file, returning a hash of the configuration settings.

=head2 shorten_path

    shorten_path($string);

Shortens the given file path by replacing the DATAPATH environment variable.

=head2 find_highest_numbered_file

    find_highest_numbered_file($directory);

Finds the highest-numbered file in the given directory.

=head2 set_env_variable

    set_env_variable($variable, $value);

Sets an environment variable to a given value.

=head2 get_env_variable

    get_env_variable($variable);

Gets the value of an environment variable.

=head2 check_env_variable

    check_env_variable($variable, $print);

Checks if an environment variable is set and optionally prints its value.

=head2 set_and_check_env_variable

    set_and_check_env_variable($variable, $value, $print);

Sets an environment variable and then checks if it is set.

=head2 calculate_directories

    calculate_directories($config_ref);

Calculates various directories based on the base directory and batch name.

=head2 set_python_path

    set_python_path($pythonpath, $print);

Sets the PYTHONPATH environment variable.

=head2 print_banner

    print_banner($batch_name, $CR, $reduction, $n_fluxons_wanted, $recompute_string);

Prints a banner with various details.

=head2 search_files_in_directory

    search_files_in_directory($directory, $known_string, $extension);

Searches for files in a directory that match a known string and file extension.

=head2 check_second_file_presence

    check_second_file_presence($file_path);

Checks for the presence of a second file related to the given file path.

=head1 AUTHOR

Gilly <gilly@swri.org> (and others!)

=head1 SEE ALSO

L<PDL::AutoLoader>, L<PDL>, L<Time::Piece>

=cut

package pipe_helper;
use strict;
use warnings;
use Exporter qw(import);
use Flux::World    qw(read_world);
our @EXPORT_OK =
  qw(shorten_path find_highest_numbered_file find_highest_numbered_file_with_string set_env_variable get_env_variable check_env_variable configs_update_magdir
  set_and_check_env_variable calculate_directories set_python_path set_paths print_banner search_files_in_directory check_second_file_presence configurations load_highest_numbered_world);
use File::Basename qw(dirname basename);
use PDL::AutoLoader;
use PDL;
use Time::Piece;
no warnings 'redefine';

=head2 shorten_path

Shortens the given file path by replacing the DATAPATH environment variable.

=cut


sub pdl_range {
    my ($start, $stop, $step) = @_;
    $step //= 1;  # Default step is 1 if not provided

    # Calculate the number of elements needed
    my $size = int(($stop - $start) / $step);
    $size++ if (($stop - $start) % $step) > 0;

    # Create the sequence
    return $start + $step * sequence($size);
}


sub configurations {
    my ( $adapt, $debug, $config_name, $config_filename ) = @_;
    $adapt           //= 0;
    $config_name     //= "DEFAULT";
    $config_filename //= "config.ini";
    $debug           //= 0;
    $debug           //= 0;

    use Config::IniFiles;
    use Cwd;
    use File::Spec::Functions;
    use File::Find;
    use File::Temp qw/ tempfile tempdir /;

    # Define the path of the configuration file
    my $config_path =
      catfile( $ENV{'FL_MHDLIB'}, "fluxpipe", $config_filename );

    # Check if the file exists at the defined path
    unless ( -e $config_path ) {
        my $found = 0;
        find(

            sub {
                if ( $_ eq $config_filename ) {
                    $config_path = $File::Find::name;
                    $found       = 1;
                }
            },

            getcwd()

        );
        die "Configuration file not found." unless $found;
    }

    # Create a temporary file to store the config data without comments
    my ( $fh, $temp_filename ) = tempfile();

    # Remove inline comments and write to temporary file
    open( my $in, '<', $config_path )
      or die "Could not open '$config_path' for reading: $!";
    while (<$in>) {
        s/#.*$//;            # Remove inline comments
        s/\s+$//;            # Remove trailing whitespace
        print $fh "$_\n";    # Append a newline character, then print to file
    }
    close $in;
    close $fh;

    # Read the configuration file
    my $cfg = Config::IniFiles->new( -file => $temp_filename );

    # Select the correct section
    $config_name = $cfg->val( 'DEFAULT', 'config_name' )
      if ( $config_name eq 'DEFAULT' );

    # Load all parameters from the DEFAULT section
    my %the_config = ();
    for my $key ( $cfg->Parameters('DEFAULT') ) {
        $the_config{$key} = $cfg->val( 'DEFAULT', $key );
    }

# If a different section is specified, load its parameters, overwriting defaults where applicable
    if ( $config_name ne 'DEFAULT' ) {
        $config_name = $cfg->val( 'DEFAULT', 'config_name' )
          if ( $config_name eq 'DEFAULT' );
        for my $key ( $cfg->Parameters($config_name) ) {
            $the_config{$key} = $cfg->val( $config_name, $key );
        }
    }

    # Perform additional processing on the configuration settings

    $the_config{'adapt'}       = $adapt;
    $the_config{'abs_rc_path'} = glob( $the_config{'rc_path'} );
    $the_config{"run_script"} =
      catfile( $the_config{'fl_mhdlib'}, $the_config{"run_script"} );

    #if the first character of the rotations is a [ then it is a list of rotations
    if ( substr( $the_config{'rotations'}, 0, 1 ) eq "[" ) {
        $the_config{'rotations'}    =~ s/[\[\]]//g;
        $the_config{'rotations'} = PDL->new( split( /\s*,\s*/, $the_config{'rotations'} ) );
    } else {
        #if the first character of the rotations is a ( then it is a start, stop, step
        if ( substr( $the_config{'rotations'}, 0, 1 ) eq "(" ) {
            $the_config{'rotations'} =~ s/[\(\)]//g;
            my ($start, $stop, $step) = split( /\s*,\s*/, $the_config{'rotations'} );
            # my $pdl = pdl_range(10, 20, 2);
            $the_config{'rotations'} = pdl_range($start, $stop, $step);
        } else {
            #otherwise it is a single rotation
            $the_config{'rotations'} = PDL->new( $the_config{'rotations'} );
        }
    }

    if ( substr( $the_config{'flow_method'}, 0, 1 ) eq "[" ) {
        $the_config{'flow_method'}    =~ s/[\[\]]//g;
        $the_config{'flow_method'} = [split( /\s*,\s*/, $the_config{'flow_method'} )] ;
    } else {
        $the_config{'flow_method'} = [$the_config{'flow_method'}];
    }

    # Remove brackets from rotations and fluxon_count
    $the_config{'fluxon_count'} =~ s/[\[\]]//g;
    $the_config{'adapts'}       =~ s/[\[\]]//g;

    # Create PDL objects

    $the_config{'fluxon_count'} =
      PDL->new( split( /\s*,\s*/, $the_config{'fluxon_count'} ) );
    $the_config{'adapts'} =
      PDL->new( split( /\s*,\s*/, $the_config{'adapts'} ) );

    $the_config{'n_jobs'} =
      $the_config{'rotations'}->nelem *
      $the_config{'fluxon_count'}->nelem *
      $the_config{'adapts'}->nelem;

    # Calculate directories
    calculate_directories( \%the_config );

    my $magfile;
    my $reduction    = $the_config{'mag_reduce'};
    my $magdir       = $the_config{'mag_dir'};
    my $adapt_select = $the_config{'adapt_select'};

    if ( $the_config{'adapt'} ) {
        $magfile = "CR%s\_rf$adapt_select\_adapt.fits";
    }
    else {
        $magfile = "CR%s\_r$reduction\_hmi.fits";
    }

    $the_config{'magfile'} = $magfile;
    $the_config{'magpath'} = "$magdir/$magfile";

    # $the_config{'flocfile'} = $flocfile;
    # $the_config{'flocpath'} = "$flocdir/$flocfile";

    if ($debug) {

        #Print the content of the configuration hash for debugging.
        print "Configuration file values:\n--------------------------------\n";
        foreach my $key ( keys %the_config ) {
            print "$key: $the_config{$key}\n";
        }
        print "--------------------------------\n\n";
    }
    return %the_config;
}

sub shorten_path_real {
    my ($string) = @_;
    my $datapath = $ENV{'DATAPATH'};
    if ($datapath) {
        $string =~ s/\Q$datapath\E/\$DATAPATH/g;
    }
    return $string;
}

sub shorten_path {
    my ($string) = @_;
    return $string;
}

sub configs_update_magdir {
    my ($configs_ref) = @_;    # get the hash reference

    my $magdir           = $configs_ref->{'mag_dir'};
    my $adapt_select     = $configs_ref->{'adapt_select'};
    my $CR               = $configs_ref->{'CR'};
    my $batchdir         = $configs_ref->{'batch_dir'};
    my $flocdir          = "$batchdir/data/cr$CR/floc";
    my $n_fluxons_wanted = $configs_ref->{'n_fluxons_wanted'};
    my $reduction        = $configs_ref->{'mag_reduce'};
    my $magfile;
    my $flocfile;

    if ( $configs_ref->{'adapt'} ) {

        # Run the adapt maps
        $magfile = "CR$CR\_rf$adapt_select\_adapt.fits";
        $flocfile =
          "floc_cr$CR\_rf$adapt_select\_f$n_fluxons_wanted\_adapt.dat";
    }
    else {
        # Run the hmi maps
        $magfile  = "CR$CR\_r$reduction\_hmi.fits";
        $flocfile = "floc_cr$CR\_r$reduction\_f$n_fluxons_wanted\_hmi.dat";
    }

    my $magpath  = "$magdir/$magfile";
    my $flocpath = "$flocdir/$flocfile";

    # $configs_ref->{'magfile'}  = $magfile;
    $configs_ref->{'magpath'}  = $magpath;
    $configs_ref->{'flocdir'}  = $flocdir;
    $configs_ref->{'flocfile'} = $flocfile;
    $configs_ref->{'flocpath'} = $flocpath;
}

=head2 find_highest_numbered_file

Finds the highest-numbered file in the given directory.

=cut

sub find_highest_numbered_file {
    my ($directory) = @_;
    opendir( my $dir_handle, $directory ) or die "Cannot open directory: $!";
    my @files = grep { !/^\.{1,2}$/ } readdir($dir_handle);
    closedir $dir_handle;

    my $highest_numbered_file;
    my $highest_number = -1;
    my $found          = 0;
    for my $file_name (@files) {

        # Match file names containing "_relaxed", a number, and ".flux"
        if ( $file_name =~ /_relaxed_s(\d+)\.flux/ ) {
            my $number = $1;
            if ( $number > $highest_number ) {
                $highest_number        = $number;
                $highest_numbered_file = $file_name;
            }
            $found = 1;
        }
    }
    return $highest_numbered_file ? "$directory$highest_numbered_file" : 0,
      $highest_number;
}

sub find_highest_numbered_file_with_string {
    my ($directory, $search_string) = @_;
    opendir(my $dir_handle, $directory) or die "Cannot open directory $directory: $!";
    my @files = grep { !/^\.{1,2}$/ } readdir($dir_handle);
    closedir $dir_handle;

    my $highest_numbered_file;
    my $highest_number = -1;
    my $found = 0;
    for my $file_name (@files) {

        # Match file names containing the search string, a number, and ".flux"
        if ( $file_name =~ /$search_string(\d+)\.flux/ ) {
            my $number = $1;
            if ( $number > $highest_number ) {
                $highest_number        = $number;
                $highest_numbered_file = $file_name;
            }
            $found = 1;
        }
    }
    return $highest_numbered_file ? "$directory/$highest_numbered_file" : 0,
      $highest_number;
}




=head2 set_env_variable

Sets an environment variable to a given value.

=cut

sub set_env_variable {
    my ( $variable, $value ) = @_;
    $ENV{$variable} = $value;
    return $ENV{$variable};
}

=head2 get_env_variable

Gets the value of an environment variable.

=cut

sub get_env_variable {
    my ($variable) = @_;
    return $ENV{$variable};
}

=head2 check_env_variable

Checks if an environment variable is set and optionally prints its value.

=cut

sub check_env_variable {
    my ( $variable, $print ) = @_;
    my $value = $ENV{$variable};
    if ( defined $value ) {
        if ( defined $print ) {
            if ($print) {
                print "\$$variable: \t$value\n";
            }
        }
        return $value;
    }
    else {
        print "\$$variable is not set\n";
        exit();
    }
}

=head2 set_and_check_env_variable

Sets an environment variable and then checks if it is set.

=cut

sub set_and_check_env_variable {
    my ( $variable, $value, $print ) = @_;
    set_env_variable( $variable, $value );
    return check_env_variable( $variable, $print );

    # return $value;
}

=head2 calculate_directories

Calculates various directories based on the base directory and batch name.

=cut

sub calculate_directories {
    my ($config_ref) = @_;

    my $data_dir =
      $config_ref->{'data_dir'};    # Assuming you have this in your config
    my $batch_name = $config_ref->{'batch_name'};

    # $basedir    =~ s/^\s+|\s+$//g;
    $batch_name =~ s/^\s+|\s+$//g;
    if ( $config_ref->{'adapt'} && index( $batch_name, "adapt" ) == -1 ) {
        $batch_name = $batch_name . "_adapt";
    }

    use File::Spec::Functions;





    my $fluxdir = $config_ref->{'fl_mhdlib'};
    my $pipedir = catdir( $fluxdir, "fluxpipe", "fluxpipe" );
    my $pdldir  = catdir( $fluxdir, "pdl",      "PDL" );

    # # Use the provided data_dir if defined, otherwise calculate it
    # my $datdir =
    #   defined($data_dir) ? $data_dir : catdir( $basedir, "fluxon-data" );

    my $magdir   = catdir( $data_dir, "magnetograms" );
    my $batchdir = catdir( $data_dir, "batches", $batch_name );
    my $logfile  = catfile( $batchdir, "pipe_log.txt" );


    #Remove the tilde from the path
    use File::HomeDir;
    my $home_dir = $ENV{'HOME'};
    $data_dir =~ s{^~}{$home_dir};
    $fluxdir =~ s{^~}{$home_dir};
    $pdldir =~ s{^~}{$home_dir};
    $magdir =~ s{^~}{$home_dir};
    $batchdir =~ s{^~}{$home_dir};
    $logfile =~ s{^~}{$home_dir};

    # use File::Glob ':glob';
    # # Replace ~ with the home directory
    # my $home_dir = bsd_glob("~");
    # $data_dir = bsd_glob($data_dir);
    # $fluxdir = bsd_glob($fluxdir);
    # $pdldir = bsd_glob($pdldir);
    # $magdir = bsd_glob($magdir);
    # $batchdir = bsd_glob($batchdir);
    # $logfile = bsd_glob($logfile);


    # Update the original config hash
    $config_ref->{'pipe_dir'}  = $pipedir;
    $config_ref->{'pdl_dir'}   = $pdldir;
    $config_ref->{'datdir'}    = $data_dir;
    $config_ref->{'data_dir'}  = $data_dir;
    $config_ref->{'mag_dir'}   = $magdir;
    $config_ref->{'batch_dir'} = $batchdir;
    $config_ref->{'logfile'}   = $logfile;

    set_and_check_env_variable( 'DATAPATH', $data_dir, 0 );
}

=head2 set_python_path

Sets the PYTHONPATH environment variable.

=cut

sub set_python_path {
    my ( $pythonpath, $print ) = @_;
    set_and_check_env_variable( 'PYTHONPATH', $pythonpath, $print );
    return $pythonpath;
}

=head2 print_banner

Prints a banner with various details.

=cut

sub print_banner {
    my ( $batch_name, $CR, $reduction, $n_fluxons_wanted, $recompute_string,
        $adapt, $flow_method)
      = @_;
    print "\n\n\n\n\n\n|\n|\n|\n|\n|\n|\n|\n|\n|\n|";
    print "\n\n";
    print
"--------------------------------------------------------------------------------------------------\n";
    print
"FLUXPipe: Indicate a Carrington Rotation and this script will run the entire Flux Pipeline for it.\n";
    print
"--------------------------------------------------------------------------------------------------\n";
    print "\n\n";

    check_env_variable( 'DATAPATH', 1 );
    print
"\nBatch: $batch_name, CR: $CR, Reduction: $reduction, Fluxons: $n_fluxons_wanted, Adapt: $adapt, Wind: $flow_method\n";

    my $time  = localtime;
    my $ftime = $time->strftime('%m-%d-%Y %H:%M:%S');

    print "\n\n";
    print
"\t>>>>>>>>>>>>>>>>>>>>> Recompute = $recompute_string <<<<<<<<<<<<<<<<<<<<<<";
    print "\n\tStarting FLUXPipe at $ftime ";
    print
      "\n\t>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
    print "\n\n";
    print
"--------------------------------------------------------------------------------------------------\n";
    return 1;
}

=head2 search_files_in_directory

Searches for files in a directory that match a known string and file extension.

=cut

sub search_files_in_directory {
    my ( $directory, $known_string, $extension ) = @_;

    # Escape the known string to avoid regex special characters
    my $escaped_string = quotemeta $known_string;

    # Generate the regular expression pattern
    my $pattern = $escaped_string . '.*' . $extension;

    opendir( my $dh, $directory ) or die "Failed to open directory: $!";
    while ( my $file = readdir($dh) ) {
        next if ( $file =~ /^\./ );    # Skip hidden files/directories
        next unless ( $file =~ /$pattern/ );
        print "$file\n";               # Process the matching file
        closedir($dh);

        # print $file;
        return $file;
    }
}

=head2 check_second_file_presence

Checks for the presence of a second file related to the given file path.

=cut

sub check_second_file_presence {
    my ($file_path) = @_;

    # Get the directory and base name of the first file
    my $directory = dirname($file_path);
    my $file_name = basename($file_path);

    # Generate the pattern for the second file
    my $second_file_pattern = $file_name;
    print "    File name: $file_name\n";
    $second_file_pattern =~ s/(\.[^.]+)$/_relaxed_.*${1}/;

    opendir( my $dh, $directory ) or die "Failed to open directory: $!";
    while ( my $file = readdir($dh) ) {
        next if ( $file =~ /^\./ );       # Skip hidden files/directories
        next if ( $file =~ /\.png$/ );    #skip png files

        if ( $file =~ /^$second_file_pattern$/ ) {
            closedir($dh);
            return 1, $file;              # Second file found
        }

        # print("$file is wrong\n");
    }
    closedir($dh);

    return 0, 0;    # Second file not found
}

# Check if the environment variable is set
sub set_paths {
    my ($do_plot) = @_;
    if ( defined $ENV{'FL_PREFIX'} ) {
        my $envpath =
          "$ENV{'FL_PREFIX'}/fluxpipe/fluxpipe/helpers/perl_paths.pm";

        # Check if the file exists and is readable
        if ( -e $envpath && -r _ ) {
            require $envpath;
            print "\n\n";
        }
        else {
            warn "File does not exist or is not readable: $envpath\n\n";
        }
    }
    else {
        warn "Environment variable FL_PREFIX is not set.\n\n";
    }

    # print the lists of directories
    if ($do_plot) {
        my @PDL_INC;
        my @PDLLIB;
        print "\n\nINC has:\n ";
        print map { " $_\n" } @INC;
        print "--------------------------\n";

        print "\nPDL_INC has:\n ";
        print map { " $_\n" } @PDL_INC;
        print "--------------------------\n";

        print "\nPDLLIB has:\n ";
        print map { " $_\n" } @PDLLIB;
        print "--------------------------\n\n";

        # Print each command-line argument
        foreach my $arg (@ARGV) {
            print "Argument: $arg\n";
        }
        print "\n\n";
    }
    return;
}

use File::Spec;
use List::Util qw(max);


sub load_highest_numbered_world {
    my ($datdir, $batch_name, $CR, $n_fluxons_wanted, $inst) = @_;

    my $world_out_dir = File::Spec->catdir($datdir, "batches", $batch_name, "data", "cr${CR}", "world");
    my $file_pattern = qr/cr${CR}_f${n_fluxons_wanted}_${inst}_relaxed_s(\d+)\.flux$/;
    my $original_pattern = qr/cr${CR}_f${n_fluxons_wanted}_${inst}\.flux$/;

    my $max_d = -1;
    my $selected_file_path;
    my $original_file_path;
    print "Searching for files in $world_out_dir\n";
    print "File pattern: $file_pattern\n";
    print "Original pattern: $original_pattern\n";


    opendir(my $dh, $world_out_dir) or die "Cannot open directory: $!";
    while (my $file = readdir($dh)) {
        print $file . "\n";
        if ($file =~ /$file_pattern/) {
            my $d_value = $1;
            if ($d_value > $max_d) {
                $max_d = $d_value;
                $selected_file_path = File::Spec->catfile($world_out_dir, $file);
            }
        }
        if ($file =~ /$original_pattern/) {
            $original_file_path = File::Spec->catfile($world_out_dir, $file);
        }
    }
    closedir($dh);

    print "Selected file: $selected_file_path\n";
    print "Original file: $original_file_path\n";

    if (defined $selected_file_path & defined $original_file_path) {
        my $this_world_relaxed = read_world($selected_file_path);
        my $this_world_original = read_world($original_file_path);
        my @fluxons = $this_world_relaxed->fluxons;

        if (scalar @fluxons == 0) {
            # Consider logging a warning or handling this case differently as needed
            return die "World loaded, but contains no fluxons.";
        }

        return $this_world_relaxed, $this_world_original; # Successful load
    } else {
            die "No matching files found."; # No file found
    }
}

1;
