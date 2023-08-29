=head1 NAME

YourScriptName - Utility Functions for File and Environment Management

=head1 SYNOPSIS

    use pdl_pipe_helper;
    shorten_path($string);
    find_highest_numbered_file($directory);
    set_env_variable($variable, $value);
    # ... and so on

=head1 DESCRIPTION

This Perl script provides utility functions for managing files and environment variables.
It includes functions for shortening file paths, finding the highest-numbered file in a directory,
setting and checking environment variables, and more.

=head1 FUNCTIONS

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

    calculate_directories($basedir, $batch_name, $print);

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

Your Name <your.email@example.com>

=head1 SEE ALSO

L<PDL::AutoLoader>, L<PDL>, L<Time::Piece>

=cut



# use strict;
use warnings;
use PDL::AutoLoader;
use PDL;
use Time::Piece;

# our @PDLLIB;
# my @INC;

sub shorten_path {
    my ($string) = @_;
    my $datapath = $ENV{'DATAPATH'};
    if ($datapath) {
        $string =~ s/\Q$datapath\E/\$DATAPATH/g;
    }
    return $string;
}


sub find_highest_numbered_file {
    my ($directory) = @_;
    opendir(my $dir_handle, $directory) or die "Cannot open directory: $!";
    my @files = grep { !/^\.{1,2}$/ } readdir($dir_handle);
    closedir $dir_handle;

    my $highest_numbered_file;
    my $highest_number = -1;
    my $found = 0;
    for my $file_name (@files) {
        # Match file names containing "_relaxed", a number, and ".flux"
        if ($file_name =~ /_relaxed_s(\d+)\.flux/) {
            my $number = $1;
            if ($number > $highest_number) {
                $highest_number = $number;
                $highest_numbered_file = $file_name;
            }
            $found = 1;
        }
    }
    return $highest_numbered_file ? "$directory$highest_numbered_file" : 0, $highest_number;
}


sub set_env_variable {
    my ($variable, $value) = @_;
    $ENV{$variable} = $value;
    return $ENV{$variable};
}


sub get_env_variable {
    my ($variable) = @_;
    return $ENV{$variable};
}


sub check_env_variable {
    my ($variable, $print) = @_;
    my $value = $ENV{$variable};
    if (defined $value) {
        if (defined $print) {
            if ($print) {
            print "\$$variable: \t$value\n";}}
        return $value;
    } else {
        print "\$$variable is not set\n";
        exit();
    }
}


sub set_and_check_env_variable {
    my ($variable, $value, $print) = @_;
    set_env_variable($variable, $value);
    return check_env_variable($variable, $print);
    # return $value;
}


sub calculate_directories {
(my $basedir, my $batch_name, my $print) = @_;

    # Calculate Directory Structure
    my $fluxdir =  "$basedir/fluxon-mhd";
        my $pipedir =  "$fluxdir/flux-pipe";
        my $pdldir =   "$fluxdir/pdl/PDL";

    my $datdir =   "$basedir/fluxon-data";
        my $magdir =   "$datdir/magnetograms";
        my $batchdir = "$datdir/batches/$batch_name";
            my $logfile = "$batchdir/pipe_log.txt";

    set_and_check_env_variable('FLUXPATH', $fluxdir, $print);
    set_and_check_env_variable('PIPEPATH', $pipedir, $print);
    set_and_check_env_variable('BATCHPATH', $batchdir, $print);
    set_and_check_env_variable('DATAPATH', $datdir, $print);

    return ($pipedir, $pdldir, $datdir, $magdir, $batchdir, $logfile);
}


sub set_python_path {
    my ($pythonpath, $print) = @_;
    set_and_check_env_variable('PYTHONPATH', $pythonpath, $print);
    return $pythonpath;
}


sub print_banner {
    my ($batch_name, $CR, $reduction, $n_fluxons_wanted, $recompute_string) = @_;
    print "|\n|\n|\n|\n|\n|\n|\n|\n|\n|";
    print "\n\n";
    print "--------------------------------------------------------------------------------------------------\n";
    print "FLUXPipe: Indicate a Carrington Rotation and this script will run the entire Flux Pipeline for it.\n";
    print "--------------------------------------------------------------------------------------------------\n";
    print "\n\n";

    check_env_variable('DATAPATH', 1);
    print "\nBatch: $batch_name, CR: $CR, Reduction: $reduction, Fluxons: $n_fluxons_wanted";

    my $time = localtime;
    my $ftime = $time->strftime('%m-%d-%Y %H:%M:%S');

    print"\n\n";
    print "\t>>>>>>>>>>>>>>>>>>>>> Recompute = $recompute_string <<<<<<<<<<<<<<<<<<<<<<";
    print "\n\tStarting FLUXPipe at $ftime ";
    print "\n\t>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
    print "\n\n";
    print "--------------------------------------------------------------------------------------------------\n";
    return 1;
}


sub search_files_in_directory {
    my ($directory, $known_string, $extension) = @_;

    # Escape the known string to avoid regex special characters
    my $escaped_string = quotemeta $known_string;

    # Generate the regular expression pattern
    my $pattern = $escaped_string . '.*' . $extension;

    opendir(my $dh, $directory) or die "Failed to open directory: $!";
    while (my $file = readdir($dh)) {
        next if ($file =~ /^\./);  # Skip hidden files/directories
        next unless ($file =~ /$pattern/);
        print "$file\n";  # Process the matching file
    }
    closedir($dh);
}


sub check_second_file_presence {
    my ($file_path) = @_;

    # Get the directory and base name of the first file
    my $directory = dirname($file_path);
    my $file_name = basename($file_path);

    # Generate the pattern for the second file
    my $second_file_pattern = $file_name;
    print "    File name: $file_name\n";
    $second_file_pattern =~ s/(\.[^.]+)$/_relaxed_.*${1}/;

    opendir(my $dh, $directory) or die "Failed to open directory: $!";
    while (my $file = readdir($dh)) {
        next if ($file =~ /^\./);  # Skip hidden files/directories
        next if ($file =~ /\.png$/); #skip png files

        if ($file =~ /^$second_file_pattern$/) {
            closedir($dh);
            return 1, $file;  # Second file found
        }
        # print("$file is wrong\n");
    }
    closedir($dh);

    return 0, 0;  # Second file not found
}

1;