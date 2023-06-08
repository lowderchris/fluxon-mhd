

use strict;
use warnings;
use PDL::AutoLoader;
use PDL;
use Time::Piece;


sub shorten_path_old{
    my ($start_path, $levels) = @_;
    $levels = $levels || 5;
    my $out_string = (split('/', $start_path))[-$levels .. -1] ? join('/', (split('/', $start_path))[-$levels .. -1]) : $start_path;
    return "DATAPATH/" . $out_string;
}

sub shorten_path {
    my ($string) = @_;
    my $datapath = $ENV{'DATAPATH'};
    if ($datapath) {
        $string =~ s/\Q$datapath\E/\$DATAPATH/g;
    }
    return $string;
}

sub find_file_with_string {
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
    return $highest_numbered_file ? "$directory$highest_numbered_file" : 0;
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
    
    set_and_check_env_variable('FLUXPATH', $fluxdir, 0);
    set_and_check_env_variable('PIPEPATH', $pipedir, $print);
    set_and_check_env_variable('BATCHPATH', $batchdir, $print);
    # print "\n";
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

    # check_env_variable('PYTHONPATH', 1);
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

}


# system("cd $pipedir");

1;