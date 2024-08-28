
=head1 NAME

make_world_from_pfss - Take the PFSS field lines and make a FLUX world from them.
=cut

package make_world_from_pfss;
use strict;
use warnings;
use Exporter qw(import);
our @EXPORT_OK = qw(make_world_from_pfss);
use make_world_sphere qw(make_world_sphere);
use PDL::AutoLoader;
use PDL;
use PDL::Transform;
use PDL::NiceSlice;
use PDL::Options;
use Flux;
use PDL::IO::Misc;
use File::Path;
use Time::HiRes qw(clock_gettime);
use File::Basename qw(fileparse);
use File::Basename;
use pipe_helper qw(shorten_path);
use plot_world qw(plot_world);
use Flux::World    qw(read_world);
use Flux::World qw(str2world);
use PDL::Graphics::Gnuplot qw(gpwin);
use File::Path  qw(mkpath);


# use Flux::World "str2world";
# str2world();

=head1 SYNOPSIS

    use make_world_from_pfss;
    make_world_from_pfss($world_out_dir, $floc_path, $open_file, $closed_file,
    $force_make_world, $CR, $datdir, $batch, $N_actual, $n_fluxons_wanted, $lim, $lim2, $n_outliers);

=head1 DESCRIPTION

This script uses PDL to manipulate and visualize magnetic field lines in the "world" representation from FLUX.
It reads field lines from files, transforms them into a specific format, and then saves and visualizes them.

=head1 FUNCTIONS

=head2 make_world_from_pfss

    make_world_from_pfss($world_out_dir, $floc_path, $open_file, $closed_file,
    $force_make_world, $CR, $datdir, $batch, $N_actual, $n_fluxons_wanted, $lim, $lim2, $n_outliers);

This function does the following:

=over

=item * Reads the PFSS field lines from the given files.

=item * Transforms the field lines into the FLUX world format.

=item * Saves the transformed field lines to a file.

=item * Optionally visualizes the field lines.

=back

=head3 PARAMETERS

=over

=item * C<$world_out_dir>: Directory where the output world file will be saved.

=item * C<$floc_path>: File containing the field lines locations.

=item * C<$open_file>: File containing open field lines.

=item * C<$closed_file>: File containing closed field lines.

=item * C<$force_make_world>: Flag to force the creation of a new world.

=item * C<$CR>: Carrington Rotation

=item * C<$datdir>: Data Directory

=item * C<$batch>: Batch Name

=item * C<$N_actual>: Actual number of fluxons.

=item * C<$n_fluxons_wanted>: Number of fluxons wanted.

=item * C<$lim>: Limit for the plot range, usually zoomed in.

=item * C<$lim2>: Another limit for the plot range, usually wider angle.

=item * C<$n_outliers>: Number of outliers.

=back

=head3 OUTPUT

This function does not return any value.

=head1 AUTHOR

Gilly <gilly@swri.org> (and others!)

=head1 SEE ALSO

L<PDL>, L<PDL::Transform>, L<PDL::NiceSlice>, L<PDL::Options>, L<Flux>, L<PDL::IO::Misc>, L<File::Path>, L<Time::HiRes>, L<File::Basename>

=cut


sub make_world_from_pfss {
    my ($datdir, $batch_name, $CR, $reduction, $n_fluxons_wanted, $adapt, $force_make_world, $lim, $lim2, $configs) = @_;

    # Define the output directory and floc path
    my $world_out_dir = "$datdir/batches/$batch_name/data/cr$CR/world";
    my $floc_path = "$datdir/batches/$batch_name/data/cr$CR/floc";
    my $file_end;

    if ($adapt) {
        $reduction = "f".$configs->{'adapt_select'};
        $file_end = "adapt"
    } else {
        # $reduction = "R";
        $file_end = "hmi"
    }

    my $open_file = "$floc_path/floc_open_cr$CR\_r$reduction\_f$n_fluxons_wanted\_$file_end.dat";
    my $closed_file = "$floc_path/floc_closed_cr$CR\_r$reduction\_f$n_fluxons_wanted\_$file_end.dat";
    my $world_out_path = $world_out_dir .'/cr'.$CR.'_f'.$n_fluxons_wanted."\_$file_end.flux";

    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    print "(pdl) Converting PFSS Fieldlines into FLUX Fluxons\n";
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

    # Make the world
    my $need_world=0;
    if (! -f $world_out_path ) {
        $need_world=1;
        }

    if ($force_make_world || $need_world){
        ## Read python-processed data ###############################
        my ($oflnum, $oflx, $olat, $olon, $orad) = rcols $open_file;
        my ($cflnum, $cflx, $clat, $clon, $crad) = rcols $closed_file;

        my $open_name_short = shorten_path($open_file, 5);
        my $closed_name_short = shorten_path($closed_file, 5);


        ## Define a blank world ###############################
        my @flines = ();
        my $world = ();
        my $xform = !t_spherical() x t_scale([3.14159/180,3.14159/180,1]);


        ## Generate open fieldlines ####################################
        print "\n \n\tGenerating open fluxons from \n\t\t$open_name_short\n\t ";
        my ($ofln, $ofl) = rle($oflnum);
        $ofln = $ofln->cumusumover;
        $ofln = append(0, $ofln);
        my $open_count = 1;

        # my $extend_down = 0;
        # my $extend_up = 0;
        # my $num_replications = 7;
        # my $power = -3;

        for my $i(0..nelem($ofl)-1){
            # Extract latitude, longitude, and radial arrays for the current flux tube
            my $flxlat = $olat($ofln($i):$ofln($i+1)-1);
            my $flxlon = $olon($ofln($i):$ofln($i+1)-1);
            my $flxrad = $orad($ofln($i):$ofln($i+1)-1);
            my $open = pdl($flxlon, $flxlat, $flxrad)->transpose;


            # my $print_debug = 0;
            # if ($extend_down) {
            #     if ($print_debug) {print "EXTEND DOWN\n";}
            #     # Deep copy of the first column
            #     my @replicated_columns;
            #     for (1..$num_replications) {
            #         push @replicated_columns, $open->(:, 0)->copy;
            #     }
            #     my $replicated_first_col = cat(@replicated_columns)->squeeze;
            #     if ($print_debug) {print "Replicated first column: $replicated_first_col\n";}

            #     # Calculate the second lowest value in the radial component
            #     my $sorted_radial = qsort($open(2, :));
            #     if ($print_debug) {print "Sorted radial values: $sorted_radial\n";}

            #     # Ensure there are at least two unique values
            #     my $unique_values = $sorted_radial->uniq;
            #     if (nelem($unique_values) < 2) {
            #         die "Not enough unique values in the radial component to determine the second lowest value.\n";
            #     }

            #     my $second_min_radial = $unique_values(1);  # The second lowest unique value
            #     if ($print_debug) {print "Second lowest radial value: $second_min_radial\n";}

            #     # Calculate the logarithmically spaced values
            #     my $log_min = log10(1 + 10**$power);
            #     my $log_max = log10($second_min_radial);
            #     if ($print_debug) {print "Log min: $log_min, Log max: $log_max\n";}

            #     my $log_values = ($log_min + (pdl(sequence($num_replications)) / ($num_replications - 1)) * ($log_max - $log_min));
            #     my $log_spaced_values = 10**$log_values;
            #     if ($print_debug) {print "Log spaced values: $log_spaced_values\n";}

            #     # Ensure the dimensions match for assignment
            #     $log_spaced_values = $log_spaced_values->reshape($num_replications, 1);
            #     if ($print_debug) {print "Shape of replicated_first_col before assignment: " . $replicated_first_col->shape . "\n";}
            #     if ($print_debug) {print "Shape of log_spaced_values: " . $log_spaced_values->shape . "\n";}

            #     # Directly modify the radial component of the replicated columns
            #     for my $j (0..($num_replications-1)) {
            #         $replicated_first_col(2, $j) .= $log_spaced_values(($j));
            #     }
            #     if ($print_debug) {print "Modified replicated_first_col: $replicated_first_col\n";}

            #     # Remove the last row from replicated_first_col
            #     $replicated_first_col = $replicated_first_col(:, 0:-2);
            #     if ($print_debug) {print "Replicated first column after removing last row: $replicated_first_col\n";}

            #     # Glue the replicated columns to the beginning of the original PDL object
            #     $open = $replicated_first_col->glue(1, $open)->squeeze;
            #     if ($print_debug) {print "Extended open: $open\n";}

            # }

            # my $power_up = 2;
            # $print_debug = 0;
            # if ($extend_up) {
            #     if ($print_debug) {print "EXTEND UP\n";}
            #     # Deep copy of the last column
            #     my @replicated_columns;
            #     for (1..$num_replications) {
            #         push @replicated_columns, $open->(:, -1)->copy;
            #     }
            #     my $replicated_last_col = cat(@replicated_columns)->squeeze;
            #     if ($print_debug) {print "Replicated last column: $replicated_last_col\n";}

            #     # Calculate the highest value in the radial component
            #     my $max_radial = $open(2, :)->max;
            #     if ($print_debug) {print "Maximum radial value: $max_radial\n";}

            #     # Calculate the logarithmically spaced values
            #     my $log_min = log10($max_radial);
            #     my $log_max = log10(1 + 10**$power_up);
            #     if ($print_debug) {print "Log min: $log_min, Log max: $log_max\n";}

            #     my $log_values = ($log_min + (pdl(sequence($num_replications)) / ($num_replications - 1)) * ($log_max - $log_min));
            #     my $log_spaced_values = 10**$log_values;
            #     if ($print_debug) {print "Log spaced values: $log_spaced_values\n";}

            #     # Ensure the dimensions match for assignment
            #     $log_spaced_values = $log_spaced_values->reshape($num_replications, 1);
            #     if ($print_debug) {print "Shape of replicated_last_col before assignment: " . $replicated_last_col->shape . "\n";}
            #     if ($print_debug) {print "Shape of log_spaced_values: " . $log_spaced_values->shape . "\n";}

            #     # Directly modify the radial component of the replicated columns
            #     for my $j (0..($num_replications-1)) {
            #         $replicated_last_col(2, $j) .= $log_spaced_values(($j));
            #     }
            #     if ($print_debug) {print "Modified replicated_last_col: $replicated_last_col\n";}

            #     # Remove the first row from replicated_last_col
            #     $replicated_last_col = $replicated_last_col(:, 1:);
            #     if ($print_debug) {print "Replicated last column after removing first row: $replicated_last_col\n";}

            #     # Glue the replicated columns to the end of the original PDL object
            #     $open = $open->glue(1, $replicated_last_col)->squeeze;
            #     if ($print_debug) {print "Extended open: $open\n";}
            # }


            # # Print radial values before sorting
            # if ($print_debug) {print "Radial values before sorting: ", $open(2, :), "\n";}

            # # Sort the array by the radial component to ensure monotonic increase
            # my $sorted_indices = $open(2, :)->flat->qsorti;
            # if ($print_debug) {print "Sorted indices: $sorted_indices\n";}

            # # Reorder the array based on the sorted indices
            # $open = $open->dice_axis(1, $sorted_indices);

            # # Remove the first element, which is a 1.0 value
            # $open = $open->(:, 1:-1);


            if ($oflx->($ofln($i))<0){
                $open = $open->(:,-1:0:-1);
            }

            my $line = ($open)->apply($xform);

            if (0) {print "open: $open\n";}
            if (0) {print "line: $line\n";
                    <STDIN>;}
            push @flines,$line->copy;
            $open_count++;
        }
        print "\tDone! $open_count open fluxons generated.\n";


        ## Generate closed fieldlines  ###################################
        print "\n \n\tGenerating closed fluxons from \n\t\t$closed_name_short\n\t ";
        my ($cfln, $cfl) = rle($cflnum);
        $cfln = $cfln->cumusumover;
        $cfln = append(0, $cfln);
        my $closed_count = 1;
        for my $j(0..nelem($cfl)-1){
            my $flxlat = $clat($cfln($j):$cfln($j+1)-1);
            my $flxlon = $clon($cfln($j):$cfln($j+1)-1);
            my $flxrad = $crad($cfln($j):$cfln($j+1)-1);
            my $closed = pdl($flxlon, $flxlat, $flxrad)->transpose;
            if ($cflx->($cfln($j))<0){
                $closed = $closed->(:,-1:0:-1);
            }
            my $line = ($closed)->apply($xform);
            push(@flines,$line->copy);
            $closed_count++;
        }
        my $closed_nflux = $closed_count * 2 ;
        print("\tDone! $closed_count closed fluxons generated (from $closed_nflux footpoints)\n");

        my $total_fluxons = $open_count + $closed_count;
        print "\n\t    Total fluxons created: $total_fluxons\n";






    ## Generate the world  ############################################
        print "\n\tGenerating the World...\n \n";
        my $fbg = make_world_sphere(@flines, {rmax=>21.5});
        # my $fbg = make_world_sphere(@flines, {rmax=>215});

        $world = str2world($fbg);


        ## Save the initial world state  ###############################
        print "\n\tSaving the World...";
        my $flen = $oflnum->max() + $cflnum->max() + 4; #@flines+0;
        if (! -d $world_out_dir ) {mkpath($world_out_dir) or die "Failed to create directory: $world_out_dir $!\n";}
        $world->write_world($world_out_path);
        my $short_world_out_path = shorten_path($world_out_path);
        print "\n\t    Saved to $short_world_out_path\n";


        ## Plot the initial World State ################################
        plot_world($world, $datdir, $batch_name, $CR, $reduction, $n_fluxons_wanted, $adapt, $force_make_world, $lim, $lim2, $configs, "initial");


    } else {
        print "\n\n \tSkipped! World already exists:\n";
        my $short_world_out_path = shorten_path($world_out_path);
        print "\t    $short_world_out_path\n";

        print "\t Plotting World!\n";
        my $world = read_world($world_out_path);
        plot_world($world, $datdir, $batch_name, $CR, $reduction, $n_fluxons_wanted, $adapt, $force_make_world, $lim, $lim2, $configs, "initial");

    }



        print "\n\t\t\t```````````````````````````````\n \n\n\n";
    return $open_file, $closed_file, $world_out_path;
}


1;