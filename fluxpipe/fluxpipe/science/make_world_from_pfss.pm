
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
    my $world_out_dir = "$datdir/batches/$batch_name/cr$CR/world";
    my $floc_path = "$datdir/batches/$batch_name/cr$CR/floc";
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
        for my $i(0..nelem($ofl)-1){
            # print "$i" . "\n";
            my $flxlat = $olat($ofln($i):$ofln($i+1)-1);
            my $flxlon = $olon($ofln($i):$ofln($i+1)-1);
            my $flxrad = $orad($ofln($i):$ofln($i+1)-1);
            my $open = pdl($flxlon, $flxlat, $flxrad)->transpose;

            if ($oflx->($ofln($i))<0){
                $open = $open->(:,-1:0:-1);
            }
            my $line = ($open)->apply($xform);
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
        $world = str2world($fbg);


        ## Save the initial world state  ###############################
        print "\n\tSaving the World...";
        my $flen = $oflnum->max() + $cflnum->max() + 4; #@flines+0;
        if (! -d $world_out_dir ) {mkpath($world_out_dir) or die "Failed to create directory: $world_out_dir $!\n";}
        $world->write_world($world_out_path);
        my $short_world_out_path = shorten_path($world_out_path);
        print "\n\t    Saved to $short_world_out_path\n";


        ## Display ####################################################
        print "\n \n\tPlotting the World...";
        # Set Ranges of Plots
        my $range_i = [-$lim, $lim, -$lim, $lim, -$lim, $lim];
        my $range_f = [-$lim, $lim, -$lim, $lim, -$lim, $lim];
        my $range_f2 = [-$lim2, $lim2, -$lim2, $lim2, -$lim2, $lim2];

        # number of open fluxons, number of fluxons, number of fluxons requested
        my $top_dir = $datdir."/batches/$batch_name/imgs/initial/";
        # my $wide_dir = $datdir."/batches/$batch/imgs/world/wide/";
        # my $narrow_dir = $datdir."/batches/$batch/imgs/world/narrow/";



        if (! -d $top_dir ) {mkpath($top_dir) or die "Failed to create directory: $top_dir $!\n";}
        my $wide_dir = $world_out_dir ."wide/";
        if (! -d $wide_dir ) {mkpath($wide_dir) or die "Failed to create directory: $wide_dir $!\n";}
        my $narrow_dir = $world_out_dir ."narrow/";
        if (! -d $narrow_dir ) {mkpath($narrow_dir) or die "Failed to create directory: $narrow_dir $!\n";}

        my $ext = 'png';
        my $renderer = $ext.'cairo';
        # my $filename
        my $world_png_path = $narrow_dir."cr$CR\_f". $n_fluxons_wanted. "_initial_pfss.$ext";
        my $world_png_path2= $wide_dir."cr$CR\_f". $n_fluxons_wanted. "_initial_pfss_wide.$ext";
        my $world_png_path_top = $top_dir   ."cr$CR\_f". $n_fluxons_wanted. "_initial_pfss.$ext";

        # my $window00 = gpwin($renderer,size=>[9,9], dashed=>0, output=> $world_png_path);
        # $world->render( {'window'=>$window00, range=>$range_i});
        my $window000 = gpwin($renderer,size=>[9,9], dashed=>0, output=> $world_png_path_top);
        $world->render( {'window'=>$window000, range=>$range_i});
        # my $window01 = gpwin($renderer,size=>[9,9], dashed=>0, output=> $world_png_path2);
        # $world->render( {'window'=>$window01, range=>$range_f2});
        print "Done!\n";

    } else {
        print "\n\n \tSkipped! World already exists:\n";
        my $short_world_out_path = shorten_path($world_out_path);
        print "\t    $short_world_out_path\n";
    }
        print "\n\t\t\t```````````````````````````````\n \n\n\n";
    return $open_file, $closed_file, $world_out_path;
}
1;