
=head1 NAME

get_wind - Calculate Solar Wind Plasma Parameters for a Given Carrington Rotation (CR)

=head1 SYNOPSIS

    use pipe_helper;
    my ($out_b, $out_fr, $out_wind) = get_wind($this_world_relaxed, $datdir, $batch_name, $CR, $N_actual, $recompute, $n_want, $pythondir);

=cut

package get_wind;
use strict;
use warnings;
use Exporter qw(import);
our @EXPORT_OK = qw(get_wind);
use pipe_helper                     qw(shorten_path);
use File::Path                      qw(mkpath);
use map_fluxon_b                    qw(map_fluxon_b);
use map_fluxon_fr                   qw(map_fluxon_fr);
use map_fluxon_flow_parallel_master qw(map_fluxon_flow_parallel_master);

=head1 DESCRIPTION

This subroutine calculates various solar wind plasma parameters such as the radial magnetic field, radial expansion factor, and radial wind speed. It performs these calculations based on the provided Carrington Rotation (CR) and other parameters.

=head1 PARAMETERS

=over 4

=item * C<$this_world_relaxed> - The relaxed world object containing fluxons.

=item * C<$datdir> - Data directory path.

=item * C<$batch_name> - Name of the batch.

=item * C<$CR> - Carrington Rotation number.

=item * C<$N_actual> - Actual number of footpoints.

=item * C<$recompute> - Flag to indicate whether to recompute the parameters or not.

=item * C<$n_want> - Number of footpoints wanted.

=item * C<$pythondir> - Python directory path (currently not used).

=back

=head1 WORKFLOW

1. Checks for the existence of output directories and files.
2. If necessary, creates output directories.
3. If the output files do not exist or if recompute is true, performs the following calculations:
    - Updates neighbors in the relaxed world.
    - Calculates the radial magnetic field.
    - Calculates the radial expansion factor.
    - Calculates the radial wind speed.

=head1 OUTPUT

Returns the paths to the output files for radial magnetic field (C<$out_b>), radial expansion factor (C<$out_fr>), and radial wind speed (C<$out_wind>).

=head1 EXCEPTIONS

Dies if it fails to create the output directory.

=head1 AUTHOR

Gilly <gilly@swri.org> (and others!)

=head1 SEE ALSO

L<pipe_helper>, L<map_fluxon_b>, L<map_fluxon_fr>, L<map_fluxon_flow_parallel_master>

=cut

sub get_wind {

    my ( $this_world_relaxed, $datdir, $batch_name, $CR, $N_actual, $recompute,
        $n_want, $pythondir )
      = @_;

    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    print "(pdl) Calculating Solar Wind Plasma Parameters for CR$CR\n";
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

    my $do_wind_calc = 0;

    my $wind_out_dir   = $datdir . "/batches/$batch_name/cr" . $CR . '/wind';
    my $prefix         = "$wind_out_dir/cr$CR\_f$n_want";
    my $out_b          = "$prefix\_radial_bmag.dat";
    my $out_fr         = "$prefix\_radial_fr.dat";
    my $out_wind       = "$prefix\_radial_wind.dat";
    my $short_out_wind = shorten_path($out_wind);
    my $skipstring =
      "\n\tWind Calculation Skipped! \n\t\tFound $short_out_wind\n\n";

    # Make the directory if necessary
    if ( !-d $wind_out_dir ) {
        $do_wind_calc = 1;
        mkpath($wind_out_dir)
          or die "Failed to create directory: $wind_out_dir $!\n";
    }

    # Check if the files exist
    if ( !-f $out_b or !-f $out_fr or !-f $out_wind or $recompute ) {
        $do_wind_calc = 1;
    }

    # Perform the calculation if necessary
    if ( $do_wind_calc || $recompute ) {

        #Initialize the world
        print "\n\tUpdating neighbors...";
        $this_world_relaxed->update_force(0);
        my @fluxons = $this_world_relaxed->fluxons;
        print "Done!\n";

        # Calculate the radial magnetic field
        print "\n\tRadial Magnetic Field (B) Calculation...";
        map_fluxon_b( $out_b, \@fluxons );
        print "Done!\n";

        # print "\t\t...done with radial B!";

        # Calculate the radial expansion factor
        print "\n\tRadial Expansion Factor (Fr) Calculation...";
        map_fluxon_fr( $out_fr, \@fluxons );
        print "Done!";

        # Calculate the radial wind speed
        no warnings 'misc';
        print "\n\n\tRadial Wind Speed Calculation...\n";
        my $do_wind_map = 0 || $recompute;

        # $do_wind_map=1; #OVERRIDE WIND MAP

        if ( !-e $out_wind ) { $do_wind_map = 1; }
        if ($do_wind_map) {
            {
                no warnings;    # Suppress warnings
                eval {
                    map_fluxon_flow_parallel_master( $out_wind, \@fluxons );
                };
            };

        }
        else {
            print $skipstring;
        }
        use warnings 'misc';
    }
    else {
        print $skipstring;
    }

    print "\t\t\t```````````````````````````````\n\n\n\n";

    return $out_b, $out_fr, $out_wind;
}
1;