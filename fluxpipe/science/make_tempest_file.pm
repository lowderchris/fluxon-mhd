package make_tempest_file;
use strict;
use warnings;
use Exporter qw(import);
our @EXPORT_OK = qw(make_tempest_file);
use PDL;
use PDL::NiceSlice;
use PDL qw(squeeze);
use PDL::GSL::INTERP;
use Math::Interpolator;
use PDL::Options;
use Math::Interpolate qw(linear_interpolate);
use Math::Interpolator::Linear;

=head2 map_fluxon_b

=for ref

Given a world and list of fluxons, generates a mapping of the representative magnetic field from these fluxons.

=cut

sub interpolate_1d_array {
    my ($x, $y, $x_new) = @_;

    # print pdl($x_new) . "\n";
    my $x_new_pdl = pdl($x_new);
    my @y_new;

    my $loops = $x_new_pdl->nelem;
    print $loops . "\n";

    for my $i (0..$loops) {

        push @y_new, linear_interpolate($x_new_pdl->(($i)), $x, $y);

    }
    print pdl(\@y_new) . "\n";
    return @y_new;
}




sub make_tempest_file {
    # """This is the function that should produce readable output for the Tempest code.
    # It should be able to take in a list of fluxons and output a file that can be read
    # by the Tempest code."""

    my $output_filename = shift;
    my $fluxon_list = shift;

    my @fluxons = @{$fluxon_list};

    die "No fluxons found. Cannot proceed without data." if scalar(@fluxons) == 0;
    print "Writing Tempest file to $output_filename\n";
    # Open an output file
    # if the file exists, delete it
    if (-e $output_filename) {
        unlink $output_filename;
    }

    open(my $fh, '>', $output_filename) or die "Could not open file '$output_filename': $!";

    my $first_ind = 1;

    # Loop through open fluxons and generate magnetic field profiles
    for my $fluxon_id (0..scalar(@fluxons)-1) {
        my $fluxon = $fluxons[$fluxon_id];

        # Check if the fluxon represents open fieldlines
        my $start_open = ($fluxon->{fc_start}->{label} == -1);
        my $end_open = ($fluxon->{fc_end}->{label} == -2);

        # Skip processing if it's a plasmoid or not an open fieldline
        if ($fluxon->{plasmoid} || ($start_open + $end_open != 1)) {
            next;
        }
        my $r0 = my $solar_radius = 1; #696e6;
        my $x, my $y, my $z;
        my $radius, my $radius_vec;
        my $theta, my $phi;
        my $bmag, my $bmag_2;
        my $magnetic_area;

        # Extract coordinates
        if ($end_open) {
            $x = squeeze($fluxon->dump_vecs->(0));
            $y = squeeze($fluxon->dump_vecs->(1));
            $z = squeeze($fluxon->dump_vecs->(2));
            $magnetic_area = pdl(map {$_->{A}} ($fluxon->vertices)) * $r0 * $r0;
            $magnetic_area(-1) .= $magnetic_area(-2);
        } else {
            $x = squeeze($fluxon->dump_vecs->(0,-1:0:-1));
            $y = squeeze($fluxon->dump_vecs->(1,-1:0:-1));
            $z = squeeze($fluxon->dump_vecs->(2,-1:0:-1));
            $magnetic_area = pdl(map { $_->{A}} ($fluxon->vertices))->(-1:0:-1) * $r0 * $r0;
            $magnetic_area(0) .= $magnetic_area(1);
        }
        $radius_vec = ($x**2 + $y**2 + $z**2)->sqrt * $solar_radius;
        $radius = 0.5 * $radius_vec->range([[0],[1]],[$radius_vec->dim(0)],'e')->sumover;
        $theta = acos($z/$radius_vec*$solar_radius);
        $phi = atan2($y, $x);

        # Calculate magnetic field magnitude
        my $bfield = $fluxon->bfield;
        # $bmag = squeeze(sqrt($bfield(0,:)**2 + $bfield(1,:)**2 + $bfield(2,:)**2));
        my $bmag_coord = squeeze(sqrt($bfield(0,:)**2 + $bfield(1,:)**2 + $bfield(2,:)**2));
        my $bmag_value = squeeze(sqrt($bfield(3,:)**2 + $bfield(4,:)**2 + $bfield(5,:)**2));

        # if (!$end_open) {
        #     $bmag = $bmag->(-1:0:-1,:);
        #     # $radius = $radius->(-1:0:-1);
        #     # $bmag_2 = $bmag_2->(-1:0:-1,:);
        # }


        if (!$end_open) {
            $bmag_coord = $bmag_coord->(-1:0:-1,:);
            $bmag_value = $bmag_value->(-1:0:-1,:);
        }


        # Create a pdl radius array with logarithmic spacing from 1 to 21.5
        my @bmag_array = $bmag_value->list;
        # my @radius_array = $radius->list;
        my @radius_array = $bmag_coord->list;
        my @radius_common_log = pdl(exp(sequence(102)/100*3.0445));
        my @bmag_common_log = interpolate_1d_array(\@radius_array, \@bmag_array, \@radius_common_log);

        # Interpolate values for test points





        # my $ipl = Math::Interpolator::Linear->new(@radius_array, @bmag_array);

        # $y = $ipl->y($x);
        # $x = $ipl->x($y);


        print "radius: ". pdl(\@radius_array) . "\n";
        print "bmag: ". pdl(\@bmag_array) . "\n\n";
        print "radius common log: ". pdl(\@radius_common_log) . "\n\n";
        print "bmag common log: ". pdl(\@bmag_common_log) . "\n";
        # Print outputs for verification


        if ($first_ind == 1) {
            our $model_count = 1;
            our $model_length = scalar(@radius_common_log);
            printf "%05d %05d\n", $model_count, $model_length;
            printf $fh "%05d %05d\n", $model_count, $model_length;
            $first_ind = 0;

            for my $i (0..$radius->getdim(0)-1) {
                printf $fh "%.8f\n", radius_common_log->at(($i))
            }
        }
        for my $i (0..$radius->getdim(0)-1) {
            printf $fh "%.8f\n", $bmag($i)->sclr
        }
    }

    # Close output file
    close $fh;
    print "Finished writing Tempest file.\n";

    # # Call a Python script
    # system("/opt/homebrew/anaconda3/envs/fluxenv/bin/python /Users/cgilbert/vscode/fluxons/fluxon-mhd/fluxpipe/fluxpipe/plotting/plot_bmag_all.py");


}

1;
