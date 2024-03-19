package map_fluxon_b;
use strict;
use warnings;
use Exporter qw(import);
our @EXPORT_OK = qw(map_fluxon_b);
use PDL::NiceSlice;
use PDL qw(squeeze);

=head2 map_fluxon_b

=for ref

Given a world and list of fluxons, generates a mapping of the representative magnetic field from these fluxons.

=cut
use PDL::Options;

sub map_fluxon_b {
    my $output_filename = shift;
    my $fluxon_list = shift;
    my @fluxons = @{$fluxon_list};

    # Initialize storage arrays
    my @fluxon_positions = ();
    my @beginning_phis = ();
    my @beginning_thetas = ();
    my @ending_phis = ();
    my @ending_thetas = ();
    my @beginning_magnetic_fields = ();
    my @ending_magnetic_fields = ();
    my @area_beginning = ();
    my @area_ending = ();

    # Loop through open fluxons
    for my $fluxon_id (0..scalar(@fluxons)-1) {
        my $fluxon = $fluxons[$fluxon_id];

        # Check if the fluxon represents open fieldlines
        my $start_open = ($fluxon->{fc_start}->{label} == -1);
        my $end_open = ($fluxon->{fc_end}->{label} == -2);

        # Skip processing if it's a plasmoid or not an open fieldline
        if ($fluxon->{plasmoid} || ($start_open + $end_open != 1)) {
            next;
        }

        my $solar_radius = 696e6;
        my $radius;
        my $theta;
        my $phi;
        my $magnetic_area;
        my $x, my $y, my $z;
        my $radius_vec;

        # Extract coordinates
        if ($end_open) {
            $x = squeeze($fluxon->dump_vecs->(0));
            $y = squeeze($fluxon->dump_vecs->(1));
            $z = squeeze($fluxon->dump_vecs->(2));
            $magnetic_area = pdl(map {$_->{A}} ($fluxon->vertices)) * $solar_radius * $solar_radius;
            $magnetic_area(-1) .= $magnetic_area(-2);
        } else {
            $x = squeeze($fluxon->dump_vecs->(0,-1:0:-1));
            $y = squeeze($fluxon->dump_vecs->(1,-1:0:-1));
            $z = squeeze($fluxon->dump_vecs->(2,-1:0:-1));

            $magnetic_area = pdl(map {$_->{A}} ($fluxon->vertices))->(-1:0:-1) * $solar_radius * $solar_radius;
            $magnetic_area(0) .= $magnetic_area(1);
        }

        $radius_vec = ($x**2 + $y**2 + $z**2)->sqrt * $solar_radius;
        $radius = 0.5 * $radius_vec->range([[0],[1]],[$radius_vec->dim(0)],'e')->sumover;
        $theta = acos($z/$radius_vec*$solar_radius);
        $phi = atan2($y, $x);




        # Calculate magnetic field magnitude
        my $bfield = $fluxon->bfield();
        my $bmag = cat squeeze(sqrt($bfield(0,:)**2 + $bfield(1,:)**2 + $bfield(2,:)**2)),
                       squeeze(sqrt($bfield(3,:)**2 + $bfield(4,:)**2 + $bfield(5,:)**2));


        if ($end_open) {
            $bmag = $bmag->(-1:0:-1,:);
        }

        ## First / last two points without b-field
        $bmag = $bmag->slice('2:-2,:');

        # for my $id (0..5) {
        #     # printf "\n bmag: %6f,", $bmag($id)->sclr;
        #     printf "\n b0: %.6f, b1: %6f, b2 %6f", $bfield(0, $id)->sclr, $bfield(1, $id)->sclr, $bfield(2, $id)->sclr;
        #     printf " b3: %.6f, b4: %6f, b5 %6f", $bfield(3, $id)->sclr, $bfield(4, $id)->sclr, $bfield(5, $id)->sclr;
        #     print "\n";
        # }


        my $magnetic_beginning = squeeze($bmag(0,1));
        my $magnetic_middle = squeeze($bmag($bmag->shape->(0)/2,1));
        my $magnetic_ending = squeeze($bmag(-1,1));

        # Append values to storage arrays
        push(@fluxon_positions, $fluxon_id);
        push(@beginning_phis, squeeze($phi(0)));
        push(@beginning_thetas, squeeze($theta(0)));
        push(@ending_phis, squeeze($phi(-1)));
        push(@ending_thetas, squeeze($theta(-1)));
        push(@beginning_magnetic_fields, $magnetic_beginning);
        push(@ending_magnetic_fields, $magnetic_ending);
        push(@area_beginning, squeeze($magnetic_area(0)));
        push(@area_ending, squeeze($magnetic_area(-1)));
    }

    # Output data to disk
    wcols pdl(@fluxon_positions), pdl(@beginning_phis), pdl(@beginning_thetas), pdl(@ending_phis), pdl(@ending_thetas), squeeze(pdl(@beginning_magnetic_fields)), squeeze(pdl(@ending_magnetic_fields)), squeeze(pdl(@area_beginning)), squeeze(pdl(@area_ending)), $output_filename;

    system("/opt/homebrew/anaconda3/envs/fluxenv/bin/python /Users/cgilbert/vscode/fluxons/fluxon-mhd/fluxpipe/fluxpipe/plotting/plot_bmag.py");

}

__END__
