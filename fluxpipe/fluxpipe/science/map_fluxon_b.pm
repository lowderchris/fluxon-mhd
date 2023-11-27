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

    # Loop through open fluxons and generate wind profiles
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

        # Extract coordinates
        if ($end_open) {
            my $x = squeeze($fluxon->dump_vecs->(0));
            my $y = squeeze($fluxon->dump_vecs->(1));
            my $z = squeeze($fluxon->dump_vecs->(2));
            my $radius_vec = ($x**2 + $y**2 + $z**2)->sqrt * $solar_radius;
            $radius = 0.5 * $radius_vec->range([[0],[1]],[$radius_vec->dim(0)],'e')->sumover;
            $theta = acos($z/$radius_vec*$solar_radius);
            $phi = atan2($y, $x);
            $magnetic_area = pdl(map {$_->{A}} ($fluxon->vertices)) * $solar_radius * $solar_radius;
            $magnetic_area(-1) .= $magnetic_area(-2);
        } else {
            my $x = squeeze($fluxon->dump_vecs->(0,-1:0:-1));
            my $y = squeeze($fluxon->dump_vecs->(1,-1:0:-1));
            my $z = squeeze($fluxon->dump_vecs->(2,-1:0:-1));
            my $radius_vec = ($x**2 + $y**2 + $z**2)->sqrt * $solar_radius;
            $radius = 0.5 * $radius_vec->range([[0],[1]],[$radius_vec->dim(0)],'e')->sumover;
            $theta = acos($z/$radius_vec*$solar_radius);
            $phi = atan2($y, $x);
            $magnetic_area = pdl(map {$_->{A}} ($fluxon->vertices))->(-1:0:-1) * $solar_radius * $solar_radius;
            $magnetic_area(0) .= $magnetic_area(1);
        }

        my $magnetic_field_0 = $fluxon->bfield();
        my $magnetic_magnitude = cat squeeze(sqrt($magnetic_field_0(0,:)**2 + $magnetic_field_0(1,:)**2 + $magnetic_field_0(2,:)**2)), squeeze(sqrt($magnetic_field_0(3,:)**2 + $magnetic_field_0(4,:)**2 + $magnetic_field_0(5,:)**2));

        if ($magnetic_magnitude(0,0) > $magnetic_magnitude(-1,0)) {
            $magnetic_magnitude = $magnetic_magnitude->(-1:0:-1,:);
        }

        $magnetic_magnitude = $magnetic_magnitude->slice('2:-2,:');

        my $magnetic_beginning = squeeze($magnetic_magnitude(0,1));
        my $magnetic_middle = squeeze($magnetic_magnitude($magnetic_magnitude->shape->(0)/2,1));
        my $magnetic_ending = squeeze($magnetic_magnitude(-1,1));

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
}

__END__
