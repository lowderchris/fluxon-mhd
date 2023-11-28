package map_fluxon_b_all;
use strict;
use warnings;
use Exporter qw(import);
our @EXPORT_OK = qw(map_fluxon_b_all);
use PDL::NiceSlice;
use PDL qw(squeeze);

=head2 map_fluxon_b

=for ref

Given a world and list of fluxons, generates a mapping of the representative magnetic field from these fluxons.

=cut
use PDL::Options;

sub map_fluxon_b_all {
    my $output_filename = shift;
    my $fluxon_list = shift;

    my @fluxons = @{$fluxon_list};

    die "No fluxons found. Cannot proceed without data." if scalar(@fluxons) == 0;

    # Open an output file
    open(my $fh, '>', $output_filename) or die "Could not open file '$output_filename': $!";

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
        my $r0 = 696e6;
        my $solar_radius = 696e6;
        my $x, my $y, my $z;
        my $radius, my $radius_vec;
        my $theta, my $phi;
        my $magnetic_magnitude;
        my $A;

        # Extract coordinates
        if ($end_open) {
            $x = squeeze($fluxon->dump_vecs->(0));
            $y = squeeze($fluxon->dump_vecs->(1));
            $z = squeeze($fluxon->dump_vecs->(2));
            $A = pdl(map {$_->{A}} ($fluxon->vertices)) * $r0 * $r0;
            $A(-1) .= $A(-2);
        } else {
            $x = squeeze($fluxon->dump_vecs->(0,-1:0:-1));
            $y = squeeze($fluxon->dump_vecs->(1,-1:0:-1));
            $z = squeeze($fluxon->dump_vecs->(2,-1:0:-1));
            $A = pdl(map { $_->{A}} ($fluxon->vertices))->(-1:0:-1) * $r0 * $r0;
            $A(0) .= $A(1);
        }
        $radius_vec = ($x**2 + $y**2 + $z**2)->sqrt * $solar_radius;
        $radius = 0.5 * $radius_vec->range([[0],[1]],[$radius_vec->dim(0)],'e')->sumover;
        $theta = acos($z/$radius_vec*$solar_radius);
        $phi = atan2($y, $x);

        # Calculate magnetic field magnitude
        my $magnetic_field_0 = $fluxon->bfield;
        $magnetic_magnitude = squeeze(sqrt($magnetic_field_0(0,:)**2 + $magnetic_field_0(1,:)**2 + $magnetic_field_0(2,:)**2));
        # if ($magnetic_magnitude(0,0) > $magnetic_magnitude(-1,0)) {
        if ($end_open) {
            $magnetic_magnitude = $magnetic_magnitude->(-1:0:-1,:);
        }

        for my $i (0..$radius->getdim(0)-1) {
            printf $fh "%05d %.8f %.8f %.8f %.8e %.8f %.8f %.8e %.8e\n", $fluxon_id, $x($i)->sclr,
            $y($i)->sclr, $z($i)->sclr, $radius($i)->sclr, $theta($i)->sclr, $phi($i)->sclr, $magnetic_magnitude($i)->sclr, $A($i)->sclr;
        }
    }

    # Close output file
    close $fh;


    # Call a Python script
    system("/opt/homebrew/anaconda3/envs/fluxenv/bin/python /Users/cgilbert/vscode/fluxons/fluxon-mhd/fluxpipe/fluxpipe/plotting/plot_bmag_all.py");


}

1;
