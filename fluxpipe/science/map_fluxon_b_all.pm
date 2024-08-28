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
    printf $fh "fnum \t\t x \t\t y \t\t z \t\t\t radius \t\t theta \t\t phi \t\t A \t\t\t b_mag \t\t b_mag2 \t\t\t b_x \t\t b_y \t\t b_z\n";

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
        my $bmag_coord, my $bmag_value;
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
        $bmag_coord = squeeze(sqrt($bfield(0,:)**2 + $bfield(1,:)**2 + $bfield(2,:)**2));
        $bmag_value = squeeze(sqrt($bfield(3,:)**2 + $bfield(4,:)**2 + $bfield(5,:)**2));

        if (!$end_open) {
            $bmag_coord = $bmag_coord->(-1:0:-1,:);
            $bmag_value = $bmag_value->(-1:0:-1,:);
        }

        # Fix the bottom magnetic field vertex
        $bmag_value->((1)) .= $bmag_value->((2)) * $magnetic_area->((2)) / $magnetic_area->((1));
        $bmag_value->((0)) .= $bmag_value->((1)) * $magnetic_area->((1)) / $magnetic_area->((0));

        # print("\n\nmap_fluxon_b_all.pm: A:\n");
        # print $magnetic_area . "\n\n";
        # print "Bmag2:" . $bmag_value;
        # print "\n";
        # <STDIN>;

        for my $i (0..$radius->getdim(0)-1) {
            printf $fh "%05d %.8f %.8f %.8f %.8e %.8f %.8f %.8e %.8e %.8e %.8e %.8e %.8e\n", $fluxon_id,
            $x($i)->sclr, $y($i)->sclr, $z($i)->sclr,
            $radius($i)->sclr, $theta($i)->sclr, $phi($i)->sclr,
            $magnetic_area($i)->sclr, $bmag_coord($i)->sclr, $bmag_value($i)->sclr,
            $bfield(0,$i)->sclr, $bfield(1,$i)->sclr, $bfield(2,$i)->sclr;

        }

    }
    # Close output file
    close $fh;


    # # Call a Python script
    # system("/opt/homebrew/anaconda3/envs/fluxenv/bin/python /Users/cgilbert/vscode/fluxons/fluxon-mhd/fluxpipe/fluxpipe/plotting/plot_bmag_all.py");


}

1;
