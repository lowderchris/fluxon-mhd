=head1 NAME

gen_fluxon_wsaflow - Comprehensive Module for Solar Wind Solution and Fluxon Flow Analysis

=cut

package gen_fluxon_ghostsflow;
use strict;
use warnings;
use Exporter qw(import);
our @EXPORT_OK = qw(gen_fluxon_ghostsflow);

use PDL::Graphics::Gnuplot;
use PDL;
use PDL::NiceSlice;
use PDL::Options;
use PDL::ImageND;
use PDL ('which_min'); # This imports which_min directly into your current package
use POSIX;
use Math::RungeKutta;
use PDL::GSL::INTEG;
use PDL qw(squeeze);
use PDL::GSL::INTERP;

$PDL::verbose = 0;

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 USAGE


=head1 DEPENDENCIES

This module depends on the following Perl modules:
=over 4

=back

=head1 DIAGNOSTICS

Set the C<$verbose> flag to 1 for diagnostic output.

=head1 AUTHOR

Gilly <gilly@swri.org> (and others!)

=head1 SEE ALSO

Other relevant modules and documentation.

=cut


our $FIRST = 1;
my $verbose = 0;

# Interpolate an array given X and Y values along with a test point
sub interpolate_1d {
    my ($x, $y, $x0) = @_;
    my $y0;
    my $len = $x->nelem;
    my $ind = abs($x - $x0)->minimum_ind;
    if ($ind == 0) {
        $y0 = $y->(($ind));
    } elsif ($ind == $len-1) {
        $y0 = $y->(($ind));
    } else {
        my $x1 = $x->(($ind));
        my $x2 = $x->(($ind+1));
        my $y1 = $y->(($ind));
        my $y2 = $y->(($ind+1));
        $y0 = $y1 + ($y2 - $y1) * ($x0 - $x1) / ($x2 - $x1);
    }
    return $y0;
}

sub gen_fluxon_ghostsflow {
    my $fluxon = shift;
    my $fid = shift;

    # Check if the start and end    points of the fluxon are open
    my $is_start_open = ($fluxon->{fc_start}->{label} == -1);
    my $is_end_open = ($fluxon->{fc_end}->{label} == -2);

    # If the fluxon is labeled as a plasmoid, or if both ends are open, return undefined
    return undef if($fluxon->{plasmoid} || ($is_start_open + $is_end_open != 1));


    (my $r_vr_scaled, my $r_fr_scaled, my $theta, my $phis) =
            gen_fluxon_ghostsflow_physics($fluxon, $fid);

    # Return the final fluxon array, fluxon radius, and magnetic field components
    return ($r_vr_scaled, $r_fr_scaled, $theta, $phis);
}

sub read_data {
    my ($filename) = @_;
    # Open the file or die with an error message
    open my $fh, '<', $filename or die "Cannot open file $filename: $!";

    # Skip the header lines
    for (1..10) {
        <$fh>;
    }

    # Initialize arrays to store the data
    my (@r_Rsun, @rho, @u, @Valf, @T);

    # Read the file line by line
    while (my $line = <$fh>) {
        # Split the line into columns
        my ($dummy, $r_Rsun_val, $rho_val, $u_val, $Valf_val, $T_val) = split /\s+/, $line;

        # Store the data in arrays
        push @r_Rsun, $r_Rsun_val;
        push @rho, $rho_val;
        push @u, $u_val;
        push @Valf, $Valf_val;
        push @T, $T_val;
    }
    # Close the file handle
    close $fh;

    # Convert the arrays to PDL objects for numerical operations
    my $r_Rsun_pdl = pdl(\@r_Rsun);
    my $rho_pdl = pdl(\@rho);
    my $u_pdl = pdl(\@u) / 1e5;
    my $Valf_pdl = pdl(\@Valf) / 1e5;
    my $T_pdl = pdl(\@T);

    # Return the PDL objects
    return ($r_Rsun_pdl, $rho_pdl, $u_pdl, $Valf_pdl, $T_pdl);
}



# Function to plot data
sub plot_data {
    my ($r_Rsun, $rho, $u, $Valf, $T) = @_;
    $r_Rsun = $r_Rsun -1;
    my $plot = PDL::Graphics::Gnuplot->new(persist => 1);

    # Plot mass density
    $plot->plot({title=>"Polar Coronal Hole Background Properties", logscale=>'xy', xlabel=>'r/Rsun', ylabel=>'Values'}, {legend=>"Mass Density (g/cm^3)"}, $r_Rsun, $rho);
    $plot->replot({legend=>"Mass Density (g/cm^3)", with =>'points'}, $r_Rsun, $rho);

    # Plot radial wind flow speed
    $plot->replot({legend=>"Radial Wind Flow Speed (km/s)"}, $r_Rsun, $u);
    $plot->replot({legend=>"Radial Wind Flow Speed (km/s)", with =>'points'}, $r_Rsun, $u);

    # Plot Alfven speed
    $plot->replot({legend=>"Alfven Speed (km/s)"}, $r_Rsun, $Valf);
    $plot->replot({legend=>"Alfven Speed (km/s)", with =>'points'}, $r_Rsun, $Valf);

    # Plot temperature
    $plot->replot({legend=>"Temperature (K)"}, $r_Rsun, $T);
    $plot->replot({legend=>"Temperature (K)", with =>'points'}, $r_Rsun, $T);

    <STDIN>;
}


sub interpolate_1d_array {
    my ($x, $y, $x_new) = @_;
    my $interp = PDL::GSL::INTERP->new("linear", $x, $y);
    my $y_new = $interp->interp($x_new);
    return $y_new;
}

sub gen_fluxon_ghostsflow_physics {
    my $me = shift;
    my $fid = shift;
    my $u_opt = shift // {};


    # Define and append optional option input
    my %opt = parse( {
        steps => 500,
        g0 => 280,      # m/s^2
        r0 => 696e6,    # m
        #v0 => 10,       # km/sec
        v0 => rand(10)+10,       # km/sec
        cs => 180,      # km/sec at 2e6K and fully ionized protons
                        },
                        $u_opt
        );

    # Pass options into local variables
    my $g0 = $opt{g0};
    my $r0 = $opt{r0};
    my $cs = $opt{cs} * 1e3; ## convert to m/s from km/s

    # Check for start and ending points being open
    our $en_open = ($me->{fc_end}->{label}==-2);

    # Calculate array of sphiserical coordinate positions and areas along the fluxon
    # Work along the correct direction depending on which end is open
    if($en_open) {
        my $x = squeeze($me->dump_vecs->(0));
        my $y = squeeze($me->dump_vecs->(1));
        my $z = squeeze($me->dump_vecs->(2));
        our $r1 = ($x**2 + $y**2 + $z**2)->sqrt * $opt{r0};
        our $r = 0.5 * $r1->range([[0],[1]],[$r1->dim(0)],'e')->sumover;
        our $thetas = acos($z/$r1*$opt{r0});
        our $phis = atan2($y, $x);
        our $A = pdl(map {$_->{A}} ($me->vertices));
        our $bfield = pdl(map {$_->{b_vec}} ($me->vertices));
        our $bmag = pdl(map {$_->{b_mag}} ($me->vertices));
        our $T = pdl(map {$_->{T}} ($me->vertices));
        our $rho = pdl(map {$_->{rho}} ($me->vertices));

    } else {
        my $x = squeeze($me->dump_vecs->(0,-1:0:-1));
        my $y = squeeze($me->dump_vecs->(1,-1:0:-1));
        my $z = squeeze($me->dump_vecs->(2,-1:0:-1));
        our $r1 = ($x**2 + $y**2 + $z**2)->sqrt * $opt{r0};
        our $r = 0.5 * $r1->range([[0],[1]],[$r1->dim(0)],'e')->sumover;
        our $thetas = acos($z/$r1*$opt{r0});
        our $phis = atan2($y, $x);
        our $A = pdl(map { $_->{A}} ($me->vertices))->(-1:0:-1);
        our $bfield = pdl(map {$_->{b_vec}} ($me->vertices))->(-1:0:-1);
        our $bmag = pdl(map {$_->{b_mag}} ($me->vertices))->(-1:0:-1);
        our $T = pdl(map {$_->{T}} ($me->vertices))->(-1:0:-1);
        our $rho = pdl(map {$_->{rho}} ($me->vertices))->(-1:0:-1);

    }

    our $A;
    our $bfield;
    our $bmag;
    our $r1;
    our $r;
    our $thetas;
    our $phis;
    our $T;
    our $rho;
    my $rn = $r1/$opt{r0};

    # Fix the bottom r vertex
    my $diff0 = abs(1 - $rn->((0)));
    $rn(0) .= 1 + $diff0;

    my $zn = $rn - 1;
    my $len = $rn->nelem - 1;


    # Get rid of end anomalies
    # $A->((0)) .= $A->((1));
    $A->((-1)) .= $A->((-2));

    # Fix the bottom magnetic field vertex
    $bmag->((1)) .= $bmag->((2)) * $A->((2)) / $A->((1));
    $bmag->((0)) .= $bmag->((1)) * $A->((1)) / $A->((0));

    ## Declare Variables
    my $vnow;
    my $dr;
    my @rv;
    my $rnow;
    my @num;
    my @denom;
    my $step;
    my $vint;
    my $abserr;
    my $ierr;

    my $filename = "fluxon-data/gilly_background_cvb07.dat";
    # my $filename = '/Users/cgilbert/vscode/fluxons/fluxon-data/gilly_background_cvb07.dat';

    my ($r_zeph, $rho_zeph, $u_zeph, $v_alph_zeph, $T_zeph) = read_data($filename);


    # my $rho_interp  = $rho_zeph->   interpolate($rn, $r_zeph);
    my ($rho_interp , $err1) =  interpolate($rn, $r_zeph, $rho_zeph);
    my ($u_interp   , $err2) =  interpolate($rn, $r_zeph, $u_zeph);
    my ($Valf_interp, $err3) =  interpolate($rn, $r_zeph, $v_alph_zeph);
    my ($T_interp   , $err4) =  interpolate($rn, $r_zeph, $T_zeph);

    # plot_data($r_zeph, $rho_zeph, $u_zeph, $v_alph_zeph, $T_zeph);
    # plot_data($rn, $rho_interp, $u_interp, $Valf_interp, $T_interp);



    # FIND THE SUN SURFACE
    my $first_ind = 1;
    my $r_sun = $rn->(($first_ind));
    my $B_sun = $bmag->(($first_ind));
    my $A_sun = $A->(($first_ind));

    my $phisi0 = $phis->(($first_ind));
    my $sin_theta0 = sin($thetas)->(($first_ind));
    my $r00 = $rn->(($first_ind));

    # FIND THE SOURCE SURFACE
    my $r_ss = 2.5; # Solar radii
    my $B_ss = interpolate_1d($rn, $bmag, $r_ss);
    my $A_ss = interpolate_1d($rn, $A, $r_ss);
    my $f_ss = abs($A_ss / $A_sun) * ($r_sun * $r_sun) / ($r_ss * $r_ss);

    # my $f_all= abs($B_sun / $bmag) * ($r_sun * $r_sun) /   ($rn * $rn);
    my $f_all =    abs($A / $A_sun) * ($r_sun * $r_sun) /   ($rn * $rn);

    # FIND THE TOP OF THE DOMAIN
    my $last_ind    = $phis->nelem - 1;
    my $phisi1      = $phis->(($last_ind));
    my $r11         = $rn->(($last_ind));
    my $sin_theta1  = sin($thetas)->(($last_ind));

    if ($phisi0 < 0) {$phisi0 += 2 * 3.1415926;}
    if ($phisi1 < 0) {$phisi1 += 2 * 3.1415926;}


    sub densfac_func{
        my $B_sun = shift;
        my $zn = shift;

        my $D0 = 0.85;
        my $D_ifty = ($B_sun / 6.9)**(1.18); #B in gauss from 3 - 50 gauss or so

        my $D = $D0 + 0.5*($D_ifty - $D0) * tanh(($zn-0.37)/0.26);

        return $D;
    }

    my $densfac = densfac_func($B_sun, $zn);
    my $sqrt_densfac = $densfac**(0.5);
    my $u_scaled = 2.5*$u_interp/ $sqrt_densfac;

    my $u_top = $u_scaled(-1);
    # print("Densfac = $sqrt_densfac\n");
    # print("u_scaled = $u_scaled\n");
    # print("B_sun = $B_sun, u_top = $u_top\n");

    # my $distance_degrees = interpolate_2d_lonlat($distance_array_degrees, $phisi0, $sin_theta0);

    # Calculate the wind speed at the given location
    # my $flow_field = wind_speed($An, $distance_degrees);
    # my $speed = interpolate_1d($rn, $flow_field, 5.0);

    # my $speed = wind_speed($A_ss, $distance_degrees);

    # write the $fid, $f_ss, and $distance_degrees to a file
    # open my $fh, '>>', 'seq_a ss_theta.csv' or di e "Cannot open data.cs v: $!";
    # print $fh "$fluxon_id, $A_ss, $distance_degrees, $speed\n";
    # close $fh;

    # print "\nfss = $f_ss\n";
    # print "Distance = $dista nce_degrees\n";
     # print "Flow Field =  $flow_field\n\n";

    use strict;
    use warnings;
    use PDL;
    use PDL::Graphics::Gnuplot;

    # Sample data
    my $ones = ones($r->dims);
    # print "\nghostsflow bmag:\n";
    # print $bmag;
    # print "\n\n";

    # <STDIN>;

    if (0) {
        # # Plot the data with labels
        my $plot = PDL::Graphics::Gnuplot->new(persist => 1);
        $plot->plot({title=>"Fluxon #$fid", logscale=>'xy', xlabel=>'z = r/R - 1', font=>'48'}, {legend=>"|B|"},$zn, $bmag);
        $plot->replot({legend=>"|B|", with =>'points'},$zn, $bmag);

        # $plot->replot({legend=>"bmag"},$zn, $bmag);
        # $plot->replot({legend=>"f_ss"},$rn, $ones * $f_ss);
        # $plot->replot({legend=>"f_all"} ,$rn, $f_all);
        $plot->replot({legend=>"f_r"},$zn, $f_all);
        $plot->replot({legend=>"f_r", with=>"points"},$zn, $f_all);

        $plot->replot({legend=>"A"},$zn, $A);
        $plot->replot({legend=>"A", with=>"points"},$zn, $A);

        $plot->replot({legend=>"d_f"},$zn, $densfac);
        $plot->replot({legend=>"d_f", with=>"points"},$zn, $densfac);

        $plot->replot({legend=>"u_w"},$zn, $u_scaled);
        $plot->replot({legend=>"u_w", with=>"points"},$zn, $u_scaled);

        # $plot->replot({legend=>"flow", with=>"points"},$rn->((-1)), $flow_field->((-1)));
        # $plot->replot({legend=>"phis", with=>"points"},$zn, $phis);

        # $plot->replot(legend=>"phis0", with=>"points", pointsize=>3, $r00, $phisi0);
        # $plot->replot(legend=>"phis1", with=>"points", pointsize=>3, $r11, $phisi1);

        # $plot->replot({legend=>"stheta"},$zn, sin($thetas));
        # $plot->replot({legend=>"stheta", with=>"points"}, $zn, sin($thetas));

        # $plot->replot(legend=>"sth0", with=>"points", pointsize=>3 ,$r00, $sin_theta0);
        # $plot->replot(legend=>"sth0", with=>"points", pointsize=>3 ,$r11, $sin_theta1);
        # $plot->replot({legend=>"distance"},$rn, $ones * $distance_degrees);

        # my $plot = PDL::Graphics::Gnuplot->new(persist => 1);
        # my $top_speed = $flow_field->((-1));
        # $plot->plot({title=>"ID: $fluxon_id, Top Speed: $top_speed", logscale => 'xy'}, {legend=>"rn", with=>"points"},$rn, $rn);

        # Block execution
        print "Press ENTER to continue...\n";
        <STDIN>;
    }

    ## Calculate Return Values ##
    # my $speed_tall = $u_interp / $densfac**(0.5);
    # This array is (r, v) in units of (r_sun, km/s)

    # our $r_v_scaled = pdl($rn, $ones * ($flow_field->at($len)));
    our $r_v_scaled = pdl($rn, $u_scaled);
    # if (1-$en_open) {
    #     # print $r_v_scaled;
    # } else {
    #     our $r_v_scaled = pdl($rn, $ones * $flow_field->at(0));
    #     # print $r_v_scaled;
    # }

    # our $r_v_scaled;
    # This array is (r, fr) in units of (r_sun, unitless)
    my $r_fr_scaled = pdl($rn, $f_all);
    # print $f_all;

    # Return the constructed arrays
    return ($r_v_scaled, $r_fr_scaled, $thetas, $phis);
}


use PDL;
use PDL::NiceSlice;

sub gradient_descent {
    my ($image, $latitude, $longitude) = @_;
    my $threshold = 0.1;

    # Convert geographisic coordinates to image coordinates
    my $row = ($latitude);   # Assuming latitude and longitude are already in pixel coordinates
    my $col = ($longitude);

    # Check initial value
    my $initial_value = $image->at($row, $col);

    # Determine whether to perform descent or ascent
    my $is_ascent = $initial_value < $threshold;

    # Gradient descent or ascent
    my $max_steps = 500;  # Maximum number of steps to prevent infinite loops
    my $step_size = 1;    # Initial step size in pixels

    for (my $step = 0; $step < $max_steps; $step++) {
        my ($best_value, $best_row, $best_col) = ($initial_value, $row, $col);

        foreach my $r_offset (-$step_size, 0, $step_size) {
            foreach my $c_offset (-$step_size, 0, $step_size) {
                next if $r_offset == 0 && $c_offset == 0;  # Skip the current point

                my $new_row = $row + $r_offset;
                my $new_col = $col + $c_offset;

                # Ensure new coordinates are within bounds
                next if $new_row < 0 || $new_row >= $image->dim(0);
                next if $new_col < 0 || $new_col >= $image->dim(1);

                my $new_value = $image->at($new_row, $new_col);

                # Update best value found depending on descent or ascent
                if (($is_ascent && $new_value > $best_value) || (!$is_ascent && $new_value < $best_value)) {
                    ($best_value, $best_row, $best_col) = ($new_value, $new_row, $new_col);
                }
            }
        }

        # # Dynamically adjust the step size
        # if (abs($best_value - $initial_value) < 0.01) {
        #     $step_size += 1;
        # } else {
        #     $step_size = 1;  # Reset step size to initial value
        # }

        # Check if the best value found crosses the threshold
        if ((!$is_ascent && $best_value < $threshold) || ($is_ascent && $best_value >= $threshold)) {
            my $distance = sqrt(($best_row - $row)**2 + ($best_col - $col)**2);
            return ($best_row, $best_col, $distance);
        }

        # Update row and col for the next iteration
        $initial_value = $best_value;
        ($row, $col) = ($best_row, $best_col);
    }

    # Return NaN if threshold not crossed within max steps
    return ('NaN', 'NaN', 'NaN');
}


use PDL;
use PDL::Graphics::Gnuplot;

sub plot_image_with_points {
    my ($image_path, $original_point, $discovered_point, $fluxon_id) = @_;
    our $fh;
    # print $fluxon_id;
    # Save data and image path to a CSV file
    if ($fluxon_id == 0){
        open $fh, '>', 'data.csv' or die "Cannot open data.csv: $!";
        print $fh "x1,y1,x2,y2,image_path=$image_path\n";
    } else {
        open $fh, '>>', 'data.csv' or die "Cannot open data.csv: $!";
    }
    print $fh "$original_point->[1],$original_point->[0],$discovered_point->[1],$discovered_point->[0]\n";
    close $fh;


    # print "Data saved to data.csv\n";


}

sub do_image_plot{
    # Call the Python script for plotting
    my $python_script = 'plot_distance.py';
    my $exit_status = system("python $python_script");
    if ($exit_status == 0) {
        print "Python script executed successfully.\n";
    } else {
        print "Python script failed to execute.\n";
    }

    # <STDIN>;
}
# our $count = 0;

1;
__END__


    # # # Define a normalized expansion factor
    # my $Arat = $A / $A->((1));
    # my $Rrat = ($r->((1))) * ($r->((1))) / $r / $r;
    # my $Brat = $bmag / $bmag->((1));
    # our $fr1 = $Arat * $Rrat;

    # # Initialize Gnuplot object
    # use PDL::Graphics::Gnuplot;

    # gpwin('pngcairo')->plot(
    #     { title => 'Field vs. r', xlabel => 'r', ylabel => 'Field Value', logscale => 'xy' },
    #     with => 'lines', $r1/$opt{r0}, $Arat,
    #     with => 'lines', $r1/$opt{r0}, $fr1,
    #     # with => 'lines', $r, $r,
    #     with => 'lines', $r1/$opt{r0}, $bmag,
    # );

    # my $ch_map = rfits($ch_map_path);
    # my ($theta_bound, $phisi_bound, $distance) = gradient_descent($ch_map, $th0_ind, $phis0_ind);

    # $distance = $distance * 3.14159265358 / 180;
    # plot_image_with_points($ch_map_path, [$th0_ind, $phis0_ind], [$theta_bound, $phisi_bound], $fluxon_id);