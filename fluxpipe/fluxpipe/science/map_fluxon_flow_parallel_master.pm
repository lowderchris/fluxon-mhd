=head1 NAME

Fluxon Flow Mapper - Parallelized Solar Wind Flow Mapping Along Fluxon Structures

=cut

package map_fluxon_flow_parallel_master;
use strict;
use warnings;
use Exporter qw(import);
our @EXPORT_OK = qw(map_fluxon_flow_parallel_master);
use PDL;
use PDL::NiceSlice;
use PDL::Options;
use Parallel::ForkManager;
use PDL::IO::Storable;
use File::Path qw(make_path);
use Time::HiRes qw(clock_gettime);
use Data::Dumper;
use Chart::Gnuplot;
use gen_fluxon_tflow qw(gen_fluxon_tflow);
use gen_fluxon_schonflow qw(gen_fluxon_schonflow);
use gen_fluxon_wsaflow qw(gen_fluxon_wsaflow do_image_plot);
use gen_fluxon_cranmerflow qw(gen_fluxon_cranmerflow);
use pipe_helper qw(configurations);


my %configs = configurations();

# Control Flags
$PDL::verbose = 0;
my $concurrency = $configs{concurrency} // 4;
my $temp_dir = "temp_dir";


=head1 DESCRIPTION

This Perl module provides functionality to map the solar wind flow along fluxon structures in a parallelized manner. It uses the Parallel::ForkManager library to manage concurrent processes and PDL (Perl Data Language) for numerical computations.

=head1 AUTHOR

Gilly <gilly@swri.org> and others!

=head1 LICENSE

This library is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut

=head2 map_fluxon_flow_parallel_master

=for ref

Given a world and list of fluxons, generates a mapping of the solar wind flow along these fluxon structures.

=cut

sub map_fluxon_flow_parallel_master {
    my $output_file_name = shift;
    my $fluxon_list = shift;
    my $distance_array_degrees = shift;
    my @fluxons = @{$fluxon_list};
    make_path($temp_dir) unless -d $temp_dir;
    my $max_processes = shift || $concurrency;    # Set the number of parallel processes
    my $highest_fluxon_done = 0;
    my $before_time = clock_gettime();

    my $n_choke = 12;
    my $choke = 0;
    my $do_plot_charts = 0;


    # # Create a simple plot window
    # my $win = PDL::Graphics::Simple->new();

    # # Plot the image
    # $win->imag($distance_array_degrees);

    # Save the image as img.png
    # $win->close("img.png");



    use PDL;
    use PDL::Primitive;
    use PDL::Func;


    sub interpolate_2d_lonlat {
        our ($image, $long_i, $latt_i) = @_;

        $image = pdl($image);
        # Define your image dimensions
        my ($img_width, $img_height) = $image->dims; # 900 by 360

        # Create grids for latitude and longitude
        my $long_vals = pdl(sequence($img_width)/($img_width-1) * 2 * 3.14159265);  # 0 to 2*pi
        my $latt_vals = pdl(sequence($img_height)/($img_height-1) * 2 - 1);  # -1 to 1

        my ($ind_long) = minimum_ind(abs($long_vals - $long_i));
        my ($ind_latt) = minimum_ind(abs($latt_vals - $latt_i));
        my $imval = $image->at($ind_long, $ind_latt);
        return $imval;
    }

    # Count the number of open and closed fluxons
    my $num_closed_fluxons = 0;
    my @open_fluxon_indices = ();
    my $num_open_fluxons = 0;
    for my $fluxon_index (0..scalar(@fluxons)-1) {
        my $fluxon = $fluxons[$fluxon_index];

        # Check for open fieldline, and skip if not
        my $start_open = ($fluxon->{fc_start}->{label} == -1);
        my $end_open = ($fluxon->{fc_end}->{label} == -2);
        if ($fluxon->{plasmoid} || ($start_open + $end_open != 1)) {
            $num_closed_fluxons++;
        } else {
            $num_open_fluxons++;
            push(@open_fluxon_indices, $fluxon_index);
        }
    }

    my $iteration_count = -1;
    my $fork_manager = Parallel::ForkManager->new($max_processes, $temp_dir);

    my @results = ();

    $fork_manager->run_on_finish(sub {
        $iteration_count++;
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $result) = @_;
        if ($exit_code == 0) {
            push(@results, $result);
            my $this_fluxon_id = $result->{fluxon_position} + 1;

            if ($this_fluxon_id > $highest_fluxon_done) {
                $highest_fluxon_done = $this_fluxon_id;

                my $after_time = clock_gettime();
                my $elapsed_time = $after_time - $before_time;
                my $rounded_elapsed_time = rint($elapsed_time * 10) / 10;
                my $remaining_fluxons = $num_open_fluxons - $highest_fluxon_done;
                my $time_each_fluxon = $elapsed_time / $highest_fluxon_done;
                my $time_remaining = $remaining_fluxons * $time_each_fluxon;
                my $rounded_remaining_time = rint($time_remaining * 10) / 10;
                my $rounded_remaining_minutes = rint($rounded_remaining_time / 6) / 10;

                print "\r", "\t\t\tCalculated fluxons: ", $highest_fluxon_done, ' of ', $num_open_fluxons, ", ", $rounded_elapsed_time, "(s) elapsed, ", $rounded_remaining_time, "(s) [$rounded_remaining_minutes mins] remaining......";
            }
        }
    });

    my $num_total_fluxons = scalar(@fluxons);
    my $num_active_fluxons = 2 * ($num_total_fluxons - $num_open_fluxons) + $num_open_fluxons;
    print "\n\t\tRunning with $max_processes Cores on $num_open_fluxons open Fluxons out of $num_total_fluxons total Fluxons from $num_active_fluxons footpoints!\n\n\t\tBeginning Calculation...\n\n";

    $PDL::verbose = 0;
    no warnings;

    my $max_fluxon_id = scalar(@open_fluxon_indices) - 1;
    if ($choke) {
        if ($n_choke < $max_fluxon_id) {
            $max_fluxon_id = $n_choke;
        }
    }

    my %configs = configurations();
    my $flow_method = $configs{flow_method} || "parker";

    print "\t\tThe flow method is ''$flow_method.''\n\n\n";


    for my $fluxon_id (0..$max_fluxon_id - 1) {
    # for my $fluxon_id (12..13 - 1) {
        $fork_manager->start and next;

        my $fluxon = $fluxons[$fluxon_id];
        my ($r_vr_scaled, $r_fr_scaled, $thetas, $phis);

        my $dist_interp;
        # Which method should we use to calculate the flow?
        if ($flow_method eq "wsa"){
            ($r_vr_scaled, $r_fr_scaled, $thetas, $phis) = gen_fluxon_wsaflow($fluxon, $distance_array_degrees, $fluxon_id);
        } elsif ($flow_method eq "parker"){
            ($r_vr_scaled, $r_fr_scaled, $thetas, $phis) = gen_fluxon_tflow($fluxon);
        } elsif ($flow_method eq "schonfeld"){
            ($r_vr_scaled, $r_fr_scaled, $thetas, $phis) = gen_fluxon_schonflow($fluxon);
        } elsif ($flow_method eq "psw"){
            ($r_vr_scaled, $r_fr_scaled, $thetas, $phis) = gen_fluxon_pswflow($fluxon);
        } elsif ($flow_method eq "tempest"){
            ($r_vr_scaled, $r_fr_scaled, $thetas, $phis) = gen_fluxon_tempestflow($fluxon);
        } elsif ($flow_method eq "cranmer"){
            ($r_vr_scaled, $r_fr_scaled, $thetas, $phis) = gen_fluxon_cranmerflow($fluxon);
        } else {
            die "Invalid flow method: $flow_method";
        }


        my $result = {
            fluxon_position  => pdl($fluxon_id),

            phi_base    => squeeze(pdl($phis(0)    )),
            theta_base  => squeeze(pdl($thetas(0)  )),
            phi_end     => squeeze(pdl($phis(-1)   )),
            theta_end   => squeeze(pdl($thetas(-1) )),

            radial_velocity_base    => squeeze($r_vr_scaled(1, 0   )),
            radial_velocity_end     => squeeze($r_vr_scaled(1, -1  )),

            # flux_expansion_base     => squeeze($r_fr_scaled(1, 1   )),
            # flux_expansion_end      => squeeze($r_fr_scaled(-2, 1  )),
            # flux_expansion_middle   => squeeze($r_fr_scaled(-1, 1  )),
            flux_expansion_base     => squeeze($r_fr_scaled(1, 0   )),
            flux_expansion_end      => squeeze($r_fr_scaled(1, -1  )),
            # flux_expansion_middle   => squeeze($r_fr_scaled(1, 1  )),
        };


        $fork_manager->finish(0, $result);
    }

    $fork_manager->wait_all_children;

    my $after_time = clock_gettime();
    my $elapsed_time = $after_time - $before_time;
    my $rounded_elapsed_time = rint($elapsed_time * 10) / 10;

    print "\n\n\t\tWind Calculation Complete\n\n";

    # if ($flow_method eq "wsa"){
    #     do_image_plot();
    # }

    @results = sort { $a->{fluxon_position} <=> $b->{fluxon_position} } @results;

        wcols
        # print
        pdl(map { $_->{fluxon_position} } @results),


        pdl(map { $_->{phi_base}  } @results),
        pdl(map { $_->{theta_base} } @results),
        pdl(map { $_->{phi_end}    } @results),
        pdl(map { $_->{theta_end}  } @results),


        squeeze(pdl(map { $_->{radial_velocity_base} } @results)),
        squeeze(pdl(map { $_->{radial_velocity_end} } @results)),

        squeeze(pdl(map { $_->{flux_expansion_base} } @results)),
        squeeze(pdl(map { $_->{flux_expansion_end} } @results)),
        # squeeze(pdl(map { $_->{flux_expansion_middle} } @results)),
        $output_file_name;

    # Delete the temporary directory.
    rmtree($temp_dir);
}

if ($0 eq __FILE__) {
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
    use pipe_helper qw(configurations);

    my %configs = configurations();
    my $datdir = $configs{datdir};

    my $cr = 2160;
    my $batch_name = "fluxon_paperfigs";

    # Pathing
    my $world_out_dir = $datdir."$batch_name/cr".$cr.'/rlx/';
    my $full_world_path = $world_out_dir . "cr2160_relaxed_s4000.flux";
    my $wind_out_dir = $datdir."$batch_name/cr".$cr.'/wind';
    my $wind_out_file = "$wind_out_dir/radial_wind.dat";

    # Loading the world
    my $this_world_relaxed = read_world($full_world_path);
    $this_world_relaxed->update_force(0);
    my @fluxons = $this_world_relaxed->fluxons;

    map_fluxon_flow_parallel_master($wind_out_file, \@fluxons);
}

1;
__END__
