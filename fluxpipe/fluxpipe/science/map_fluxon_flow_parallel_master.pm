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
use gen_fluxon_ghostsflow qw(gen_fluxon_ghostsflow);
use gen_fluxon_tempestflow qw(gen_fluxon_tempestflow);
use make_tempest_file qw(make_tempest_file);
use pipe_helper qw(configurations);
use PDL::Graphics::Gnuplot;


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

# Subroutine to save full velocity profiles
sub save_full_velocity_profiles {
    my ($results, $file_path) = @_;

    # Ensure results are sorted by fluxon_position
    my @sorted_results = sort { $a->{fluxon_position} <=> $b->{fluxon_position} } @$results;

    open(my $fh, '>', $file_path) or die "Could not open file '$file_path' $!";
    print $fh "fluxon_position,radius,velocity\n";

    foreach my $result (@sorted_results) {
        my $fluxon_id = $result->{fluxon_position}->sclr;
        my @r = listy($result->{r});  # Convert PDL object to Perl array
        my @vr = listy($result->{vr});  # Convert PDL object to Perl array

        # Debug: Print to ensure arrays are being correctly accessed
        if (scalar(@r) != scalar(@vr)) {
            warn "Mismatched array lengths for fluxon ID $fluxon_id: r length=" . scalar(@r) . ", vr length=" . scalar(@vr);
            next;
        }


        for my $i (0 .. $#r) {
            print $fh join(",", $fluxon_id, $r[$i], $vr[$i]), "\n";
        }
    }

    close $fh;


}

sub listy {
    my $pdl_obj = shift;
    return $pdl_obj->list;  # Convert PDL object to Perl array
}



sub map_fluxon_flow_parallel_master {
    my $output_file_name = shift;
    my $fluxon_list = shift;
    # my $distance_array_degrees = shift;
    my $flow_method = shift;
    my $CR = shift;
    my $n_want = shift;

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

    for my $fluxon_index (0..scalar(@fluxons)-1) {
        my $me = $fluxons[$fluxon_index];

        # Check for open fieldline, and skip if not
        my $start_open = ($me->{fc_start}->{label} == -1);
        my $end_open = ($me->{fc_end}->{label} == -2);

        # Calculate array of sphiserical coordinate positions and areas along the fluxon
        # Work along the correct direction depending on which end is open
        if($end_open) {
            my $x = squeeze($me->dump_vecs->(0));
            my $y = squeeze($me->dump_vecs->(1));
            my $z = squeeze($me->dump_vecs->(2));
            our $r1 = ($x**2 + $y**2 + $z**2)->sqrt * 696e6,;
            our $r = 0.5 * $r1->range([[0],[1]],[$r1->dim(0)],'e')->sumover;
            our $bmag = pdl(map {$_->{b_mag}} ($me->vertices));

        } else {
            my $x = squeeze($me->dump_vecs->(0,-1:0:-1));
            my $y = squeeze($me->dump_vecs->(1,-1:0:-1));
            my $z = squeeze($me->dump_vecs->(2,-1:0:-1));
            our $r1 = ($x**2 + $y**2 + $z**2)->sqrt * 696e6,;
            our $r = 0.5 * $r1->range([[0],[1]],[$r1->dim(0)],'e')->sumover;
            our $bmag = pdl(map {$_->{b_mag}} ($me->vertices))->(-1:0:-1);
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

                print "\r", "\t\t\tCalculated $flow_method: ", $highest_fluxon_done, ' of ', $num_open_fluxons, ", ", $rounded_elapsed_time, "(s) elapsed, ", $rounded_remaining_time, "(s) [$rounded_remaining_minutes mins] remaining......";
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

    print "\t\tThe flow method is ''$flow_method.''\n\n\n";


    # Precomputations

    if ($flow_method eq "tempest"){
        # # Make the tempest file
        my $fmap_command =
        "$configs{python_dir} fluxon-mhd/fluxpipe/fluxpipe/science/tempest.py --cr $CR --nwant $n_want";
        # print "\n\n$fmap_command\n\n";
        system($fmap_command) == 0 or ( die "Python script returned error $?", exit );
    } elsif ($flow_method eq "cranmer"){
        my $fmap_command =
        "$configs{python_dir} fluxon-mhd/fluxpipe/fluxpipe/science/cranmer_wind.py --cr $CR --nwant $n_want";
        # print "\n\n$fmap_command\n\n";
        system($fmap_command) == 0 or ( die "Python script returned error $?", exit );

    } elsif ($flow_method eq "wsa"){
        # run the python file footpoint_distances.py
        system("$configs{python_dir}  fluxon-mhd/fluxpipe/fluxpipe/science/footpoint_distances_2.py --cr $CR");
        my $distance_file = $configs{data_dir} . "/batches/" . $configs{batch_name} . "/data/cr" . $CR . "/floc/distances.csv";
        open my $fh, '<', $distance_file or die "Could not open '$distance_file': $!";
        # Read the file line by line and split each line
        my @rows;
        while (my $line = <$fh>) {
            chomp $line;  # Remove newline character
            my @values = split /, /, $line;  # Split the line into values
            push @rows, pdl(@values);  # Convert the list of values into a PDL piddle and store it
        }
        close $fh;

        # Combine all rows into a 2D PDL array
        our $distance_array_degrees = cat(@rows);
    }


    # Loop through each fluxon
    for my $fluxon_id (0..$max_fluxon_id - 1) {
    # for my $fluxon_id (12..13 - 1) {
        $fork_manager->start and next;

        my $fluxon = $fluxons[$fluxon_id];
        my ($r_vr_scaled, $r_fr_scaled, $thetas, $phis);

        my $dist_interp;

        # Which method should we use to calculate the flow?
        if ($flow_method eq "wsa"){
            our $distance_array_degrees;
            ($r_vr_scaled, $r_fr_scaled, $thetas, $phis) = gen_fluxon_wsaflow($fluxon, $distance_array_degrees, $fluxon_id);
            our $vr = $r_vr_scaled(:, 1);
            our $fr = $r_fr_scaled(:, 1);
            our $rn = $r_fr_scaled(:, 0);

        } elsif ($flow_method eq "parker"){
            ($r_vr_scaled, $r_fr_scaled, $thetas, $phis) = gen_fluxon_tflow($fluxon);
            our $vr = $r_vr_scaled(1, :)->transpose;
            our $rn = $r_vr_scaled(0, :)->transpose;
            our $fr = $r_fr_scaled(:, 1);
            our $fn = $r_fr_scaled(:, 0);

        } elsif ($flow_method eq "schonfeld"){
            ($r_vr_scaled, $r_fr_scaled, $thetas, $phis) = gen_fluxon_schonflow($fluxon);

        } elsif ($flow_method eq "psw"){
            ($r_vr_scaled, $r_fr_scaled, $thetas, $phis) = gen_fluxon_pswflow($fluxon);

        } elsif ($flow_method eq "tempest"){
            ($r_vr_scaled, $r_fr_scaled, $thetas, $phis) = gen_fluxon_tempestflow($fluxon, $fluxon_id, $output_file_name);
            our $vr = $r_vr_scaled(:, 1);
            our $fr = $r_fr_scaled(:, 1);
            our $rn = $r_fr_scaled(:, 0);

        } elsif ($flow_method eq "cranmer"){
            ($r_vr_scaled, $r_fr_scaled, $thetas, $phis) = gen_fluxon_cranmerflow($fluxon, $fluxon_id);
            our $vr = $r_vr_scaled(:, 1);
            our $fr = $r_fr_scaled(:, 1);
            our $rn = $r_fr_scaled(:, 0);

        } elsif ($flow_method eq "ghosts") {
            ($r_vr_scaled, $r_fr_scaled, $thetas, $phis) = gen_fluxon_ghostsflow($fluxon, $fluxon_id);
            our $vr = $r_vr_scaled(:, 1);
            our $fr = $r_fr_scaled(:, 1);
            our $rn = $r_fr_scaled(:, 0);
        } else {
            die "Invalid flow method: $flow_method";
        }

        our $vr;
        our $fr;
        our $rn;
        my $zn = $rn - 1;
        my $ones = ones($rn->dims);

        my $bot_ind = 1;
        my $top_ind = -2;


        if (0) {
            my $plot = PDL::Graphics::Gnuplot->new(persist => 1);
            $plot->plot({title=>"Fluxon #$fluxon_id", logscale=>'xy', xlabel=>'z = r/R - 1', font=>'48'}, {legend=>"vr"},$zn, $vr);
            $plot->replot({legend=>"vr", with =>'points'},$zn, $vr);


            $plot->replot({legend=>"fr"},$zn, $fr);
            $plot->replot({legend=>"fr", with=>"points"},$zn, $fr);

            $plot->replot({legend=>"theta"},$zn, $thetas);
            $plot->replot({legend=>"theta", with=>"points"},$zn, $thetas);

            $plot->replot({legend=>"phi"},$zn, $phis);
            $plot->replot({legend=>"phi", with=>"points"},$zn, $phis);

            $plot->replot({legend=>"z"},$zn, $zn);
            $plot->replot({legend=>"z", with=>"points"},$zn, $zn);

            my $z0 = $zn(0);

            print "\nZ0 = $z0\n\n";

            <STDIN>;
        }

        my $result = {
            fluxon_position  => pdl($fluxon_id),

            r => $rn,
            vr => $vr,

            phi_base    => squeeze(pdl($phis(0)    )),
            theta_base  => squeeze(pdl($thetas(0)  )),
            phi_end     => squeeze(pdl($phis($top_ind)   )),
            theta_end   => squeeze(pdl($thetas($top_ind) )),

            radial_velocity_base    => squeeze(pdl($vr($bot_ind))),
            radial_velocity_end     => squeeze(pdl($vr($top_ind))),

            flux_expansion_base     => squeeze(pdl($fr($bot_ind))),
            flux_expansion_end      => squeeze(pdl($fr($top_ind))),
        };

        if (0) {
            foreach my $key (keys %$result) {
                my $value = $result->{$key};
                print "$key: $value\n";
            }
        print "\n\n";
        }

        $fork_manager->finish(0, $result);
    }

        $fork_manager->wait_all_children;

    # Extract the directory part of $output_file_name
    my ($output_filename, $output_dir) = fileparse($output_file_name);

    # Create a new directory for full velocity profile results within the output directory
    my $new_results_dir = File::Spec->catdir($output_dir, "full_velocity_profiles");
    make_path($new_results_dir) unless -d $new_results_dir;

    # Modify the filename to indicate it contains full velocity profiles
    my $results_file = File::Spec->catfile($new_results_dir, "results_${flow_method}_full_velocity.dat");
    save_full_velocity_profiles(\@results, $results_file);

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
