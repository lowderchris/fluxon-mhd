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
# use PDL::Graphics::Gnuplot qw(gnuplot_closeall);
# gnuplot_closeall();


# Control Flags
$PDL::verbose = 0;
# no warnings 'cleanup';
# my $do_max = 10000;

=head1 DESCRIPTION

This Perl module provides functionality to map the solar wind flow along fluxon structures in a parallelized manner. It uses the Parallel::ForkManager library to manage concurrent processes and PDL (Perl Data Language) for numerical computations.

=head1 AUTHOR

Gilly <gilly@swri.org> and others!

=head1 LICENSE

This library is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut


=head2 map_fluxon_flow_parallel_master

=over

=item B<Description>

Given a world and list of fluxons, this function generates a mapping of the solar wind flow along these fluxon structures. It uses parallel processing to speed up the calculations.

=item B<Parameters>

=over

=item C<$file_name>: The name of the output file where the results will be saved.

=item C<$fluxons>: A reference to an array containing the fluxons to be mapped.

=item C<$max_processes>: (Optional) The maximum number of parallel processes to use. Defaults to the value of C<$concurrency>.

=back

=item B<Returns>

The function writes the results to disk and returns nothing.

=back

=cut




my $concurrency = 11;
my $temp_dir = "fluxons/fluxon-mhd/temp";

#TODO Make the gen_fluxon_flow and gen_fluxon_tflow functions get loaded only once.

=head2 map_fluxon_flow_parallel_master

=for ref

Given a world and list of fluxons, generates a mapping of the solar wind flow along these fluxon structures.

=cut

sub map_fluxon_flow_parallel_master {
    my $file_name      = shift;
    my $fluxons = shift;
    my @fluxons = @{$fluxons};
    make_path($temp_dir) unless -d $temp_dir;
    my $max_processes = shift || $concurrency;    # Set the number of parallel processes
    my $highest_fid_done = 0;
    my $before = clock_gettime();


    my $n_choke = 12;
    my $choke = 0;
    my $do_plot_charts = 0;
    # my $chartpath = "/Users/cgilbert/vscode/fluxons/Fluxon-Scripts-Gilly/fluxons/profiles/";

    # Count the number of open and closed fluxons
    my $n_closed = 0;
    my @open_inds = ();
    my $n_open = 0;
    for my $fid ( 0 .. scalar(@fluxons) - 1 ) {
        ## Condition the Input ####
        my $me = $fluxons[$fid];

        # Check for open fieldline, and skip if not
        my $st_open = ( $me->{fc_start}->{label} == -1 );
        my $en_open = ( $me->{fc_end}->{label} == -2 );
        if ( $me->{plasmoid} || ( $st_open + $en_open != 1 ) ) {
            $n_closed = $n_closed + 1;
        } else {
            $n_open = $n_open + 1;
            push(@open_inds, $fid);
        }
    }

    my $itercount = -1;
    my $fork_manager = Parallel::ForkManager->new( $max_processes, $temp_dir );

    # Generate a blank storage array
    my @results = ();

    # Set the finish callback to update $highest_fid_done and print it
    $fork_manager->run_on_finish(sub {
        $itercount = $itercount + 1;
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $result) = @_;
        if ($exit_code == 0) {
            push(@results, $result);
            my $this_fid = $result->{fps}+1;
            # print "\nExit_code: $exit_code, Result: $this_fid, Highest: $highest_fid_done, Itercount: $itercount\n";

            if ($this_fid > $highest_fid_done) {
                $highest_fid_done = $this_fid;

                # Print timing information
                my $after   = clock_gettime();
                my $elapsed = $after - $before;
                my $round_elapsed = rint($elapsed*10) / 10;
                my $remaining = $n_open - $highest_fid_done;
                my $time_each = $elapsed / $highest_fid_done;
                my $time_remaining = $remaining * $time_each;
                my $round_remain = rint($time_remaining*10) / 10;
                my $round_minutes = rint($round_remain / 6)/10;

                print "\r", "\t\t\tCalculated fluxons: ", $highest_fid_done, ' of ', $n_open, ", ", $round_elapsed, "(s) elapsed, ", $round_remain, "(s) [$round_minutes mins] remaining......";
            }
        }
    });

    #     # create chart object
    # my $chart = Chart::Gnuplot->new(
    #     output => "/Users/cgilbert/vscode/Fluxon-Scripts-Gilly/this_line_all.png",
    #     title => "My Plot",
    #     xlabel => "radius",
    #     ylabel => "wind",
    # );




    my $n_tot = scalar(@fluxons);
    my $n_act = 2*($n_tot - $n_open) + $n_open;
    print "\n\t\tRunning with $max_processes Cores on $n_open open Fluxons out of $n_tot total Fluxons from $n_act footpoints!\n\n\t\tBeginning Calculation...\n\n";



    $PDL::verbose = 0;
    no warnings;

    my $max_fid = scalar(@open_inds) - 1;
    if ($choke){
        if ($n_choke < $max_fid){
            $max_fid = $n_choke;
        }
    }
    for my $fid ( 0 .. $max_fid - 1 ) {
        ## Inside the Loop ##################################################
        # Start a new process and proceed to the next iteration
        $fork_manager->start and next;

        ## Condition the Input ####
        my $me = $fluxons[$fid];

        # ## Run the Main Algorithm ####
        # # Find the transonic wind solution
        (my $farr, my $fr, my $bth, my $bph) = gen_fluxon_tflow($me);


        # if ($do_plot_charts) {
        #     my @rr =  $farr->( 0 , : )->list;
        #     my @ww =  $farr->( 1 , : )->list;
        #     # if $chartpath does not exist then make it
        #     if ( ! -d $chartpath ) {
        #         make_path($chartpath);
        #     }

        #         # create chart object
        #     my $chart = Chart::Gnuplot->new(
        #         output => $chartpath . "this_line_op$n_open\_$fid.png",
        #         title => "My Plot",
        #         xlabel => "radius",
        #         ylabel => "wind",
        #         bg => "white",
        #         logscale => "xy",  # set both x and y axes to logarithmic scale
        #     );

        #     # define and plot data set
        #     my $dataSet = Chart::Gnuplot::DataSet->new(
        #         xdata => \@rr,
        #         ydata => \@ww,
        #         style => "lines",
        #     );

        #     # # create a data set for the horizontal line
        #     # my $horizontalLine = Chart::Gnuplot::DataSet->new(
        #     #     xdata => [$rr[0], $rr[-1]],  # x-values from the start to the end of your data
        #     #     ydata => [$sound_speed, $sound_speed],  # constant y-value
        #     #     style => "lines",
        #     #     linetype => "dashed",  # set line type to dashed
        #     # );

        #     # plot the horizontal line
        #     $chart->plot2d($dataSet);
        #     # $chart->plot2d($dataSet, $horizontalLine);
        # }

        ## Report the Results ####
        # Create a result hash
        my $result = {
            fps  => $fid,                       # Fluxon ID
            phb  => squeeze( $bph ( 0, 0 ) ),   # Phi at the base
            thb  => squeeze( $bth ( 0, 0 ) ),   # Theta at the base
            phe  => squeeze( $bph ( 0, 1 ) ),   # Phi at the end
            the  => squeeze( $bth ( 0, 1 ) ),   # Theta at the end
            vrb  => $farr ( 1, 0 ),             # Radial velocity at the base
            vre  => $farr ( 1, -1 ),            # Radial velocity at the end
            frb  => $fr   ( 1,  1 ),            # flux expansion at the base
            fre  => $fr   ( -2, 1 ),            # flux expansion at the end
            fre2 => $fr   ( -1, 1 ),            # flux expansion at the middle
        };

        # Finish the process and return the result
        $fork_manager->finish( 0, $result );
    }


    ## Outside of the Loop ## #############################################
    $fork_manager->wait_all_children;    # Wait for all processes to finish

    # $chart->plot2d($dataSet);

    my $after   = clock_gettime();
    my $elapsed = $after - $before;
    my $round_elapsed = rint($elapsed*10) / 10;
    # print "\n\n\t\tWind calculation took: $round_elapsed(s) with $max_processes cores\n \n";
    print "\n\n\t\tWind Calculation Complete\n\n";
    ## Output the Results ##
    # Sort the results by fps
    @results = sort { $a->{fps} <=> $b->{fps} } @results;

    # Write the results to disk
    wcols pdl( map { $_->{fps} } @results ),
        pdl( map { $_->{phb} } @results ),
        pdl( map { $_->{thb} } @results ),
        pdl( map { $_->{phe} } @results ),
        pdl( map { $_->{the} } @results ),
        squeeze( pdl( map { $_->{vrb} } @results ) ),
        squeeze( pdl( map { $_->{vre} } @results ) ),
        squeeze( pdl( map { $_->{frb} } @results ) ),
        squeeze( pdl( map { $_->{fre} } @results ) ),
        squeeze( pdl( map { $_->{fre2} } @results ) ),
        $file_name;

    # delete the temp_dir
    rmtree($temp_dir);
}


if ($0 eq __FILE__) {
    # use warnings;
    use PDL::AutoLoader;
    use PDL;
    use PDL::Transform;
    use PDL::NiceSlice;
    use PDL::Options;
    # # use PDL::ImageND;
    use Flux;
    use PDL::IO::Misc;
    use File::Path;
    use Time::HiRes qw(clock_gettime);
    use File::Basename qw(fileparse);
    use pipe_helper qw(configurations);

    my %configs = configurations();
    my $datdir = $configs{datdir};

    # my $datdir = "/Users/cgilbert/vscode/fluxons/Fluxon-Scripts-Gilly/";
    # push(@PDLLIB,"+".$datdir);
    # push(@PDLLIB,"+/Users/cgilbert/vscode/fluxons/Fluxon-Scripts-Gilly");
    # push(@PDLLIB,"+/Users/cgilbert/vscode/fluxons/fluxon-mhd/pdl/PDL");
    # push(@INC, "+/Users/cgilbert/opt");
    # push(@INC, "/Users/cgilbert/.cpan/build");




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

            # Pass the functions to the child processes using run_on_start callback
            # $fork_manager->run_on_start(sub {
            #     my ($pid, $ident) = @_;
            #     # Access the function references from the global variable
            #     my $gen_fluxon_tflow = $FUNCTION_REFERENCES->{gen_fluxon_tflow};
            #     my $gen_fluxon_flow = $FUNCTION_REFERENCES->{gen_fluxon_flow};
            #     # Now $gen_fluxon_tflow and $gen_fluxon_flow contain the functions in the child process
            # });

                    # if ( $choke && ($fid > $n_choke) ) {
        #     # For only doing a subset of the inputs
        #     # Finish the process and move to the next iteration
        #     $fork_manager->finish(1);
        #     next;
        # }
        # # Check for open fieldline, and skip if not
        # my $st_open = ( $me->{fc_start}->{label} == -1 );
        # my $en_open = ( $me->{fc_end}->{label} == -2 );
        # if ( $me->{plasmoid} || ( $st_open + $en_open != 1 ) ) {
        #     $fork_manager->finish(1);    # Finish the process and move to the next iteration
        #     # finish( 0, $result )
        #     next;
        # }

                # # Print timing information
        # my $after   = clock_gettime();
        # my $elapsed = $after - $before;
        # my $round_elapsed = rint($elapsed*10) / 10;
        # my $remaining = $n_open - $fid;
        # my $time_each = $elapsed / $fid;
        # my $time_remaining = $remaining * $time_each;
        # print "\r", '  Calculated fluxons: ', $fid, ' of ', $n_open, ", ", $round_elapsed, "(s) elapsed, ", $time_remaining, "(s) remaining.";


        # if (grep { $_ == $fid } @open_inds) {}
        # else {
        #     # Finish the process and move to the next iteration
        #     print "FAILLLL";
        #     $fork_manager->finish(1);
        #     next;
        # }