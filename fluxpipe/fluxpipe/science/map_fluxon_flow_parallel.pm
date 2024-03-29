package map_fluxon_flow_parallel;

use strict;
use warnings;

use Exporter qw(import);
our @EXPORT_OK = qw(map_fluxon_flow_parallel);

use PDL::NiceSlice;
use PDL::Options;
use Parallel::ForkManager;

=head2 map_fluxon_flow

=for ref

Given a world and list of fluxons, generates a mapping of the solar wind flow along these fluxon structures.

=cut

sub map_fluxon_flow_parallel {
    my $fn = shift;
    my $fluxons = shift;
    my @fluxons = @{$fluxons};

    my $max_processes = 4; # Set the number of parallel processes you want to run
    my $fork_manager = new Parallel::ForkManager($max_processes);

# Generate some blank storage arrays
my @fps = ();   # Fluxon array position
my @phb = ();   # Begining coordinates
my @thb = ();
my @phe = ();   # End coordinates
my @the = ();
my @vrb = ();   # Beginning and end wind velocity
my @vre = ();
my @frb = ();   # Beginning and end expansion factors
my @fre = ();
my @fre2 = ();

# Loop through open fluxons and generate wind profiles
for my $fid(0..scalar(@fluxons)-1){
    $fork_manager->start and next; # Start a new process and proceed to the next iteration

    if ($fid>0){
    print "\r", '  Calculating fluxon:', $fid, '/', scalar(@fluxons);
    }
    my $me = $fluxons[$fid];

    # Check for open fieldlines
    my $st_open = ($me->{fc_start}->{label}==-1);
    my $en_open = ($me->{fc_end}->{label}==-2);
    if ($me->{plasmoid} || ($st_open + $en_open != 1)){
        $fork_manager->finish; # Finish the process and move to the next iteration
        next;
    }

    # Find the transonic wind solution
    (my $farr, my $fr, my $bth, my $bph) = gen_fluxon_tflow($me);

    # Append any needed values and move along home
    # Note to distinguish between fluxons ends...
    push(@fps, $fid);
    push(@phb, squeeze($bph(0,0)));
    push(@thb, squeeze($bth(0,0)));
    push(@phe, squeeze($bph(0,1)));
    push(@the, squeeze($bth(0,1)));
    push(@vrb, $farr(1,0));
    push(@vre, $farr(1,-1));
    push(@frb, $fr(1,1));
    push(@fre, $fr(-2,1));
    push(@fre2, $fr(-1,1));

    # Output this data to disk
    wcols pdl(@fps), pdl(@phb), pdl(@thb), pdl(@phe), pdl(@the), squeeze(pdl(@vrb)), squeeze(pdl(@vre)), squeeze(pdl(@frb)), squeeze(pdl(@fre)), squeeze(pdl(@fre2)), $fn;

    $fork_manager->finish; # Finish the process
}

$fork_manager->wait_all_children; # Wait for all processes to finish

# Output this data to disk
wcols pdl(@fps), pdl(@phb), pdl(@thb), pdl(@phe), pdl(@the), squeeze(pdl(@vrb)), squeeze(pdl(@vre)), squeeze(pdl(@frb)), squeeze(pdl(@fre)), squeeze(pdl(@fre2)), $fn;

}

__END__
