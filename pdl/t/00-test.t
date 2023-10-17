use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Constants qw/PI/;
use File::Temp qw/tempfile/;

use Test::More tests=>37;

BEGIN {use_ok('Flux');}

my $so = read_world('menagerie/single-open-line.flux');

#check that the world and its components are the correct classes
isa_ok($so, 'Flux::World');
isa_ok($so->{concentrations}, 'Flux::Concentration');
isa_ok($so->{lines}, 'Flux::Fluxon');
isa_ok($so->{vertices}, 'Flux::Vertex');

ok(all($so->{fc_ob}->{x}==zeroes(3)),"Open beginning location");
ok(all($so->{fc_oe}->{x}==zeroes(3)),"Open ending location");
is($so->{fc_ob}->{locale_radius},6,"Open boundary size");
is($so->{auto_open},0,"Open boundary default no auto-open");

#should not barf on rendering a simple world
#First get a list of available terminals:
my @gpterms = @Alien::Gnuplot::terms;
my @pref_terms = ('pngcairo','pdfcairo','qt','wxt','x11','postscript','svg');

my $dev;
foreach (@pref_terms){
    if (grep $_,@gpterms){
	$dev=$_;
	last;
    }
}

 SKIP: {
     skip "No suitable gnuplot terminal found for rendering",1 unless $dev;
     eval {$so->render({dev=>$dev,range=>[-1,1,-1,6,0,1]});};
     is($@,'','rendering does not throw an error') or diag("You might need to change the gnuplot terminal type. Availble terminals are:\n" . join(' ',@Alien::Gnuplot::terms));
};

my ($zn,$xn,$yn) = (11,3,5);
my $z = zeroes($zn)->xlinvals(0,1);
my $x = zeroes($xn)->xlinvals(-1,1);
my $y = zeroes($yn)->xlinvals(-1,1);
my @lines = ();
#this is a stupid way to do it but my index tricks are rusty.
foreach my $yi($y->list){
    foreach my $xi($x->list){
	my $coords = zeroes(3,$z->nelem);
	$coords((0)).=$xi;
	$coords((1)).=$yi;
	$coords((2)).=$z;
	push @lines,$coords;
    }
}
require_ok('PDL/make_world.pdl');
my $string = make_world(@lines,{photosphere=>[1,0,0,0,0,0,1],
				photosphere2=>[1,0,0,1,0,0,-1]});
ok(length($string)>0,'make_world returns a string of some length');
my $world = str2world($string);
isa_ok($world,'Flux::World');
my (undef,$tmpfile)=tempfile();
eval {write_world($world,$tmpfile);} ;
is($@,'','write_world does not throw an error');
ok(length($world->string)>0,'world string');
ok(length("$world")>0,'world stringification');
ok(length($world->summary)>0,'world summary');
(undef,$tmpfile)=tempfile();
$tmpfile .= '.gz';
eval {write_world($world,$tmpfile);} ;
is($@,'','gzip write_world does not throw an error');
eval {read_world($tmpfile);};
is($@,'','gzip read_world does not throw an error');

my @concentrations = $world->concentrations;
my @fluxons = $world->fluxons;
ok(scalar(@fluxons)>0,'list of all fluxons');
my @vertex_ids = $world->vertex_ids;
ok(scalar(@vertex_ids)>0,'list of all vertex ids');
my @vertices = $world->vertices(@vertex_ids[0..9]);
is(scalar(@vertices),10,'vertices');

#insert new flux concentrations with a fluxon between them
my $source_loc = pdl(5,1,0);
my $sink_loc = $source_loc + pdl(1,2,0); #this new concentration is in the same plane as the source, unlike the rest of the concentrations in the $world.
my $src = $world->emerge($source_loc,$sink_loc,1,($sink_loc-$source_loc)/2+$source_loc+pdl(0,0,.9));
$source_loc += pdl(1,0,0);
$sink_loc = $source_loc + pdl(1,2,0);
$src = $world->emerge($source_loc,$sink_loc,1,3);

$world->render;
my @conc2 = $world->concentrations;
is(scalar(@conc2),scalar(@concentrations)+4,'emerge added two flux concentrations each');
my @flux2 = $world->fluxons;
is(scalar(@flux2),scalar(@fluxons)+2,'emerge added one fluxon each');

#fix_curvature should eliminate a bunch of vertices since these are all straight fluxons
$world->update_neighbors(1);
ok($world->fix_curvature(PI/4,0)!=0,'fix_curvature removes or adds at least one vertex in a simple world');
$world->update_neighbors(1);

#some extra options to render, mostly just to make sure it doesn't crash.
#NEED TO CHECK WHETHER THE NEIGHBOR AND HULL PLOTS ACTUALLY MAKE ANY SENSE
eval{$world->render({neighbors=>1});};
is($@,'','render with neighbors turned on');
eval{$world->render({hull=>1});};
is($@,'','render with hull turned on');

#check different values for concurrency
require_ok('PDL/simple_relaxer.pdl');

$world->{'concurrency'}=0;
eval {simple_relaxer($world,0,1,{movie_n=>0,disp_n=>0});};
is($@,'','simple_relaxer with concurrency=0');

$world->{concurrency}=1;
eval {simple_relaxer($world,0,2,{movie_n=>0,disp_n=>0});};
is($@,'','simple_relaxer with concurrency=1');

$world->{concurrency}=2;
eval {simple_relaxer($world,0,3,{movie_n=>0,disp_n=>0});};
is($@,'','simple_relaxer with concurrency=2');

###Test interpolation routines
##start with coordinate interpolation, exactly on a vertex
{
    my $in_loc = $world->vertex($vertex_ids[$#vertex_ids/2])->{'x'};
    my $out_loc;
    eval{$out_loc = $world->interpolate_value('x',$in_loc,1,0,1);};
    is($@,'','interpolation of a coordinate value at a vertex location did not crash');
    ok(all(approx($out_loc,$in_loc,2**-52)),'interpolation of a coordinate value at a vertex location gives correct value');

    $in_loc = pdl(0.43,-0.25,0.6);
    $out_loc = $world->interpolate_value('x',$in_loc,1,0,0);
    ok(all(approx($out_loc,$in_loc,2**-52)),'interpolation within the cartesian world works');

    $in_loc = pdl(-1.25,0.3,0.35);
    $out_loc = $world->interpolate_value('x',$in_loc,1,0,1);
    ok(all(!isfinite($out_loc)),'interpolation outside the domain gives NaNs');


my ($xcheck,$ycheck,$zcheck,$coord,$coord3pt,$coord2pt,$out_loc3,$out_loc2);
my ($done3pt,$done2pt)=(0,0);
my $try=0;
do {
    $xcheck = rand(1)*2 -1;
    $ycheck = rand(1) * 2 -1;
    $zcheck = rand(1) ;

    $coord = pdl($xcheck,$ycheck,$zcheck);
    my @vertices = $world->closest_simplex($coord,0);
    $out_loc = $world->interpolate_value('x',$coord,0,0,1);
    if (scalar(@vertices)==3 && all(isfinite($out_loc))) {
        $done3pt=1;
        $coord3pt = $coord;
        $out_loc3 = $out_loc;
    }
    if (scalar(@vertices)==2 && all(isfinite($out_loc))) {
        $done2pt=1;
        $coord2pt = $coord;
        $out_loc2 = $out_loc;
    }

    $try++;
} until (($done3pt && $done2pt) || $try==1E4);


done_testing();


SKIP: {
          skip 'Did not find a location that gave 3 simplex points', 1 unless $done3pt;
          ok(!all($coord3pt==$out_loc3),'trusting interpolation when only 3 simplex values found gives an incorrect interpolant');
      };

SKIP: {
          skip 'Did not find a location that gave 2 simplex points (this is not unexpected)', 1 unless $done2pt;
          ok(!all($coord2pt==$out_loc2),'trusting interpolation when only 2 simplex values found gives an incorrect interpolant');
      };


}


=pod

TODO: {

    todo_skip 'this generates a segfault',0,1;
use PDL::Graphics::PGPLOT::Window;
my $w = pgwin(dev=>"/xs", size=>[5,5]);
my $points = pdl([0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[1,-1]);
Flux::World::_plot_hull($w, Flux::World::_hull_points($points));
}

=cut
