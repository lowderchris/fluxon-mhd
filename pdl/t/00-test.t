use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Constants qw/PI/;
use File::Temp qw/tempfile/;

use Test::More tests=>27;

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
my @pref_terms = ('pngcairo','pdfcairo','wxt','x11','postscript','svg');
my $dev;
foreach (@pref_terms){
    if (grep $_,@gpterms){
	$dev=$_;
	last;
    }
}

 SKIP: {
     skip "No suitable gnuplot terminal found for rendering",1 unless $dev;
     eval {$so->render({dev=>$dev});};
     is($@,'','rendering does not throw an error') or diag("You might need to change the gnuplot terminal type. Availble terminals are:\n" . join(' ',@Alien::Gnuplot::terms));
};

my ($zn,$xn,$yn) = (11,3,5);
my $z = zeroes(21)->xlinvals(0,1);
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

#some extra options to render, mostly just to make sure it doesn't crash.
#NEED TO CHECK WHETHER THE NEIGHBOR AND HULL PLOTS ACTUALLY MAKE ANY SENSE
eval{$world->render({neighbors=>1});};
is($@,'','render with neighbors turned on');
eval{$world->render({hull=>1});};
is($@,'','render with hull turned on');

