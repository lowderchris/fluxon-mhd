use strict;
use warnings;
use PDL;

use Test::More tests=>10;

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
