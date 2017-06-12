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
eval {$so->render;};
is($@,'','rendering does not throw an error');
