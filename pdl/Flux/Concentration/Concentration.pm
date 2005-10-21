=head1 NAME

Flux::Concentration - Fluxon MHD boundary condition (line tying post)

=head1 SYNOPSIS

 use PDL;
 use Flux;
  
 $world = read_world($filename);
 $fc = $world->concentration(15);
 print %$fc;

=head1 DESCRIPTION

Flux::Concentration objects are hashes tied to the FLUX_CONCENTRATION structures in the
FLUX library.  They represent endpoints of line-tied field lines.  

There are several special flux concentations in every Flux::World -
the fc_ob, fc_oe, fc_pb, and fc_pe concentrations are image
concentrations that are used for open and plasmoid boundary
conditions.  They act like normal flux concentrations for
accounting-of-flux purposes, but the location data are invalid.

VERSION 

This is version 1.1 of Concentration.pm.

=head1 FUNCTIONS

=cut





package Flux::Concentration;

require Exporter;
require DynaLoader;
our @ISA = qw(Exporter DynaLoader Flux);
our @EXPORT = ();

############################## 
# No bootstrap required just yet -- no .xs yet.
#  bootstrap Flux::Concentration;  




=pod

=head2 new_from_ptr

=for usage

$fc = new_from_ptr Flux::Concentration($ptr);

=for ref

Constructor of the perl object -- you must feed in a long-int pointer to the underlying
structure in the C arena.

=cut

sub new_from_ptr {
    my $class = shift;
    my $ptr = shift;
    my %hash;
    tie %hash,"Flux::Concentration",$ptr;
    return bless(\%hash,$class);
}


######################################################################
# TIED INTERFACE
# Mostly relies on the general utility functions in Flux....

sub TIEHASH {
    my $class = shift;
    my $ptr = shift;
    my $me = \$ptr;
    return bless($me,$class);
}

sub FETCH {
    my($me, $field)=@_;
    my $code = $Flux::codes->{concentration}->{$field};
 
    return undef unless defined($code);
    
    Flux::r_val( $me, $Flux::typecodes->{concentration}, @$code[0..1] );
}

sub STORE {
    my($me, $field,$val) = @_;
    my $code = $Flux::codes->{concentration}->{$field};
    return undef unless defined($code);
    Flux::w_val( $me, $Flux::typecodes->{concentration}, @$code[0..1], $val );
}

sub DELETE {
    print STDERR "Warning: can't delete fields from a tied CONCENTRATION hash\n";
    return undef;
}

sub CLEAR {
    print STDERR "Warning: can't clear a tied CONCENTRATION hash\n";
    return undef;
}

sub FIRSTKEY {
    return "world";
}

sub NEXTKEY {
    my ($class,$prev) = @_;
    return $Flux::ordering->{concentration}->{$prev};
    
}

sub SCALAR {
    _stringify(@_);
}



1;
