use ExtUtils::MakeMaker;
use Config;

# Find the pdlcore.h and pdl.h include files 
@inc = ();
{

  local $_;
  incfoo:foreach(@INC) {
    $dir = "$_/PDL/Core";
    if(-e "$dir/pdl.h"){  # pdlcore.h should be in the same place as pdl.h
	push(@inc,$dir);
	last incfoo;
    }
  }
  print STDERR "##########\nWarning: Flux::Core could not find pdl.h\n#########\n" unless @inc;
}

WriteMakefile( NAME => 'Flux::Core',
	       INC=> $inc,
	       CCFLAGS=> $ccflags,
	       OBJECT => '$(BASEEXT)$(OBJ_EXT) /usr/local/lib/libflux.a'
	       );

