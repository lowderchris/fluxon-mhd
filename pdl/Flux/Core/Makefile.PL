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
}

WriteMakefile( NAME => 'Flux::Core',
	       INC=>join(" ",map { "-I$_"} @inc),
	       CCFLAGS=> $Config{cc} =~ /clang/ ? '' : '-no-inline',
	       OBJECT => '$(BASEEXT)$(OBJ_EXT) /usr/local/lib/libflux.a'
	       );

