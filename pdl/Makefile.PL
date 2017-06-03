use ExtUtils::MakeMaker;


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
  print STDERR "###########\nWARNING: Flux could not find pdl.h\n#########\n" unless @inc;
}
push(@inc, ($ENV{FL_PREFIX} // "/usr/local") . "/include");
$inc = join(" ",map { "-I$_" } @inc);
$libs = '-L' . ($ENV{FL_PREFIX} // '/usr/local') . '/lib -lflux';

#### Check for clang and make ccflags accordingly.
use Config;
$using_clang = ($Config{cc} =~ /clang/);
unless($using_clang){
    $line = `cc -v 2>&1`;
    $using_clang = ($line =~ /clang/io);
}
#$ccflags = $using_clang ? '' : '-no-inline';

WriteMakefile( DIR=>[ 'Flux' ],
	       NAME=> 'Flux',
#	       OBJECT => '$(BASEEXT)$(OBJ_EXT) /usr/local/lib/libflux.a'
	       );
