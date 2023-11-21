use ExtUtils::MakeMaker;
use File::ShareDir::Install;


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

install_share dist => 'PDL'; #recover this with dist_dir from File::ShareDir.

WriteMakefile( DIR=>[ 'Flux' ],
	       NAME=> 'Flux',
	       INC=>$inc,
	       LIBS=>$libs,
	       CONFIGURE_REQUIRES=> {"ExtUtils::MakeMaker" => 6.52,
					  "File::ShareDir" => 0,
					  "File::ShareDir::Install" => 0,
	       },
	       PREREQ_PM=> { "PDL" => 2.007_13,
	       },
	       );

package MY;
use File::ShareDir::Install qw( postamble );
use File::ShareDir;

sub postamble {
    my $self = shift;
    my @ret = File::ShareDir::Install::postamble( $self );
#    print "\n\nThere is/are " . (0+@ret) . " lines in \@ret, and the contents are:\n@ret\n";
#    my $dir = File::ShareDir::dist_dir('Flux');
    push @ret, <<'EOM';
install::
	$(NOECHO) $(ABSPERLRUN) -MFlux -MFile::ShareDir -e '$$dir=File::ShareDir::dist_dir(q/Flux/);  print qq{\n\nAutoloading routines for FLUX have been installed. To achieve autoloading functionality, please add the following directory (with the plus sign) to the colon-separated list of directories in your PDLLIB environment variable: \n\n\t+$$dir \n\n\
Also make sure the following four lines are in your ~/.perldlrc file: \n\
\n\
require(q|PDL/default.perldlrc|); \n\
use PDL::AutoLoader; \n\
\$$PDL::AutoLoader::Rescan=1; \n\
1;\n\
\
};'\
EOM

    return join "\n", @ret;
}
