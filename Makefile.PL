use ExtUtils::MakeMaker;
use Config;

WriteMakefile(
               'NAME'     => 'Math::Integral::Romberg',
               'LINKTYPE' => '$(INST_PM)',
               'VERSION_FROM'  => 'Romberg.pm',
               'dist'     => { SUFFIX => 'gz', COMPRESS => 'gzip -f' },
	      'clean'   => {'FILES' => '*~'},
);
