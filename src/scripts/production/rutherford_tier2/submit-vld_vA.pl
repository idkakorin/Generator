#!/usr/bin/perl

#----------------------------------------------------------------------
# Submit standard vA event generation jobs for GENIE release validation
#
# Syntax:
#   perl submit-vld_vA.pl <options>
#
# Options:
#  --run           : Comma separated list of run numbers
#  --version       : GENIE version number
# [--production]   :
# [--cycle]        :
# [--use-valgrind] :
#
# Examples:
#  perl submit-vld_vA.pl --production 2.5.1_prelease_tests --cycle 01 --version v2.5.1 --run 1001
#  perl submit-vld_vA.pl --production 2.5.1_prelease_tests --cycle 01 --version v2.5.1 --run 1000,1001,9203
#  perl submit-vld_vA.pl --production 2.5.1_prelease_tests --cycle 01 --version v2.5.1 --run all
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#----------------------------------------------------------------------
#
# SAMPLES:
#......................................................................
#  run   | nev  |  init state      | energy  | processes
#  nu.   |      |                  | (GeV)   | enabled
#......................................................................
#  1000  | 100k |  numu    + n     |  0.5    | all       
#  1001  | 100k |  numu    + n     |    1    | all       
#  1002  | 100k |  numu    + n     |    5    | all       
#  1003  | 100k |  numu    + n     |   50    | all       
#  1101  | 100k |  numubar + p     |    1    | all       
#  1102  | 100k |  numubar + p     |    5    | all       
#  1103  | 100k |  numubar + p     |   50    | all       
#  2001  | 100k |  numu    + Fe56  |    1    | all       
#  2002  | 100k |  numu    + Fe56  |    5    | all       
#  2003  | 100k |  numu    + Fe56  |   50    | all       
#  2101  | 100k |  numubar + Fe56  |    1    | all      
#  2102  | 100k |  numubar + Fe56  |    5    | all       
#  2103  | 100k |  numubar + Fe56  |   50    | all       
#  9001  | 100k |  numu    + Fe56  |    5    | DIS charm 
#  9002  | 100k |  numu    + Fe56  |    5    | QEL charm 
#  9101  | 100k |  numu    + Fe56  |    2    | COH CC+NC 
#  9201  | 100k |  nue     + Fe56  |    1    | ve elastic
#  9202  | 100k |  numu    + Fe56  |    1    | ve elastic
#  9203  | 100k |  numu    + Fe56  |   20    | IMD
#......................................................................
#

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--run')           { $runnu         = $ARGV[$iarg+1]; }
  if($_ eq '--version')       { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')  { $use_valgrind  = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined benchmark runs #. Use the --run option]")
unless defined $runnu;
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$use_valgrind = 0                            unless defined $use_valgrind;
$production   = "benchmarks\_$genie_version" unless defined $production;
$cycle        = "01"                         unless defined $cycle;

$queue          = "prod";
$genie_inst_dir = "/opt/ppd/t2k/GENIE/";
$genie_setup    = "$genie_inst_dir/$genie_version-setup";
$jobs_dir       = "/opt/ppd/t2k/GENIE/scratch/vA-$production\_$cycle";
$xspl_file      = "/opt/ppd/t2k/GENIE/data/job_inputs/xspl/gxspl-t2k-$genie_version.xml";
$mcseed         = 210921029;

%nevents_hash = ( 
  '1000' => '100000',
  '1001' => '100000',
  '1002' => '100000',
  '1003' => '100000',
  '1101' => '100000',
  '1102' => '100000',
  '1103' => '100000',
  '2001' => '100000',
  '2002' => '100000',
  '2003' => '100000',
  '2101' => '100000',
  '2102' => '100000',
  '2103' => '100000',
  '9001' => '100000',
  '9002' => '100000',
  '9101' => '100000',
  '9201' =>  '50000',
  '9202' =>  '50000',
  '9203' =>  '50000'  
);

%nupdg_hash = ( 
  '1000' =>  '14',
  '1001' =>  '14',
  '1002' =>  '14',
  '1003' =>  '14',
  '1101' => '-14',
  '1102' => '-14',
  '1103' => '-14',
  '2001' =>  '14',
  '2002' =>  '14',
  '2003' =>  '14',
  '2101' => '-14',
  '2102' => '-14',
  '2103' => '-14',
  '9001' =>  '14',
  '9002' =>  '14',
  '9101' =>  '14',
  '9201' =>  '12',
  '9202' =>  '14',
  '9203' =>  '14'  
);

%tgtpdg_hash = ( 
  '1000' => '1000000010',
  '1001' => '1000000010',
  '1002' => '1000000010',
  '1003' => '1000000010',
  '1101' => '1000010010',
  '1102' => '1000010010',
  '1103' => '1000010010',
  '2001' => '1000260560',
  '2002' => '1000260560',
  '2003' => '1000260560',
  '2101' => '1000260560',
  '2102' => '1000260560',
  '2103' => '1000260560',
  '9001' => '1000260560',
  '9002' => '1000260560',
  '9101' => '1000260560',
  '9201' => '1000260560',
  '9202' => '1000260560',
  '9203' => '1000260560' 
);

%energy_hash = ( 
  '1000' =>  '0.5',
  '1001' =>  '1.0',
  '1002' =>  '5.0',
  '1003' => '50.0',
  '1101' =>  '1.0',
  '1102' =>  '5.0',
  '1103' => '50.0',
  '2001' =>  '1.0',
  '2002' =>  '5.0',
  '2003' => '50.0',
  '2101' =>  '1.0',
  '2102' =>  '5.0',
  '2103' => '50.0',
  '9001' =>  '5.0',
  '9002' =>  '5.0',
  '9101' =>  '2.0',
  '9201' =>  '1.0',
  '9202' =>  '1.0',
  '9203' => '20.0' 
);

%gevgl_hash = ( 
  '1000' => 'Default',
  '1001' => 'Default',
  '1002' => 'Default',
  '1003' => 'Default',
  '1101' => 'Default',
  '1102' => 'Default',
  '1103' => 'Default',
  '2001' => 'Default',
  '2002' => 'Default',
  '2003' => 'Default',
  '2101' => 'Default',
  '2102' => 'Default',
  '2103' => 'Default',
  '9001' => 'DIS-CHARM',
  '9002' => 'QEL-CHARM',
  '9101' => 'COH',
  '9201' => 'NUE-EL',
  '9202' => 'NUE-EL',
  '9203' => 'IMD'  
);

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

print "Input runs: $runnu \n";

for my $curr_runnu (keys %gevgl_hash)  {
  print "Checking benchmark run: ...... $curr_runnu \n";

  if($runnu=~m/$curr_runnu/ || $runnu eq "all") {
    print "** matched -> submitting job \n";

    #
    # get runnu-dependent info
    #
    $nev   = $nevents_hash {$curr_runnu};
    $nu    = $nupdg_hash   {$curr_runnu};
    $tgt   = $tgtpdg_hash  {$curr_runnu};
    $en    = $energy_hash  {$curr_runnu};
    $gevgl = $gevgl_hash   {$curr_runnu};

    $batch_script  = "$jobs_dir/job_vA-$curr_runnu.pbs";
    $logfile_evgen = "$jobs_dir/job_vA-$curr_runnu.evgen.log";
    $logfile_conv  = "$jobs_dir/job_vA-$curr_runnu.conv.log";
    $logfile_pbse  = "$jobs_dir/job_vA-$curr_runnu.pbs_e.log";
    $logfile_pbso  = "$jobs_dir/job_vA-$curr_runnu.pbs_o.log";

    $grep_pipe     = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
    $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
    $evgen_cmd     = "gevgen -n $nev -s -e $en -p $nu -t $tgt -r $curr_runnu | grep_pipe &> $logfile_evgen";
    $conv_cmd      = "gntpc -f gst -i gntp.$curr_runnu.ghep.root | grep -B 100 -A 30 -i \"warn\\|error\\|fatal\" &> $logfile_conv";

    # create the PBS script
    #
    open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
    print PBS "#!/bin/bash \n";
    print PBS "#PBS -o $logfile_pbso \n";
    print PBS "#PBS -e $logfile_pbse \n";
    print PBS "source $genie_setup \n"; 
    print PBS "cd $jobs_dir \n";
    print PBS "export GSPLOAD=$xspl_file \n";
    print PBS "export GEVGL=$gevgl \n";
    print PBS "export GSEED=$mcseed  \n";
    print PBS "$evgen_cmd \n";
    print PBS "$conv_cmd \n";

    print "EXEC: $evgen_cmd \n";

    # submit job
    #
    `qsub -q $queue $batch_script`;
  }
}