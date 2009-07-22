#---------------------------------------------------------------------------
# Submit GENIE/nd280 event generation job using the full geometry and a 
# JNUBEAM ntuple-based flux description.
# 
# For use at the RAL/PPD Tier2 PBS batch farm.
#
# Syntax:
# shell$ perl submit-evg_nd280.pl <options>
#
# Options:
#  --fluxrun       : Input flux run number
#  --version       : GENIE version
# [--production]   :
# [--cycle]        :
# [--use-valgrind] : 
#
# Example:
# shell$ perl submit-evg_nd280.pl --fluxrun 180 --version v2.4.0  
#                                 --production mdc0 --cycle 01
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#---------------------------------------------------------------------------

#!/usr/bin/perl

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--fluxrun')       { $fluxrun       = $ARGV[$iarg+1]; }
  if($_ eq '--version')       { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')  { $use_valgrind  = $ARGV[$iarg+1]; }
  $iarg++;
}

die("** Aborting [Undefined flux file run #. Use the --fluxrun option]")
unless defined $fluxrun;
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$use_valgrind = 0           unless defined $use_valgrind;
$production   = "nd280test" unless defined $production;
$cycle        = "01"        unless defined $cycle;

$job_pot        = "1E+18";
$mcrun_base     = 10000000;
$mcseed_base    = 210921029;
$batch_queue    = "prod";
$time_limit     = "30:00:00";
$genie_inst_dir = "/opt/ppd/t2k/GENIE";
$production_dir = "/opt/ppd/t2k/GENIE/scratch";
$inputs_dir     = "/opt/ppd/t2k/GENIE/data/job_inputs";
$genie_setup    = "$genie_inst_dir/$genie_version-setup";
$job_dir        = "$production_dir/nd280-$production\_$cycle";
$flux_dir       = "$inputs_dir/t2k_flux/07a/nd";
$flux_file_prfx = "nu.nd280.";
$flux_file_sufx = ".root";
$flux_file      = "$flux_dir/$flux_file_prfx$fluxrun$flux_file_sufx";
$flux_det_loc   = "nd5";
$geom_file      = "$inputs_dir/t2k_geom/ND280.root";
$xspl_file      = "$inputs_dir/xspl/gxspl-t2k-$genie_version.xml";
$geom_dunits    = "clhep_def_density_unit";
$geom_lunits    = "mm";
$file_prefix    = "genie_nd280";
$mcrun          = $mcrun_base  + $fluxrun; 
$mcseed         = $mcseed_base + $fluxrun;

die("** Aborting [Can not find GENIE setup script: ....... $genie_setup]") 
unless -e $genie_setup;
die("** Aborting [Can not find flux file: ................ $flux_file]") 
unless -e $flux_file;
die("** Aborting [Can not find geometry file: ............ $geom_file]") 
unless -e $geom_file;
die("** Aborting [Can not find cross sections file: ...... $xspl_file]") 
unless -e $xspl_file;

print "@@@ Will submit job with MC run number = $mcrun (seed number = $mcseed)\n";

# make the jobs directory
#
mkpath ($job_dir, {verbose => 1, mode=>0777}); 
die("$job_dir doesn't exist") unless -d $job_dir;

$batch_script  = "$job_dir/nd280job-$mcrun.pbs";
$logfile_evgen = "$job_dir/nd280job-$mcrun.evgen.log";
$logfile_conv  = "$job_dir/nd280job-$mcrun.conv.log";
$logfile_pbse  = "$job_dir/nd280job-$mcrun.pbs_e.log";
$logfile_pbso  = "$job_dir/nd280job-$mcrun.pbs_o.log";
$ghep_file     = "$file_prefix.$production\_$cycle.$mcrun.ghep.root";
$grep_pipe     = "grep -B 50 -A 30 -i \"warn\\|error\\|fatal\" ";
$valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
$evgen_cmd     = "gT2Kevgen -g $geom_file -f $flux_file,$flux_det_loc -r $mcrun -L $geom_lunits -D $geom_dunits -E $job_pot | $grep_pipe &> $logfile_evgen";
$frenm_cmd     = "mv gntp.$mcrun.ghep.root $ghep_file";
$fconv_cmd     = "gntpc -f t2k_rootracker -i $ghep_file | $grep_pipe &> $logfile_conv";

# create the PBS script
#
open(PBS, ">$batch_script") or die("Can not create the PBS batch script");

print PBS "#!/bin/bash \n";
print PBS "#PBS -N $mcrun\_nd280-$production-$cycle \n";
print PBS "#PBS -l cput=$time_limit \n";
print PBS "#PBS -o $logfile_pbso \n";
print PBS "#PBS -e $logfile_pbse \n";
print PBS "source $genie_setup \n";
print PBS "cd $job_dir \n";
print PBS "export GSPLOAD=$xspl_file \n";
print PBS "unset GEVGL \n";
print PBS "export GSEED=$mcseed \n";
print PBS "$evgen_cmd \n";
print PBS "$frenm_cmd \n";
print PBS "$fconv_cmd \n";

print "@@@ exec: $evgen_cmd \n";

# submit job
#
`qsub -q $batch_queue $batch_script`;
