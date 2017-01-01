#!/usr/local/bin/perl

$script = "ttHbb_data2mc_treeReader";       # Name of job

$workingDir = "/uscms_data/d2/dpuigh/TTH/Run2015_76x/CMSSW_7_6_3/src/TriggerRun2/TriggerAnalyzer";

$sample = 2510;
$intLumi = 2430;
$tthf   = -1;
$Njobs  = 1;
$jobN   = 1;
$useHTbins = 0;

$num = @ARGV;

if( $num >= 1 ){
  $sample = $ARGV[0];

  if( $num>=2 ){
    $Njobs = $ARGV[1];

    if( $num>=3 ){
      $jobN = $ARGV[2];

      if( $num>=4 ){
	$intLumi = $ARGV[3];
      
	if( $num>=5 ){
	  $tthf = $ARGV[4];

	  if( $num>=6 ){
	    $useHTbins = $ARGV[5];
	  }
	}
      }
    }
  }
}


open SHFILE, "> Script/condor\_$script.sh";
print SHFILE "#!/bin/sh\n";
print SHFILE "\n";
print SHFILE "echo $1\n";
print SHFILE "\n";
print SHFILE "echo \"\"\n";
print SHFILE "echo \"Using ROOT on Condor\"\n";
print SHFILE "echo \"\"\n";
#print SHFILE "cd \${_CONDOR_SCRATCH_DIR}\n";
print SHFILE "\n";
print SHFILE "sample=\$1\n";
print SHFILE "NumEvents=\$2\n";
print SHFILE "NumJobs=\$3\n";
print SHFILE "jobN=\$4+1\n";
print SHFILE "intLumi=\$5\n";
print SHFILE "tthf=\$6\n";
print SHFILE "useHTbins=\$7\n";
print SHFILE "\n";
print SHFILE "root -b -q $workingDir/macros/head.C '$workingDir/macros/$script.C+('\$sample','\$NumEvents','\$NumJobs','\$jobN','\$intLumi','\$tthf','\$useHTbins',1)'\n";
print SHFILE "\n";
close SHFILE;


open CONDORFILE, "> Condor/condor\_$script.jdl";
print CONDORFILE "# A Condor submission file\n";
print CONDORFILE "Executable              = Script/condor\_$script.sh\n";
print CONDORFILE "Universe                = vanilla\n";
print CONDORFILE "Getenv                  = true\n";
print CONDORFILE "\n";
print CONDORFILE "Arguments               = $sample -1 $Njobs $jobN $intLumi $tthf $useHTbins\n";
print CONDORFILE "Output                  = Output/condor\_$sample\_$script\_$jobN\_$tthf.out\n";
print CONDORFILE "Error                   = Error/condor\_$sample\_$script\_$jobN\_$tthf.err\n";
print CONDORFILE "Log                     = Log/condor\_$sample\_$script\_$jobN\_$tthf.log\n";
print CONDORFILE "\n";
print CONDORFILE "use_x509userproxy = true\n";
print CONDORFILE "Should_Transfer_Files   = YES\n";
print CONDORFILE "When_To_Transfer_Output = ON_EXIT\n";
print CONDORFILE "Transfer_Input_Files = CSVv2.csv, CSVv2_TagCountTT.csv, ttH_BTV_CSVv2_13TeV_2015D_20151122.csv, ttH_BTV_CSVv2_13TeV_2015D_2016_02_08.csv, ttH_BTV_CMVAv2_13TeV_2015D_2016_02_12.csv\n";
print CONDORFILE "\n";
print CONDORFILE "#+IsLocalJob             = true\n";
print CONDORFILE "#Rank                    = TARGET.IsLocalSlot\n";
print CONDORFILE "\n";
print CONDORFILE "Queue 1\n";
print CONDORFILE "\n";
close CONDORFILE;

system("chmod a+x Script/condor\_$script.sh");
print "submitting: condor_submit Condor/condor\_$script.jdl\n";
system("condor_submit Condor/condor\_$script.jdl");

