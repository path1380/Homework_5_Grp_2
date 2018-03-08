#!/usr/apps/bin/perl
#This script adjusts the annulus area calculation
#parameters to create a plot for the convergence
#as we increase the number of gauss quadrature nodes
#per element.

$cmdFile="./Templates/problemsetup.f90.Template";
$outFile="./src/problemsetup.f90";

$num_degree = 1;
$num_quad_max = 10;
$num_theta = 2;
$num_r = 2;

for($q = 2;$q <= $num_quad_max; $q = $q+1){
	# Open the Template file and the output file.
	open(FILE,"$cmdFile") || die "cannot open file $cmdFile!" ;
	open(OUTFILE,"> $outFile") || die "cannot open file!" ;
	while($line = <FILE>){
		$line =~ s/\bNNNN\b/$q/;
		$line =~ s/\bQQQQ\b/$num_degree/;
		$line =~ s/\bNTNT\b/$num_theta/;
		$line =~ s/\bNRNR\b/$num_r/;
		print OUTFILE $line;
	}
	close( OUTFILE );
	close( FILE );
	system("make -f Makefile_line");
	system("./annulus_line.x >> output.txt");
	system("make -f Makefile_line clean");
}
system("matlab \"$@\" -nosplash -nodisplay < build/annulus_line.m");
system("rm output.txt");

exit
