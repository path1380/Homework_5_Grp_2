#!/usr/apps/bin/perl
#This script adjusts the annulus area calculation
#parameters to create a plot for the convergence
#as we increase the number of gauss quadrature nodes
#per element.

$cmdFile="./Templates/problemsetup.f90.Template";
$outFile="./src/problemsetup.f90";

$num_quad_max = 30;
$num_theta = 2;
$num_r = 2;

for($q = 2;$q <= $num_quad_max; $q = $q+1){
	# Open the Template file and the output file.
	open(FILE,"$cmdFile") || die "cannot open file $cmdFile!" ;
	open(OUTFILE,"> $outFile") || die "cannot open file!" ;
	while($line = <FILE>){
		$line =~ s/\bQQQQ\b/$q/;
		$line =~ s/\bNNNN\b/$q+4/;
		$line =~ s/\bNTNT\b/$num_theta/;
		$line =~ s/\bNRNR\b/$num_r/;
		print OUTFILE $line;
	}
	close( OUTFILE );
	close( FILE );
	system("make -f Makefile_coeffs");
	system("./coeff2d.x >> output.txt");
	system("make -f Makefile_coeffs clean");
}
system("matlab \"$@\" -nosplash -nodisplay < build/coeff_plot.m");
system("rm output.txt");

exit
