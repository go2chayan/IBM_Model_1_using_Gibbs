#! /usr/bin/perl

# first argument is .afile format system output
# second argument is gold standard with "sure" and "possible" alignments

# .afile looks like this:
# begin 1
# 1 1
# 2 2 3
# 3 4
# end
# begin 2
# ...

open(OUT, shift) || die;
open(GOLD, shift) || die;

while(<OUT>){
    if (/begin (\d+)/) {
	$sent = $1;
    } elsif (/end/) {
    } else {
	($e, @c) = split;
	for (@c) {
	    $out{$sent}{$e}{$_}++;
	    $out_atomic{"$sent-$e-$_"}++;
	    $found_n++;
	}
    }
}

while(<GOLD>){
    ($sent, $e, $f, $conf) = split;
    $sent+=0;	# get rid of leading 0's
    next unless ($e || $f);	# not scoring insertions/deletions
    if ($conf eq 'S') {
	$true_n++;
	if ($out_atomic{"$sent-$e-$f"}) {
	    $recall_n++;
	}
    }
    if ($out_atomic{"$sent-$e-$f"}) {
	$prec_n++;		# either sure or possible is OK for precision
    }
}

printf "precision: %5d / %5d = %5.1f%%\n", $prec_n, $found_n, 100* $prec_n / $found_n;
printf "recall:    %5d / %5d = %5.1f%%\n", $recall_n, $true_n,  100* $recall_n / $true_n;
printf "aer:                     = %7.2f\n",    1 - ($prec_n + $recall_n) / ( $found_n + $true_n ) ;
