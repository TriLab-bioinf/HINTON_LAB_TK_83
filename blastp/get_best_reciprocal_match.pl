#!/usr/local/bin/perl

# input:
# NP_414543.1	820	WP_001264707.1	820	100.000	820	0	0	1	820	1820   	0.0	1687	100
#
  
while(<>){
    chomp;
    next if m/^#/; 
    
    @x=split /\t/; 
        
    if($cov{$x[0]} < $x[14] && 
        $x[4] > 30 && 
        $x[14] > 90 && 
        $sjt{$x[2]} < ($x[4] + $x[14])   
    ){ 
           print "$_\n"; 
           $sjt{$x[2]} = $x[4] + $x[14]; 
           
           # reset previous hit of subject
           $q = $qhit{$x[2]};
           print "$ _\n"
           delete( $shit{$q} );  

           $cov{$x[0]}   = $x[14]; 
           $shit{$x[0]}  = $x[2]; 
           $pid{$x[0]}   = $x[4]; 
           $qhit{$x[2]}  = $x[0]; 
    } 
} 

foreach $sacc (keys %shit){
        print "$sacc\t$shit{$sacc}\t$cov{$sacc}\t$pid{$sacc}\n"
} 
