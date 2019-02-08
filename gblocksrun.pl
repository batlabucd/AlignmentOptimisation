use strict;
use warnings;
my @array=(<*fas>);#list of all genes ending in "fas", as per output of prank
foreach my $file(@array){#loop through all files
    print $file."\n";
    my $gblockedfile=$file."\.out\.txt";#name of file for gblocks
    my $prot=""; my $nuc="";#protein and nucleotode file names
    if($file=~m/([\S]+\.prank)/){#get original name of prank file, asusme is form "[filename].prot.prank.aln.best.fas
      	$prot=$1;
	print "Prot file: ".$prot."\n";
    }
    if($file=~m/([\S]+)\.prot\.prank\.aln/){#get original nucleotide file
        $nuc=$1;
	print "Nuc file: ".$nuc."\n";
    }
    my $codons=$nuc."\.codons";#name of codon files
    print $nuc."\n".$prot."\n";
    if(-e $nuc){#check if nucleotode file exists
    }
    else{
	print $nuc." and ".$prot." not found. These files need to be in the same directory as ".$file."\n";
	die;
    }
    `./Gblocks $file -t=p -b1=2 -b2=43 -b3=2 -b4=2 -b5=h -d=y -k=y -p=t -e=.out`;#run gblocks initially
    my $o=$file."\.out\.txt";#filename for output file
    print $o."\n";
    open(IN, "$o");#open in Gblocks output file
    my $b2=2;
    while(<IN>){
	if($_=~m/Minimum[\s]+Number[\s]+Of[\s]+Sequences[\s]+For[\s]+A[\s]+Conserved[\s]+Position\:[\s]+([0-9]+)/){# read parameters from gblocks
	    $b2=$1;
	}
    }
    `./Gblocks $file -t=p -b1=2 -b2=$b2 -b3=2 -b4=2 -b5=h -d=y -k=y -p=t -e=.out`;#run gblocks again with maximum best estimated parameters   
    if(-e $gblockedfile){#if gblocks has ran successfully
	open(IN2, "$o");#open in blocks file
	my $fileforpal2nal=$o."\.marked";#create marked file
	open(PRINT, ">>$fileforpal2nal");#create file for pal2nal
	my @lines=<IN2>;
	@lines=split(/blue\)/, join('',@lines));#split gblock file
	@lines=split(/Parameters/,$lines[1]);
	my @data=split(/\n/,$lines[0]);
	foreach my $x(@data){
	    if($x=~m/(^Gblocks[\s]*)/){
		my $pad=$1;
		my $lin=$x;
		$lin=~s/G/ /;$lin=~s/b/ /;$lin=~s/l/ /;$lin=~s/o/ /;$lin=~s/c/ /;$lin=~s/k/ /;$lin=~s/s/ /;
		$lin=~s/\./ /g;#remove gblocks header, important for formatting output file
		print PRINT $lin."\n\n";
	    }
	    else{
		if($x=~m/\=\=/ || $x=~m/[0-9]+[\s]+[0-9]+/){#ignor this specific line
		}
		else{
		    if($x=~m/[\S]+/){#print out marked data
			print PRINT $x."\n";
		    }
		}
	    }
	}
	close IN2;
	close PRINT;
	open(FILE, "$fileforpal2nal");
	my @names=[];my $namecount=0;
	while(<FILE>){#Order DNA original file according to prank output order
	    if($_=~m/^([\S]+)/){
		my $header=$1;
		if($header~~@names){
		}
		else{
		    $names[$namecount]=$header;
		    $namecount++;
		}
	    }
	}
	close FILE;
	my $newnuc="Sorted_".$nuc;#create new 'sorted' nucleotode file, relative to the outout from prank
	if(-e $newnuc){
	    `rm $newnuc`;
	}
	#Code to order Nuc data based on order in prank
	open(OUT2, ">>$newnuc");
	open(IN2, "$nuc");
	my @nuc=(<IN2>);
	@nuc=split(/\>/,join('',@nuc));
	for(my $i=0;$i<scalar(@names);$i++){
	    my $header=$names[$i];
	    foreach my $target(@nuc){
		my $tarhead=substr($target,0,length($header));
		if($target=~m/^([\S]+)/){
		    my $spec_name_full=$1;
		    if($header eq $tarhead){
			open(NAMING, ">>names_to_fix");
			print NAMING $header."=>".$spec_name_full."\n";
			close NAMING;
			print OUT2 ">".$target."\n";
		    }
		}
	    }
	}
	close OUT2;
	`perl pal2nal.pl $fileforpal2nal $newnuc -blockonly -output paml > $codons`;#run pal2nal
	print "perl pal2nal.pl $fileforpal2nal $newnuc -blockonly -output paml > $codons\n";
	##need to fix shortened headers
	open(FIX, "names_to_fix");
	while(<FIX>){
	    if($_=~m/([\S]+)\=\>([\S]+)/){
		if ($1 eq $2){
		}
		else{
		    system("perl -pi -w -e 's/$1/$2/g;' $codons"); 
		}
	    }
	}
	system("rm names_to_fix");
    }
    else{
	print "Gblocked file ".$gblockedfile." cannot be found, gblocks has not ran for ".$file."\n"; 
	die;
    }
}
