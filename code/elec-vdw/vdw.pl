#!/usr/bin/perl

use strict;
my $lis; 
#my $infile=shift;

###################################
#chg & wdw
my %chg;
my %dvdw1;
my %dvdw;
my %evdw1;
my %evdw;
my %moltype;
my $resid;
my $atomid;
my $atomtype;
open(vdw,"./code/elec-vdw/vdw.dat")||die $!;
while(<vdw>){
        chomp $_;
	$atomtype=substr($_,0,2);

        $dvdw1{"$atomtype"}=substr($_,12,6);
        $evdw1{"$atomtype"}=substr($_,20,6);

}
close vdw;

open(elec,"./code/elec-vdw/elec.dat")||die $!;
while(<elec>){
	chomp $_;
#	my @current=split(/ +/,$_);
	 $resid=substr($_,0,3);
	 $atomid=substr($_,9,3);	
	 $atomtype=substr($_,15,2);
	$moltype{"$resid $atomid"}=substr($_,81,1);
	$chg{"$resid $atomid"}=substr($_,67,7);
	$dvdw{"$resid $atomid"}=$dvdw1{"$atomtype"};
	$evdw{"$resid $atomid"}=$evdw1{"$atomtype"};
	
}
close elec;


###################################
#get info
open(FILE,"./code/elec-vdw/lis")||die $!;
while(<FILE>){

chomp $_;
$lis=$_;

#pdb_info and elec info -in
my %cor_x;
my %cor_y;
my %cor_z;
my %type;
my %mol1;
my %mol2;
my %elec;
my %dvdwm;
my %evdwm;
my $atomno;
my $resno;




open(pdb,"./feature/$lis")||die "Can't open:$lis\n";
while(<pdb>)
{
	chomp $_;
#	my @atom=split(/ +/,$_);
#	$type{"$atom[1] $atom[4]"}=$atom[2];
    if(/ATOM/){
	 $atomno=substr($_,6,5);
	 $resno=substr($_,21,1);
	 $atomid=substr($_,13,3);
	 $resid=substr($_,17,3);
	$elec{"$atomno $resno"}=$chg{"$resid $atomid"};
	$dvdwm{"$atomno $resno"}=$dvdw{"$resid $atomid"};
	$evdwm{"$atomno $resno"}=$evdw{"$resid $atomid"};
	$type{"$atomno $resno"}=$moltype{"$resid $atomid"};

	if($type{"$atomno $resno"} eq "p"){
	$mol1{"$atomno $resno"}=substr($_,17,3);
	}else{
	$mol2{"$atomno $resno"}=substr($_,17,3);
	}

	$cor_x{"$atomno $resno"}=substr($_,30,8);
	$cor_y{"$atomno $resno"}=substr($_,38,8);
	$cor_z{"$atomno $resno"}=substr($_,46,8);
}
}

close pdb;
###############################################
#score

my $Score=0;
my $snscore=0;
my $spscore=0;
my $lnscore=0;
my $lpscore=0;
my $r_sque;
my $vdwattr=0;
my $vdwrep=0;
my $ID;
my $id;

foreach $ID(keys %mol1){
	foreach $id(keys %mol2){
	
		$r_sque=($cor_x{"$ID"}-$cor_x{"$id"})**2+($cor_y{"$ID"}-$cor_y{"$id"})**2+($cor_z{"$ID"}-$cor_z{"$id"})**2;
		my $r=$r_sque**(1/2);
		my $d=($dvdwm{"$ID"}+$dvdwm{"$id"})/2;
#ele
		if($r_sque < 25){
			if($r_sque<4) {$r_sque=4;}
			if($elec{"$ID"}*$elec{"$id"}<0){
				$snscore=$snscore+83*$elec{"$ID"}*$elec{"$id"}/$r_sque;
			}
			if($elec{"$ID"}*$elec{"$id"}>0){
				$spscore=$spscore+83*$elec{"$ID"}*$elec{"$id"}/$r_sque;
			}
		}
		if($r_sque > 25 && $r<16){
			if($elec{"$ID"}*$elec{"$id"}<0){
				$lnscore=$lnscore+83*$elec{"$ID"}*$elec{"$id"}/$r_sque; 
			}
			if($elec{"$ID"}*$elec{"$id"}>0){
				$lpscore=$lpscore+83*$elec{"$ID"}*$elec{"$id"}/$r_sque;
			}
		}
		$Score=$Score+83*$elec{"$ID"}*$elec{"$id"}/$r_sque;
#vdw
		if($r<8){
			if($r<0.89*$d){$vdwrep =$vdwrep +(10-11.2*$r/$d);     }		
			if($r>0.89*$d){$vdwattr=$vdwattr+(($evdwm{"$ID"}*$evdwm{"$id"})**(1/2))*(($d/$r)**12-2*($d/$r)**6); }
		}
	}
}

	 print "$lis       $Score $snscore $spscore $lnscore $lpscore     $vdwrep $vdwattr \n";
}	
