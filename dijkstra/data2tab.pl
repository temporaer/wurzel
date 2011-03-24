#!/usr/bin/perl -w
use Data::Dumper;

sub get_name{
	my $fn = shift;
	$fn =~ /(\d+)x(\d+)x(\d+)/;
	my $res = "$1\\times$2\\times$3";
	if ($fn =~ /\d+_(\w+)$/)
	{
		return "\$$res-\\mathrm{$1}\$";
	}
	return "\$$res\$";
}
sub get_hash{
	my $str = shift();
	my %h;
	$str =~  /Total mass: ([\d.]+)\s*$/mg;
	$h{mass} = $1/1000.0;
	$str =~ m|V/H   mass: ([\d.]+)\s*$|mg;
	$h{vh_mass} = $1;
	$str =~  /Total len: ([\d.]+)\s*$/mg;
	$h{len} = $1/1000.0;
	$str =~ m|V/H   len: ([\d.]+)\s*$|mg;
	$h{vh_len} = $1;
	return %h;
}
sub get_distance{
	my $str = shift();
	if($str =~ /Average distance: ([\d.]+) var:([\d.]+)/)
	{
		return sprintf("\$%2.2f\\pm%2.2f\$", $1,$2);
	}
	return "\$0.00\$"
}


my @bases;
#push @bases, "../data/GersteLA_256x256x410";
push @bases, "../data/GersteLA_192x192x410_normal";
#push @bases, "../data/GersteLA_128x128x410";
#push @bases, "../data/GersteLA_96x96x410";
#push @bases, "../data/GersteLA_64x64x410";
push @bases, "../data/GersteLA_192x192x410_noise2";
push @bases, "../data/GersteLA_192x192x410_noise3";
#push @bases, "../data/GersteLA_192x192x410_noise4";
#push @bases, "../data/GersteLA_192x192x410_noise5";
my @hashes;

foreach my $fn (@bases){

	my $ret = `./treecompare print $fn`;
	my %h   = get_hash($ret);
	$h{id}  = get_name($fn);

	if($fn !~ /_normal/){
	        $ret = qx(./treecompare distance $fn ../data/GersteLA_192x192x410_normal);
	}
	$h{dist} = get_distance($ret);
	print Dumper(\%h);
	push @hashes, \%h;
}

foreach my $h (@hashes){
	my %d = %$h;
	my @l = map{sprintf("\$%2.2f\$", $d{$_})}(qw/len mass vh_mass vh_len/);
	print "$d{id} & ", join(" & ",@l), " & $d{dist}\\\\\n"
}
