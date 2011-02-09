while(<>){
	last if /segID/
}

open VERT, ">data/ground_vertices.txt" or die $!;
open EDGE, ">data/ground_edges.txt"    or die $!;

$_=<>;
while(<>){
	chomp;
	last if /^$/;
	$tmp = <>;   # line containing only "0"
	$@ = split/\s+/;
	$id  = $_[0] - 1;
	$par = $_[4] - 1;
    @a = map{sprintf "%3.6f", $_}($_[1], $_[2], $_[3]);
	print VERT join(" ", (@a, int(rand()*2)+1, 1,1,1) ), "\n";
	print EDGE "$id $par\n" if $par>=0;
}
