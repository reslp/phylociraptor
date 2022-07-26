#!/bin/bash

file=$(realpath $1)

status=$(echo -e "$(pwd)\n$file" | sed 's|/| |g' | perl -ne '
$home=$_;
$from=<>;
chomp;
chomp $from;
@home=split(" ", $home);
@from=split(" ", $from);
for ($i=0; $i<@home; $i++) {
#	print "$home[$i] vs. $from[$i]\n";
	if ($home[$i] ne $from[$i]) {
		print "1"; exit; 
	}
}; print "0"')

echo $status
