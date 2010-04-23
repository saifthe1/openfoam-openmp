#!/usr/bin/perl

use Digest::MD5 qw(md5_hex);

my $string;

foreach(<*>) {
	if($_ =~ /^[0-9.e\-]+$/ and !($_ eq '0')) {
		foreach(<$_/*>) {
			local $/=undef;
			open FILE, "$_" or die "Couldn't open file: $!";
			binmode FILE;
			$string .= <FILE>;
			close FILE;
		}
	}
}

print md5_hex($string);
