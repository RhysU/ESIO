#!/usr/bin/env perl
# --------------------------------------------------------------------
# Simple support script for license management with PECOS developed
# software. Intended for use with autotools based projects. In
# particular, if a ".in" version of an input src file exists, the .in
# version will be updated in favor of named src file..
#
# Originally: August 2010 (ks)
# --------------------------------------------------------------------
# Id: update_license.pl 13675 2010-09-19 17:50:33Z karl
# --------------------------------------------------------------------


#use warnings;
use File::Compare;
use File::Basename;
use Getopt::Long;

my $license_begin_delim='----------bl-$';
my $license_end_delim='----------el-$';

# Command-line parsing

GetOptions("S=s@","c2f_comment") || die "Error using GetOptions";

if (@ARGV >= 2) {
    $license_file = shift@ARGV;
} else {
    print "\nUsage: update_license.pl [OPTION] LICENSE-FILE SOURCE-FILES...\n\n";
    exit 0;
}

if (@opt_S) {
    @top_srcdir = @opt_S;
 }

# Verify license file
# existence and cache contents.

if ( ! -s $license_file ) {
    print "\n** Error: unable to open valid license file ($license_file).\n\n";
    exit(1);
}

open ($LIC,$license_file) || die "Cannot open $license_file\n";
my @license_text = <$LIC>;
close($LIC);

# Multiple language support; in certain instances, we have multiple
# languages within the same project (eg. libgrvy and libmasa which
# provide interfaces to both C and Fortran). We would like to trigger
# the license off of a single file only and thus convert the comment
# syntax here if desired.

# Convert C/C++ comment (//) to Fortran (!!)

if ($opt_c2f_comment) {

    my @license_tmp=();

    foreach $line (@license_text) {
	if ($line =~ s/^\/\//!!/) {
	    push(@license_tmp,$line);
	} else {
	    push(@license_tmp,$line);
	}
    }

    @license_text=@license_tmp;
}

# Scan all provided input files and look for designated license
# begin/end delimiters.  When found, update with provided license file
# contents.

my $found_delim = 0;

while (@ARGV)
{
    $found_delim = 0;

    $infile_test = shift @ARGV;

    if(@opt_S) {
	$infile_test = "@top_srcdir/$infile_test";
    }

    # autoconf support - look for ".in" version of the src file

    if ( -e "$infile_test.in" ) {
	$infile = "$infile_test.in";
    } else {
	$infile = "$infile_test";
    }

    open($IN, "<$infile") || die "Cannot open $infile\n";

    my $basename = basename("$infile");
    my $dirname  = dirname("$infile");
    my $tmpfile  = "$dirname/.$basename.tmp";

    open ($TMPFILE,">$tmpfile") || die "Cannot create tmp file $tmpfile";

    while (<$IN>) {
	if(/$license_begin_delim/../$license_end_delim/) {
	    $found_delim=1;
	    if (/--bl-$/) {
		print $TMPFILE @license_text;
	    }
	} else {
	    print $TMPFILE $_;
	}
    }

    close($IN);
    close($TMPFILE);

    if( $found_delim ) {
	if ( compare($infile,$tmpfile) != 0 )  {
	    print "[license_tool]: updating license in file $infile\n";
	    rename($tmpfile,$infile) || die "Cannot rename updated file\n";
	} else {
	    unlink($tmpfile) || die "Unable to remove temporary file\n";
	}
    } else {
#	print "[license_tool]: warning: no license delimiters found in file $infile\n";
	unlink($tmpfile) || die "Unable to remove temporary file\n";
    }
}






