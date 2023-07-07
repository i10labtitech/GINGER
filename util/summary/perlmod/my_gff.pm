# Copyright (C) 2018 Itoh Laboratory, Tokyo Institute of Technology
# 
# This file is part of GINGER.
# 
# GINGER is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# GINGER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with GINGER; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

package my_gff;

use Exporter;

@ISA = qw(Exporter);
@EXPORT = qw/gff3_tree_set_relations gff3_tree_get_parent_id gff3_tree_get_root_id gff3_attr_get_id/;

sub gff3_tree_set_relations
{
	my ($in, $source, $feat, $st, $ed, $attr, $id, $parent, %tree);

	my $tree_ref = shift;
	my $file_name = shift;

	open($in, $ARGV[0]);
	while (chomp($l = <$in>)) {
		($source, $feat, $st, $ed, $attr) = (split(/\t/, $l))[1, 2, 3, 4, 8];
		if ($attr =~ /ID=([^;]+)/) {
			$id = $1;
		}
		elsif ($attr =~ /Name=([^;]+)/) {
			$id = $1;
		}
		else {
			next;
		}

		if ($attr =~ /Parent=([^;]+)/) {
			$parent = $1;
		}
		else {
			$parent = '';
		}

		$tree_ref->{$id} = [split(/,/, $parent)];
	}
	close $in;
}

sub gff3_tree_get_parent_id
{
	my $tree_ref = shift;
	my $id = shift;

	return @{$tree_ref->{$id}};
}

sub gff3_tree_get_root_id
{
	my $tree_ref = shift;
	my $id = shift;

	while ($tree_ref->{$id}[0]) {
		$id = $tree_ref->{$id}[0];
	}

	return $id;
}

sub gff3_attr_get_id
{
	my $attr = shift;
	my $id;
	if ($attr =~ /ID=([^;]+)/) {
		$id = $1;
	}
	elsif ($attr =~ /Name=([^;]+)/) {
		$id = $1;
	}
	return $id;
}

1;
