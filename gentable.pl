#!/usr/bin/perl

use strict;
use warnings;

my @MECH_LABELS = qw(S2I SC SOP I2S IDX RF IF SUB ETA RID RMDC ITP);
my @MAN_LABELS  = qw(ITA EGPT O2L FT ET DKS FS);

my $PATH_PREFIX = shift @ARGV;
my @FILES = qw(lib/block_diag.mli   lib/block_diag.ml
               lib/cov_const.mli    lib/cov_const.ml
               lib/cov_lin_one.mli  lib/cov_lin_one.ml
               lib/cov_lin_ard.mli  lib/cov_lin_ard.ml
               lib/cov_se_iso.mli   lib/cov_se_iso.ml
               lib/cov_se_fat.mli   lib/cov_se_fat.ml
               lib/fitc_gp.mli      lib/fitc_gp.ml
               lib/interfaces.ml
               lib/gpr_utils.ml
               app/ocaml_gpr.ml);

die 'Specify directory for source files' unless defined $PATH_PREFIX;
$PATH_PREFIX .= '/';

my $data = count_labels();
print_table('tblmech.tex', $data, [@MECH_LABELS, 'mech'], [@MECH_LABELS, '\textbf{Total}']);
print_table('tblman.tex', $data, [@MAN_LABELS, 'man'], [@MAN_LABELS, '\textbf{Total}']);
print_table('tblall.tex', $data, ['lines', 'mech', 'man', 'all'], ['Lines', 'Mechanical', 'Manual', '\textbf{Total}']);

sub print_table {
    my($texfile, $data, $labels, $head) = @_;

    open(my $fh, '>', $texfile) or die "Cannot open \"$texfile\"";

    print $fh '\begin{tabular}{l', ('r' x scalar(@$labels)), "}\n  \\hline\n";
    print $fh '  & ', join(' & ', @$head), " \\\\\n  \\hline\n";

    my $print_row = sub {
        my($name, $entry) = @_;
        print $fh "  $name";
        print $fh ' & ', $$entry{$_} foreach @$labels;
        print $fh " \\\\\n";
    };

    foreach my $filename (@FILES) {
        my $fn = $filename;
        $fn =~ s/_/\\_/g;
        &$print_row("$fn", $$data{$filename});
    }

    print $fh "  \\hline\n";
    &$print_row("\\textbf{Total}", $$data{'total'});
    &$print_row("\\textbf{Percentage}", $$data{'percent'});
    print $fh "  \\hline\n  \\end{tabular}\n";

    close($fh);
}

sub count_labels {
    my $data = {};
    my $total = new_entry();
    foreach my $filename (@FILES) {
        my $entry = count_labels_in_file($PATH_PREFIX . $filename);
        $$data{$filename} = $entry;
        $$total{$_} += $$entry{$_} foreach keys %$entry;
    }

    my $percent = new_entry();
    $$percent{$_} = sprintf('%.2f', 100 * $$total{$_} / $$total{'lines'}) foreach keys %$total;

    $$data{'total'} = $total;
    $$data{'percent'} = $percent;

    return $data;
}

sub count_labels_in_file {
    my $filename = shift;

    open(my $fh, '<', $filename) or die "Cannot open \"$filename\": $!";

    my $entry = new_entry();

    while (defined(my $line = <$fh>)) {
        my @labels = ();
        while ($line =~ /\(\*!([^()]*)\*\)/g) {
            push(@labels, map
                 {
                     s/\s//g;      # remove white spaces
                     s/\[\d+\]$//; # remove a reference number
                     $_;
                 } split(/,/, $1));
        }

        #my %count;
        #@labels = grep {!$count{$_}++} @labels;

        # check labels
        foreach my $lb (@labels) {
            if (!grep { $_ eq $lb } (@MECH_LABELS, @MAN_LABELS)) {
                die "Unknown label `$lb' in \"$filename\"";
            }
        }

        # count
        $$entry{'mech'}++ if union_is_nonempty(\@labels, \@MECH_LABELS);
        $$entry{'man'}++ if union_is_nonempty(\@labels, \@MAN_LABELS);
        $$entry{'all'}++ if union_is_nonempty(\@labels, [@MECH_LABELS, @MAN_LABELS]);

        my %flags = ();
        $flags{$_} = 1 foreach @labels;
        $$entry{$_}++ foreach keys %flags;

        $$entry{'lines'}++;
    }

    close($fh);

    return $entry;
}

sub union_is_nonempty {
    my($a1, $a2) = @_;
    foreach my $x (@$a1) {
        foreach my $y (@$a2) {
            return 1 if $x eq $y;
        }
    }
    return 0;
}

sub new_entry {
    my $entry = { 'lines' => 0, 'mech' => 0, 'man' => 0, 'all' => 0 };
    $$entry{$_} = 0 foreach @MECH_LABELS, @MAN_LABELS;
    return $entry;
}

__END__

my @MECH_CHANGES = qw(S2I SC SOP IDX RF IF SUB ETA RID RMDC ITP);
my @MAN_CHANGES  = qw(ITA EGPT O2L I2S FT ET DKS FS);
my @KEYWORDS     = (@MECH_CHANGES, @MAN_CHANGES);

my @FILES = qw(lib/block_diag.mli   lib/block_diag.ml
               lib/cov_const.mli    lib/cov_const.ml
               lib/cov_lin_one.mli  lib/cov_lin_one.ml
               lib/cov_lin_ard.mli  lib/cov_lin_ard.ml
               lib/cov_se_iso.mli   lib/cov_se_iso.ml
               lib/cov_se_fat.mli   lib/cov_se_fat.ml
               lib/fitc_gp.mli      lib/fitc_gp.ml
               lib/interfaces.ml
               lib/gpr_utils.ml
               app/ocaml_gpr.ml);

my $HEADER = [ 'Filename', $LINES,
               @MECH_CHANGES, '_Mech_',
               @MAN_CHANGES, '_Man_',
               '_Total_' ];
my $MASKS  = [ (map { keyword_mask($_) } @MECH_CHANGES), keyword_mask(@MECH_CHANGES),
               (map { keyword_mask($_) } @MAN_CHANGES), keyword_mask(@MAN_CHANGES),
               keyword_mask(@MECH_CHANGES, @MAN_CHANGES) ];
my @CNT_LL = count_keywords($MASKS, \@FILES);

print_table_as_gfm($HEADER, \@CNT_LL);

exit;

sub print_table_as_tex {
    my($header, $cnt_lst_lst) = @_;

    my $print_row = sub {
        my $row = shift;
        printf "| %s ", $$row[$_] foreach 0 .. $#{$row};
        print "|\n";
    };

    # Print the header
    &$print_row($header);
}

sub print_table_as_gfm {
    my($header, $cnt_lst_lst) = @_;

    # Emphasize items in columns.
    foreach my $row (@$cnt_lst_lst) {
        foreach my $i (0 .. $#$row) {
            $$row[$i] = '_' . $$row[$i] . '_' if $$header[$i] =~ /^_/;
        }
    }

    # Compute the width of each column.
    my @colw = (0) x scalar(@$header);
    foreach my $row ($header, @$cnt_lst_lst) {
        foreach my $i (0 .. $#$row) {
            my $len = length $$row[$i];
            $colw[$i] = $len if $colw[$i] < $len;
        }
    }

    my $print_row = sub {
        my $row = shift;
        printf "| %*s ", $colw[$_], $$row[$_] foreach 0 .. $#colw;
        print "|\n";
    };

    # Print the header
    &$print_row($header);

    # Print a separator
    print '|-', '-' x $colw[$_], ($_ == 0 ? '-' : ':') foreach 0 .. $#colw;
    print "|\n";

    # Print rows
    &$print_row($_) foreach @$cnt_lst_lst;
}

sub count_keywords {
    my($mask_lst, $prefix, $filename_lst) = @_;
    my @cnt_lst_lst = ();

    foreach my $filename (@$filename_lst) {
        my @cnt_lst = count_keywords_in_file($mask_lst, $prefix . $filename);
        push(@cnt_lst_lst, \@cnt_lst);
    }

    # Initialize
    my @total = (0) x (scalar(@$mask_lst) + 1);
    my @percentage = @total;

    # Total
    foreach my $i (0 .. $#cnt_lst_lst) {
        $total[$_] += $cnt_lst_lst[$i][$_] foreach 0 .. $#total;
    }
    push(@cnt_lst_lst, \@total);

    # Percentage
    $percentage[$_] = sprintf('%.1f', $total[$_] / $total[0] * 100) foreach 0 .. $#total;
    push(@cnt_lst_lst, \@percentage);

    # Labels of rows
    my @lb = (@$filename_lst, '**Total**', '**Percentage**');
    unshift(@{$cnt_lst_lst[$_]}, $lb[$_]) foreach 0 .. $#cnt_lst_lst;

    return @cnt_lst_lst;
}

sub count_keywords_in_file {
    my($mask_lst, $filename) = @_;

    open(my $fh, '<', $filename) or die "Cannot open \"$filename\"";

    my @cnt_lst = (0) x scalar(@$mask_lst);
    my $lnum = 0;
    foreach my $line (<$fh>) {
        ++$lnum;

        eval {
            foreach my $i (0 .. $#$mask_lst) {
                ++$cnt_lst[$i] if $$mask_lst[$i] & line_to_mask($line);
            }
        };
        print STDERR "$filename:$lnum: $@\n" if $@;
    }

    close($fh);

    return ($lnum, @cnt_lst);
}

sub line_to_mask {
    my $line = shift;
    my $mask = 0;

    while ($line =~ /\(\*!([^()]*)\*\)/g) {
        foreach my $str (split(/,/, $1)) {
            if ($str =~ /([0-9A-Za-z\-]+)(?:\[\d+\])?/) {
                $mask |= keyword_mask($1);
            } else {
                die "Unknown keyword `$str'";
            }
        }
    }

    return $mask;
}

sub keyword_mask {
    my $mask = 0;

    foreach my $str (@_) {
        my @i = grep { $str eq $KEYWORDS[$_] } (0 .. $#KEYWORDS);
        die "Unknown keyword `$str'" unless scalar @i; # If @i is empty, $str is an illegal keyword.
        $mask |= 1 << (shift @i);
    }

    return $mask;
}
