# -*- perl -*-

package Math::Integral::Romberg;
require Exporter;
@ISA=qw(Exporter);
@EXPORT_OK=qw(integral return_point_count);
$VERSION = "0.01";
$abort = 0;
$return_point_count = 0;

use strict;

=head1 NAME

  Math::Integral::Romberg - single-variable numerical integration

=head1 SYNOPSIS

    use Math::Integral::Romberg 'integral';
    $area = integral(\&f, $x1, $x2);    # Short form
    $area = integral                    # Long form
     (\&f, $x1, $x2, $rel_err, $abs_err, $max_split, $min_split);

=head1 DESCRIPTION

integral() numerically estimates the integral of f() using Romberg's
method, a faster relative of Simpson's method.

=head2 Parameters

=over

=item $f

A reference to the function to be integrated.

=item $x1

=item $x2

The two extreme values of the range to be integrated.  C<&$f(x)> must be
finite at $x1 and $x2.

=item $rel_err

Maximum acceptable relative error (default: 10**-15, which is close to
the accuracy limits of double-precision floating point).  Estimates of
relative and absolute error are based on a comparison of the estimate
computed using C<2**n + 1> points with the estimate computed using
C<2**(n-1) + 1> points.

Once $min_split has been reached, computation stops as soon as
relative error drops below $rel_err, absolute error drops below
$abs_err, or $max_split is reached.

=item $abs_err

Maximum acceptable absolute error (default: 10**-40).

=item $max_split

At most C<2 ** $max_split + 1> different sample x values are used to
estimate the integral of C<f()>.

=item $min_split

At least C<2 ** $min_split + 1> different sample x values are used to
estimate the integral of C<f()>.

=item $Math::Integral::Romberg::return_point_count

This value defaults to 0.  If you set it to 1, then when invoked in a
list context, integral() will return a two-element list, containing
the estimate followed by the number of different x values used to
compute the estimate.

=item $Math::Integral::Romberg::abort

This value is set to 1 if neither the $rel_err nor the $abs_err
thresholds are reached before computation stops.  Once set, this
variable remains set until you reset it to 0.

=back

=head2 About the Algorithm

Romberg integration uses progressively higher-degree polynomial
approximations each time you double the number of sample points.  For
example, it uses a 2nd-degree polynomial approximation (as Simpson's
method does) after one split (2**1 + 1 sample points), and it uses a
10th-degree polynomial approximation after five splits (2**5 + 1
sample points).  Typically, this will greatly improve accuracy
(compared to simpler methods) for smooth functions, while not making
much difference for badly behaved ones.

=head1 AUTHOR

Eric Boesch (ebo@notes.dannet.dk)

=cut

sub integral {
  my $return_pts = wantarray && $Math::Integral::Romberg::return_point_count;
  my $abort = \$Math::Integral::Romberg::abort;
  my ($f,$lo,$hi,$rel_err,$abs_err,$max_split,$min_split)=@_;
  ($lo, $hi) = ($hi, $lo) if $lo > $hi;
  
  $rel_err	||= 1e-15;
  $abs_err	||= 1e-40;
  $max_split	||= 16;
  $min_split	||= 5;
  
  my ($estimate, $split, $steps);
  my $step_len = $hi - $lo;
  my $tot = (&$f($lo) + &$f($hi))/2;
  
  # tot is used to compute the trapezoid approximations.  It is more or
  # less a total of all f() values computed so far.  The trapezoid
  # method assigns half as much weight to f(hi) and f(lo) as it does to
  # all other f() values, so f(hi) and f(lo) are divided by two here.
  
  my @row = $estimate = $tot * $step_len; # 0th trapezoid approximation.
  
  for ($split = 1, $steps=2; ; $split++, $step_len /=2, $steps *= 2) {
    my ($x, $new_estimate);

    # Don't let $step_len drop below the limits of numeric precision.
    # (This should prevent infinite loops, but not loss of accuracy.)
    if ($lo + $step_len/$steps == $lo || $hi - $step_len/$steps == $hi) {
      $$abort = 1;
      return $return_pts ? ($estimate, $steps/2 + 1) : $estimate;
    }

    # Compute the (split)th trapezoid approximation.
    for ($x = $lo + $step_len/2; $x < $hi; $x += $step_len) {
      $tot += &$f($x);
    }
    unshift @row, $tot * $step_len / 2;
    
    # Compute the more refined approximations, based on the (split)th
    # trapezoid approximation and the various (split-1)th refined
    # approximations stored in @row.
    
    my $pow4 = 4;
    
    foreach my $td ( 1 .. $split ) {
      $row[$td] = $row[$td-1] +
	($row[$td-1]-$row[$td])/($pow4 - 1);
      $pow4 *= 4;
    }
    
    # row[0] now contains the (split)th trapezoid approximation,
    # row[1] now contains the (split)th Simpson approximation, and
    # so on up to row[split] which contains the (split)th Romberg
    # approximation.
    
    # Is this estimate accurate enough?
    $new_estimate = $row[-1];
    if (($split >= $min_split &&
	 (abs($new_estimate - $estimate) < $abs_err ||
	  abs($new_estimate - $estimate) < $rel_err * abs($estimate))) ||
	($split == $max_split && ($$abort = 1))) {
      return $return_pts ? ($new_estimate, $steps + 1) : $new_estimate;
    }
    $estimate = $new_estimate;
  }
}

1;
