
package parsec::vector;

use strict;
use Carp;
use vars qw($AUTOLOAD);  # it's a package global

use overload
     'neg' => '_negate',
#       '~' => '_transpose',
#    'bool' => '_boolean',
#       '!' => '_not_boolean',
      '""' => '_stringify',
     'abs' => '_abs', # Does this work? I don't know
       '+' => '_add',
       '-' => '_subtract',
       '*' => '_multiply', # Does either dot or scalar multiply
       '/' => '_divide',
       'x' => '_cross',
#      '+=' => '_assign_add',
#      '-=' => '_assign_subtract',
#      '*=' => '_assign_multiply',
#      '==' => '_equal',
#      '!=' => '_not_equal',
#       '<' => '_less_than',
#      '<=' => '_less_than_or_equal',
#       '>' => '_greater_than',
#      '>=' => '_greater_than_or_equal',
#      'eq' => '_equal',
#      'ne' => '_not_equal',
#      'lt' => '_less_than',
#      'le' => '_less_than_or_equal',
#      'gt' => '_greater_than',
#      'ge' => '_greater_than_or_equal',
#       '=' => '_clone',
'fallback' =>   undef;

sub new {
    my $that  = shift;
    my $class = ref($that) || $that;
    my $self;
    if (@_) { # We have an argument
        if ($#_ == 2) { # it is an array of 3 nums
            my @vec = @_;
            $self = \ @vec;
        } elsif ((defined $_[0]) && ref($_[0]) &&
                   (ref($_[0]) =~ /^Math::MatrixReal$/)) {
            # Read from a matrix
            my ($rows, $cols) = $_[0]->dim;
            my $string;
            my $mat = $_[0];
            if ($rows == 3) {
                $string = "$mat";
                my @arr = $string =~ 
                    /\s*\[\s*(\S+)\s*\]\s*\[\s*(\S+)\s*\]\s*\[\s*(\S+)\s*\]/;
                $self = \ @arr;
            } elsif ($cols == 3) {
                $string = "$mat";
                my @arr = $string =~ 
                    /\s*\[\s*(\S+)\s+(\S+)\s+(\S+)\s*\]/;
                $self = \ @arr;
            } else {
                print "Matrix has wrong dimension for vector: $mat\n";
                exit;
            }
        } else {
            print "Bad argument to vector constructor @_\n";
        }
    } else {
        $self = \ (0, 0, 0);
    }
    bless $self, $class;
    return $self;
}

sub _stringify {
    my($object,$argument,$flag) = @_;
    my($j,$s);

    $s = '';
    $s .= "";
    for ( $j = 0; $j < 3; $j++ ) {
        $s .= sprintf("% #-19.16g ", $object->[$j]);
    }
    $s .= "";
    return($s);
}

sub _negate {
    my($object) = @_;
    return new parsec::vector(-$object->[0], -$object->[1], -$object->[2]);
}

sub _add {
    my($object,$argument) = @_;
    return new parsec::vector($object->[0]+$argument->[0], 
                                 $object->[1]+$argument->[1],
                                 $object->[2]+$argument->[2]);
}

sub _subtract {
    my($object,$argument) = @_;
    return new parsec::vector($object->[0]-$argument->[0], 
                                 $object->[1]-$argument->[1],
                                 $object->[2]-$argument->[2]);
}



sub _multiply
{
    my($object,$argument,$flag) = @_;
    my($name) = "'*'"; #&_trace($name,$object,$argument,$flag);
    my($temp);

    if ((defined $argument) && ref($argument) =~ /^parsec::vector/) {
        # we are multiplying two vectors
        return( _dot($argument,$object) );
    } elsif ((defined $argument) && !(ref($argument))) {
        # we do a scalar multiply
        return new parsec::vector($object->[0]*$argument, 
                                     $object->[1]*$argument,  
                                     $object->[2]*$argument);
    } else { # ERROR!!!
        if (defined $argument) {
            my $thetype = ref($argument);
            croak "parsec::vector $name: $thetype wrong argument type";
        } else {
            croak "parsec::vector $name: undefined argument";
        }
    }
}

sub _cross {
    my($object,$argument) = @_;
    return new parsec::vector($object->[1]*$argument->[2]
                                 -$object->[2]*$argument->[1], 
                                 $object->[2]*$argument->[0]
                                 -$object->[0]*$argument->[2],  
                                 $object->[0]*$argument->[1]
                                 -$object->[1]*$argument->[0]);
}

sub _dot {
    my($object,$argument) = @_;
    return ($object->[0]*$argument->[0]+ 
            $object->[1]*$argument->[1]+
            $object->[2]*$argument->[2]);
}

sub _divide {
    my($object,$argument) = @_;
    if (ref($argument) =~ /^parsec::tensor$/) {
# Divide by a tensor!!!  :)
        my $volume = $argument->det;
        $argument = ~$argument;
        return new parsec::vector(
                   ($argument->[1]x$argument->[2])*$object/$volume,
                   ($argument->[2]x$argument->[0])*$object/$volume,
                   ($argument->[0]x$argument->[1])*$object/$volume);
    } else {
        return new parsec::vector($object->[0]/$argument, 
                                   $object->[1]/$argument, 
                                   $object->[2]/$argument);
    }
}

sub _abs {
    my($object) = @_;
    return sqrt($object->[0]*$object->[0]+ 
                $object->[1]*$object->[1]+
                $object->[2]*$object->[2]);
}
1;


__END__

=head1 NAME

parsec::vector - three-vectors for use with parsec

Implements three vectors

=head1 DESCRIPTION

Creates a three vector object, which is actually an array of numbers.
Vector is intended to make threevec obsolete, as I prefer this name, and
the parsec:: avoids any name collisions with other vector types.

=head1 SYNOPSIS

=over 2

=item *

C<use parsec::vector;>

Makes the methods and overloaded operators of this module available
to your program.

=item *

C<$vector = new parsec::vector;>

Create a new vector.

=item *

C<$vector = new parsec::vector(1,2,3);>

Create a vector with value (1,2,3).

=head1 ARITHMETIC

All the ordinary arithmetic operators work as you would expect of vectors.
A couple of trickier things are possible, which I will mention here.

=item *

C<$scalar = $vector1 * $vector2>

The dot product of two vectors.

=item *

C<$scalar = $vector1 x $vector2>

The cross product of two vectors (be careful of precedence).

=item *

C<$vector1 = $tensor * $vector2>

A tensor can be multiplied by a vector.

=item *

C<$vector1 = $vector2 /  $tensor>

A vector may be divided by a tensor.  This computes the inverse of the
tensor and then operates it on the vector.  This can be handy for going
between coordinate relative to the lattice vectors and cartesian
coordinates.

=head1 FUNCTIONS

=item *

C<abs($vector)>

Returns the magnitude the vector.

=back

=head1 SEE ALSO

parsec::tensor.

=head1 VERSION

This man page documents parsec::vector version 0.1.

=head1 AUTHOR

David Roundy <droundy@physics.berkeley.edu>.

=head1 CREDITS

Thanks to Gilles Santi, for giving me the idea of writing this.

=head1 COPYRIGHT

Copyright (c) 2000, David Roundy. All rights reserved.

=head1 LICENSE AGREEMENT

This package is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.




