package parsec::tensor;

use strict;
use Carp;
use vars qw($AUTOLOAD);  # it's a package global
use parsec::vector;     # will need this

use overload
     'neg' => '_negate',
       '~' => '_transpose',
#    'bool' => '_boolean',
#       '!' => '_not_boolean',
      '""' => '_stringify',
     'det' => 'det',
#     'abs' => '_norm',
       '+' => '_add',
       '-' => '_subtract',
       '*' => '_multiply',
       '/' => '_divide',
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
        if ($#_ == 2) { # it is an array of 3 vectors
            my @vec = @_;
            if (ref($vec[0]) =~ /^parsec::vector$|^parsec::threevec$/) {
                $self = \ @vec;
            } else { # assume we have 3 numbers, and put them on the diagonal.
                my @temp;
                $temp[0] = new parsec::vector($vec[0], 0, 0);
                $temp[1] = new parsec::vector(0, $vec[1], 0);
                $temp[2] = new parsec::vector(0, 0, $vec[2]);
                $self = \ @temp;
            }
        } elsif ($#_ == 1) { # it is 2 vectors--tensor product
            my @vec = @_;
            ref($vec[0]) =~ /^parsec::vector$|^parsec::threevec$/ or
                croak "parsec::tensor->new requires vectors\n";
            $self = [$vec[1]*$vec[0]->[0],
                       $vec[1]*$vec[0]->[1],
                       $vec[1]*$vec[0]->[2]];
        } elsif ((defined $_[0]) && ref($_[0]) &&
                   (ref($_[0]) =~ /^parsec::vector$|^parsec::threevec$/)) {
            # Read from one vector (along diagonal)
            my $vec = $_[0];
            my @temp;
            $temp[0] = new parsec::vector($vec->[0], 0, 0);
            $temp[1] = new parsec::vector(0, $vec->[1], 0);
            $temp[2] = new parsec::vector(0, 0, $vec->[2]);
            $self = \ @temp;
        } elsif ((defined $_[0]) && ref($_[0]) &&
                   (ref($_[0]) =~ /^Math::MatrixReal$/)) {
            # Read from a matrix
            my ($rows, $cols) = $_[0]->dim;
            my $mat = $_[0];
            croak "I haven't yet enabled creating from Math::MatrixReal\n";
            if ($rows == 3 && $cols == 3 ) {
                1; # do nothing
            } else {
                print "Matrix has wrong dimension for tensor: $mat\n";
                exit;
            }
        } else {
            print "Bad argument to tensor constructor @_\n";
        }
    } else {
        $self = \ (new parsec::vector(1,0,0), 
                   new parsec::vector(0,1,0), 
                   new parsec::vector(0,0,1));
    }
    bless $self, $class;
    return $self;
}

sub _stringify {
    my($object,$argument,$flag) = @_;

    return "$object->[0]\n$object->[1]\n$object->[2]\n";
}

sub _negate {
    my($object) = @_;
    return new parsec::tensor(-$object->[0], -$object->[1], -$object->[2]);
}

sub _add {
    my($object,$argument) = @_;
    if (ref($argument)=~/parsec::tensor/) {
        return new parsec::tensor($object->[0]+$argument->[0], 
                                   $object->[1]+$argument->[1],
                                   $object->[2]+$argument->[2]);
    } else {
        my $thetype = ref($argument);
        croak "You can't add a tensor to $thetype";
    }
}

sub _subtract {
    my($object,$argument) = @_;
    return new parsec::tensor($object->[0]-$argument->[0], 
                               $object->[1]-$argument->[1],
                               $object->[2]-$argument->[2]);
}

sub _multiply {
    my($object,$argument,$flag) = @_;
    my($name) = "'*'"; #&_trace($name,$object,$argument,$flag);
    my($temp);

    if ((defined $argument) && 
        ref($argument) =~ /^parsec::threevec$|^parsec::vector$/) {
        # we are multiplying a tensor and a vector.  Easy! :)
        return new parsec::vector($object->[0]*$argument, 
                                   $object->[1]*$argument,  
                                   $object->[2]*$argument);
    } elsif ((defined $argument) && 
        ref($argument) =~ /^parsec::tensor$/) {
        # we are multiplying a tensor and another tensor
        return( _mul_tensor($argument,$object) );
    } elsif ((defined $argument) && !(ref($argument))) {
        # we do a scalar multiply
        return new parsec::tensor($object->[0]*$argument, 
                                   $object->[1]*$argument,  
                                   $object->[2]*$argument);
    } else {
        croak "parsec::tensor $name: wrong argument type";
    }
}

sub _mul_tensor {
    my($object,$argument) = @_;
    $argument = ~$argument;
    return new parsec::tensor(
       new parsec::vector($object->[0]*$argument->[0],
                           $object->[0]*$argument->[1],
                           $object->[0]*$argument->[2]),
       new parsec::vector($object->[1]*$argument->[0],
                           $object->[1]*$argument->[1],
                           $object->[1]*$argument->[2]),
       new parsec::vector($object->[2]*$argument->[0],
                           $object->[2]*$argument->[1],
                           $object->[2]*$argument->[2]));
}

sub _divide {
    my($object,$argument, $reversed) = @_;
    if (!$reversed) { # Ordinary scalar division (we hope).
        if (!ref($argument)) {
            return new parsec::tensor($object->[0]/$argument, 
                                       $object->[1]/$argument, 
                                       $object->[2]/$argument);
        } elsif (ref($argument)=~/^parsec::tensor/) {
            return return (1/$argument)*$object;
        } else {
            my $thetype = ref($argument);
            croak "You can't divide a tensor by a $thetype";
        }
    } else { 
# We're dividing something BY a matrix (i.e. multiply by the matrix inverse).
        if (!ref($argument)) {
            my $volume = det($object);
            if ($volume == 0.0) {
                croak "I can't invert a singular tensor";
            }
            return ~new parsec::tensor(
                   ($object->[1]x$object->[2])*$argument/$volume,
                   ($object->[2]x$object->[0])*$argument/$volume,
                   ($object->[0]x$object->[1])*$argument/$volume);
        } else {
            my $thetype = ref($argument);
            croak "Can't divide $thetype by tensor";
        }
    }
}

sub _transpose {
    my($object) = @_;
    return new parsec::tensor(
       new parsec::vector($object->[0]->[0],
                           $object->[1]->[0],
                           $object->[2]->[0]),
       new parsec::vector($object->[0]->[1],
                           $object->[1]->[1],
                           $object->[2]->[1]),
       new parsec::vector($object->[0]->[2],
                           $object->[1]->[2],
                           $object->[2]->[2]));
}

sub det {
    my($object) = @_;
    return ($object->[0]x$object->[1])*$object->[2];
}
1;


__END__

=head1 NAME

parsec::tensor - simple tensors for use with parsec

Implements tensors

=head1 DESCRIPTION

Creates a tensor object.  Overloads arithmetic operators so you can easily
do arithmetic, and also allows you to transpose and calculate the inverse
easily.  It is heavily unoptimized, but how long can it take to deal with a
3x3 matrix?

=head1 SYNOPSIS

=item *

C<use parsec::tensor;>

Makes the methods and overloaded operators of this module available
to your program.

=item *

C<$input = new parsec::tensor;>

Create a new identity tensor.

=item *

C<$input = new parsec::tensor(1,2,3);>

Create a diagonal tensor with elements C<1>, C<2> and C<3> along the
diagonal.

=item *

C<$input = new parsec::tensor($vector);>

Create a diagonal tensor with C<$vector> elements along the diagonal.

=item *

C<$input = new parsec::tensor($vector1, $vector2);>

Create a tensor having the value of the tensor product of C<$vector1>
and C<$vector2>.

=item *

C<$input = new parsec::tensor($vector1, $vector2, $vector3);>

Create a tensor with row vectors corresponding to C<$vector1>, C<$vector2>
and C<$vector3>.

=head2 FUNCTIONS

=item *

C<$volume = det($tensor);>

Calculate the determinant of a given tensor.

=head2 ARITHMETIC

All the ordinary arithmetic operators work as you would expect on tensors
and vectors.

=item *

C<$inverse = 1/$tensor;>

Calculate the inverse of a given tensor.

=item *

C<$newvector = $oldvector/$tensor;>

Multiply a vector by the inverse of a tensor.  This works out such that the
following is true:

  $vector == $tensor*$vector/$tensor;

=item *

C<$tensortensor = $oldtensor/$tensor;>

Multiply a tensor by the inverse of a tensor.  One could either multiply on
the left or on the right.  I multiply with the inverse on the left, such
that the following is true:

  $A == $B*$A/$B;

=back

=head1 SEE ALSO

parsec::vector.

=head1 VERSION

This man page documents parsec::tensor version 0.1.

=head1 AUTHOR

David Roundy <droundy@physics.berkeley.edu>.

=head1 CREDITS

Thanks to Gilles Santi, for giving me the idea of writing this.

=head1 COPYRIGHT

Copyright (c) 2000, David Roundy. All rights reserved.

=head1 LICENSE AGREEMENT

This package is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.




