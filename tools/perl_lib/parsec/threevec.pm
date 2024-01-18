package parsec::threevec;

use Carp;
use vars qw($AUTOLOAD);  # it's a package global

#@ISA = ("Math::MatrixReal");

use overload
     'neg' => '_negate',
#       '~' => '_transpose',
#    'bool' => '_boolean',
#       '!' => '_not_boolean',
      '""' => '_stringify',
     'abs' => '_norm',
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
                print "Matrix has wrong dimension for threevec: $mat\n";
                exit;
            }
        } else {
            print "Bad argument to threevec constructor @_\n";
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
    $s .= " ";
    for ( $j = 0; $j < 3; $j++ ) {
        $s .= sprintf("% #-19.16g ", $object->[$j]);
    }
    $s .= "";
    return($s);
}

sub _negate {
    my($object) = @_;
    return new parsec::threevec(-$object->[0], -$object->[1], -$object->[2]);
}

sub _add {
    my($object,$argument) = @_;
    return new parsec::threevec($object->[0]+$argument->[0], 
                                 $object->[1]+$argument->[1],
                                 $object->[2]+$argument->[2]);
}

sub _subtract {
    my($object,$argument) = @_;
    return new parsec::threevec($object->[0]-$argument->[0], 
                                 $object->[1]-$argument->[1],
                                 $object->[2]-$argument->[2]);
}



sub _multiply
{
    my($object,$argument,$flag) = @_;
    my($name) = "'*'"; #&_trace($name,$object,$argument,$flag);
    my($temp);

    if ((defined $argument) && ref($argument) =~ /^parsec::threevec/) {
        # we are multiplying two threevecs
        return( _dot($argument,$object) );
    } elsif ((defined $argument) && !(ref($argument))) {
        # we do a scalar multiply
        return new parsec::threevec($object->[0]*$argument, 
                                     $object->[1]*$argument,  
                                     $object->[2]*$argument);
    } else {
        croak "parsec::threevec $name: wrong argument type";
    }
}

sub _cross {
    my($object,$argument) = @_;
    return new parsec::threevec($object->[1]*$argument->[2]
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
    return new parsec::threevec($object->[0]/$argument, 
                                 $object->[1]/$argument, 
                                 $object->[2]/$argument);
}

sub _norm {
    my($object) = @_;
    return sqrt($object->[0]*$object->[0]+ 
                $object->[1]*$object->[1]+
                $object->[2]*$object->[2]);
}
1;


__END__

=head1 NAME

parsec::threevec - three-vectors for use with parsec

Implements three vectors

=head1 DESCRIPTION

Creates a three vector object, which is actually an array of numbers.

=head1 SYNOPSIS

=over 2

=item *

C<use parsec::threevec;>

Makes the methods and overloaded operators of this module available
to your program.

=item *

C<$input = new parsec::threevec;>

threevec constructor method.

=head1 SEE ALSO

parsec::input, parsec::OUT, parsec::PW_LOG.

=head1 VERSION

This man page documents parsec::threevec version 0.0.

=head1 AUTHOR

David Roundy <droundy@physics.berkeley.edu>.

=head1 CREDITS

Thanks to Gilles Santi, for giving me the idea of writing this.

=head1 COPYRIGHT

Copyright (c) 2000, David Roundy. All rights reserved.

=head1 LICENSE AGREEMENT

This package is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.




