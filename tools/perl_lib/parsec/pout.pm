package parsec::pout;

use Carp;
use vars qw($AUTOLOAD);  # it's a package global
use parsec::threevec;
use parsec::vector;
use parsec::tensor;

sub new {
    my $that  = shift;
    my $class = ref($that) || $that;
    my $string = "parsec\.out";
    if (@_) {
        $string = shift;
    }
    my $self = \ $string;
    bless $self, $class;
    return $self;
}

sub energy {
    my $self = shift;
    my $type = ref($self) or croak "$self is not an object";
    my $energy = 0.0;
    open(INPUT, "$$self") or
        die "Can't read $$self $!";
    while (<INPUT>) {
        if (/Total Energy[ \t]+=[ \t]+([-\d\.]+)/) {
            $energy = $1;
        }
    }
    close(INPUT);
    return $energy;
}

sub energy_vals {
    my $self = shift;
    my $type = ref($self) or croak "$self is not an object";
    my %energy_hash = ();

    open(INPUT, "$$self") or
        die "Can't read $$self $!";

    while($_=<INPUT>) {
       if(/Eigenvalue Energy\s+\=\s+(-?\d*.\d+)/) 
         { $energy_hash{'eigenvalues_energy'}=$1; }
       if(/Hartree Energy\s+=\s+(-?\d*.\d+)/) 
         { $energy_hash{'hartree'}=$1; } 
       if(/Integral_\{Vxc\*rho\}\s+\=\s+(-?\d*.\d+)/) 
         { $energy_hash{'vxc_rho'}=$1; }
       if(/Integral\{eps_xc\*rho\}\s+\=\s+(-?\d*.\d+)/) 
         { $energy_hash{'e_xc'}=$1; }
       if(/Alpha Energy\s+\=\s+(-?\d*.\d+)/) 
         { $energy_hash{'alpha'}=$1; }
       if(/Ion-Ion Energy\s+\=\s+(-?\d*.\d+)/) 
         { $energy_hash{'ion_ion_energy'}=$1; }
       if(/Electron-Ion energy\s+\=\s+(-?\d*.\d+)/) 
         { $energy_hash{'electron_ion_energy'}=$1; }
        if (/Total Energy[ \t]+=[ \t]+([-\d\.]+)/) {
           $energy_hash{'total_energy'} = $1; }
    }
    close(INPUT);
   
    return %energy_hash;
}
       

sub forcevectors {     # Read the forces
    my $self = shift;
    my $type = ref($self) or croak "$self is not an object";

    open(INPUT, "$$self") or
        die "Can't read $$self $!";
    while (<INPUT>) {
        if (/coordinates and total forces/) {
            <INPUT>; # Need to skip a line.
            my @forces =();
            while (($_ = <INPUT>) !~ /Time/) {
            if(/\d+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)/) {
                push @forces, new parsec::vector(
                     /\d+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)/); }
            }
            close(INPUT);
            return @forces; # The force vectors
        }
    }
    close(INPUT);
    return "";
}

sub eigenvalues {     # Read the forces
    my $self = shift;
    my $type = ref($self) or croak "$self is not an object";
    my $x,$y,$z,$i,$flag;
    my $k_pat_s,$k_pat_e;
    my $pat_s,$pat_e;
    my $gamma_flag;
    my $so_flag;
    my $spin_polarized_flag;


    my @eigenvals_tmp = ();

#    $k_pat_s="k-point\s+\=\s+\S+\s+\S+\s+\S+\s+\(a\.u";
    $pat_s="Fermi level at";
    $k_pat_e="k-point\\s\+\\\=\\s\+\\(\\S\+)\\s\+(\\S\+)\\s\+(\\S\+)\\s\+\\(recip";
    $gamma_flag=1;
    $so_flag=0;
    $spin_polarized_flag=0;

    open(INPUT, "$$self") or
        die "Can't read $$self $!";
    while ($_=<INPUT>) {
        if (/Number of states:\s+(\d+)/) 
        { $band_num=$1; }
        if (/kpoints selected/)
        { $gamma_flag=0; }
        if (/Spin-Orbit corrected computation!/)
        { $so_flag=1; }
        if (/Spin-polarized computation!/)
        { $spin_polarized_flag=1; }
        if (/$pat_s/) {
              $i=1;
              $flag=0;
            while (($_ = <INPUT>) !~ /Hartree potential time/ 
                   & $flag==0) {
              if(/k-point\s+\=\s+(\S+)\s+(\S+)\s+(\S+)\s+\(recip/ | $gamma_flag) {
                if($gamma_flag) {
                $x=0;
                $y=0;
                $z=0; }
                else {
                 $x=$1;
                 $y=$2;
                 $z=$3; }

		<INPUT>;
                <INPUT>;
                while (($_ = <INPUT>) !~ /\(a\.u\./
                        & $_ !~/Hartree potential time/) {
                  if(/(\d+)\s+\S+\s+(\S+)\s+(\S+)/) {
                     push @eigenvals_tmp, [$i,$x,$y,$z,$1,$2,$3];
#                     print "$flag,$i,$x,$y,$z,$1,$2,$3\n";
                     } # end of if
                  } # end of this k-point
                  if(/Hartree potential time/) 
                  { $flag=1; }
                  $i=$i+1;
              } # end of if k-point
            } # end of this round of k-points
          } # end of if start of k-point
        } # end of general loop
    close(INPUT);
    $kpt_num=$i-1;
    $offs=$#eigenvals_tmp-$kpt_num*$band_num*($spin_polarized_flag+1)+1;
    my @eigenvals = splice(@eigenvals_tmp,$offs);
    return @eigenvals;
}

1;
