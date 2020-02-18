package font_rotate;

use strict;
use warnings;
use SVG;


package GD::SVG::Image;


sub stringrotate {

    my ($self, $x, $y, $font_family, $font_style, $font_size, $text, $color_index, $rotate) = @_;
    my $img = $self->currentGroup;
    my $id = $self->_create_id($x,$y);
    my $color = $self->_get_color($color_index);
    my $result =$img->text(
	                     id=>$id,
                             'font-family'  => $font_family,
                             'font-weight'  => $font_style,        #'font-style' => $font_style,
                             'font-style'   => $font_style,
                             'font-size'    => $font_size,
	                     'fill'         => $color,
                             'fill-opacity' => 1,           # $opacity
#                            'stroke'       => 'rgb(250,123,23)',
#	                     'writing-mode' => 'tb',
	                     'transform' => "translate($x,$y) rotate($rotate)",
	                  )->cdata($text);
   return $result;
}


1;

__END__  
