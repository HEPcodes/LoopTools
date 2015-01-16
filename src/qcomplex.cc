class quad_complex {
public:
  quad_complex(Quad r = 0, Quad i = 0): re (r), im (i) { }
  quad_complex(double_complex x): re (x.real()), im (x.imag()) { }
  quad_complex& operator =  (const quad_complex&);
  quad_complex& operator =  (const double_complex&);
  quad_complex& operator += (const quad_complex&);
  quad_complex& operator += (const double_complex&);
  Quad real() const { return re; }
  Quad imag() const { return im; }
private:
  Quad re, im;
};

inline Quad imag(const quad_complex& x)
{
  return x.imag();
}

inline Quad real(const quad_complex& x)
{
  return x.real();
}

quad_complex& quad_complex::operator = (const quad_complex& y)
{
  re = y.real();
  im = y.imag();
  return *this;
}

quad_complex& quad_complex::operator = (const double_complex& y)
{
  re = y.real();
  im = y.imag();
  return *this;
}

quad_complex& quad_complex::operator += (const quad_complex& r)
{
  re += r.real();
  im += r.imag();
  return *this;
}

quad_complex& quad_complex::operator += (const double_complex& r)
{
  re += r.real();
  im += r.imag();
  return *this;
}

quad_complex operator + (const quad_complex& x, const quad_complex& y)
{
  return quad_complex(real(x) + real(y), imag(x) + imag(y));
}

quad_complex operator + (const quad_complex& x, const Quad& y)
{
  return quad_complex(real(x) + y, imag(x));
}

quad_complex operator + (const Quad& x, const quad_complex& y)
{
  return quad_complex(x + real(y), imag(y));
}

quad_complex operator - (const quad_complex& x, const quad_complex& y)
{
  return quad_complex(real(x) - real(y), imag(x) - imag(y));
}

quad_complex operator - (const quad_complex& x, const Quad& y)
{
  return quad_complex(real(x) - y, imag(x));
}

quad_complex operator * (const quad_complex& x, const quad_complex& y)
{
  return quad_complex(
    real(x)*real(y) - imag(x)*imag(y),
    real(x)*imag(y) + imag(x)*real(y) );
}

quad_complex operator * (const quad_complex& x, const Quad& y)
{
  return quad_complex(real(x)*y, imag(x)*y);
}

quad_complex operator * (const Quad& x, const quad_complex& y)
{
  return quad_complex(x*real(y), x*imag(y));
}

quad_complex operator * (const Quad& x, const double_complex& y)
{
  return quad_complex(real(y)*x, imag(y)*x);
}

quad_complex operator / (const quad_complex& x, const quad_complex& y)
{
  Quad ratio, den;

  if(fabs(real(y)) < fabs(imag(y))) {
    ratio = real(y)/imag(y);
    den = imag(y)*(1. + ratio*ratio);
    return quad_complex(
      (real(x)*ratio + imag(x))/den,
      (imag(x)*ratio - real(x))/den );
  }
  ratio = imag(y)/real(y);
  den = real(y)*(1. + ratio*ratio);
  return quad_complex(
    (real(x) + imag(x)*ratio)/den,
    (imag(x) - real(x)*ratio)/den );
}

quad_complex operator / (const quad_complex& x, const Quad& y)
{
  return quad_complex(real(x)/y, imag(x)/y);
}

quad_complex operator / (const Quad& x, const quad_complex& y)
{
  Quad ratio, den;

  if(fabs(real(y)) < fabs(imag(y))) {
    ratio = real(y)/imag(y);
    den = imag(y)*(1. + ratio*ratio);
    return quad_complex(x*ratio/den, -x/den);
  }
  ratio = imag(y)/real(y);
  den = real(y)*(1. + ratio*ratio);
  return quad_complex(x/den, -x*ratio/den);
}

quad_complex operator / (const Quad& x, const double_complex& y)
{
  return quad_complex(real(y)/x, imag(y)/x);
}

quad_complex operator - (const quad_complex& x)
{
  return quad_complex(-real(x), -imag(x));
}

quad_complex log(const quad_complex& x)
{
  return quad_complex(
    .5*log(real(x)*real(x) + imag(x)*imag(x)),
    atan2(imag(x), real(x)) );
}

double_complex q2d(const quad_complex& x)
{
  return double_complex(real(x).h(), imag(x).h());
}
