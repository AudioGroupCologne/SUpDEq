#define TYPEDEPARGS 0
#define SINGLEARGS
#define COMPLEXINDEPENDENT
#define OCTFILENAME comp_gabtight_long // change to filename
#define OCTFILEHELP "This function calls the C-library\n\
                     gd=comp_gabtight_long(g,a,M);\n\
                     Yeah."

#include "ltfat_oct_template_helper.h"
// octave_idx_type 32 or 64 bit signed integer


static inline void
fwd_gabtight_long(const Complex *g, const octave_idx_type L,
                  const octave_idx_type R, const octave_idx_type a,
                  const octave_idx_type M, Complex *gd)
{
    gabtight_long_cd(reinterpret_cast<const fftw_complex*>(g),
                     L, R, a, M,
                     reinterpret_cast<fftw_complex*>(gd));
}

static inline void
fwd_gabtight_long(const FloatComplex *g, const octave_idx_type L,
                  const octave_idx_type R, const octave_idx_type a,
                  const octave_idx_type M, FloatComplex *gd)
{
    gabtight_long_cs(reinterpret_cast<const fftwf_complex*>(g),
                     L, R, a, M,
                     reinterpret_cast<fftwf_complex*>(gd));
}

static inline void
fwd_gabtight_long(const double *g, const octave_idx_type L,
                  const octave_idx_type R, const octave_idx_type a,
                  const octave_idx_type M, double *gd)
{
    gabtight_long_d(g, L, R, a, M, gd);
}

static inline void
fwd_gabtight_long(const float *g, const octave_idx_type L,
                  const octave_idx_type R, const octave_idx_type a,
                  const octave_idx_type M, float *gd)
{
    gabtight_long_s(g, L, R, a, M, gd);
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list
octFunction(const octave_value_list& args, int nargout)
{
    MArray<LTFAT_TYPE> g = ltfatOctArray<LTFAT_TYPE>(args(0));
    const octave_idx_type L = g.rows();
    const octave_idx_type R = g.columns();
    const octave_idx_type a = args(1).int_value();
    const octave_idx_type M = args(2).int_value();

    MArray<LTFAT_TYPE> gd(dim_vector(L, R));

    fwd_gabtight_long(g.data(), L, R, a, M, gd.fortran_vec());

    return octave_value(gd);
}
