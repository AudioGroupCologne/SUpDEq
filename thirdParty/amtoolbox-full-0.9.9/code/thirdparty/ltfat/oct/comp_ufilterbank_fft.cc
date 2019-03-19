#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define OCTFILENAME comp_ufilterbank_fft // change to filename
#define OCTFILEHELP "This function calls the C-library\n\
                     Usage: c=comp_ufilterbank_fft(f,g,a);\n\
                     Yeah."

#include "ltfat_oct_template_helper.h"

static inline void
fwd_ufilterbank_fft(const Complex *f, const Complex *g,
                    const octave_idx_type L, const octave_idx_type gl,
                    const octave_idx_type W, const octave_idx_type a,
                    const octave_idx_type M, Complex *cout)
{
    ufilterbank_fft_d(reinterpret_cast<const fftw_complex*>(f),
                      reinterpret_cast<const fftw_complex*>(g),
                      L, gl, W, a, M,
                      reinterpret_cast<fftw_complex*>(cout));
}

static inline void
fwd_ufilterbank_fft(const FloatComplex *f, const FloatComplex *g,
                    const octave_idx_type L, const octave_idx_type gl,
                    const octave_idx_type W, const octave_idx_type a,
                    const octave_idx_type M, FloatComplex *cout)
{
    ufilterbank_fft_s(reinterpret_cast<const fftwf_complex*>(f),
                      reinterpret_cast<const fftwf_complex*>(g),
                      L, gl, W, a, M,
                      reinterpret_cast<fftwf_complex*>(cout));
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
    DEBUGINFO;

    MArray<LTFAT_TYPE> f = ltfatOctArray<LTFAT_TYPE>(args(0));
    MArray<LTFAT_TYPE> g = ltfatOctArray<LTFAT_TYPE>(args(1));
    const octave_idx_type a = args(2).int_value();

    const octave_idx_type L  = f.rows();
    const octave_idx_type W  = f.columns();
    const octave_idx_type gl = g.rows();
    const octave_idx_type M  = g.columns();
    const octave_idx_type N = L / a;

    dim_vector dims_out(N, M, W);
    dims_out.chop_trailing_singletons();
    MArray<LTFAT_COMPLEX> cout(dims_out);

    fwd_ufilterbank_fft(f.data(), g.data(), L, gl, W, a, M, cout.fortran_vec());

    return octave_value(cout);
}
