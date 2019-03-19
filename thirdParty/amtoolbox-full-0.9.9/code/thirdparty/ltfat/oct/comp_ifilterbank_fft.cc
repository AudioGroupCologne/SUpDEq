#define TYPEDEPARGS 0, 1
#define SINGLEARGS
#define COMPLEXARGS
#define OCTFILENAME comp_ifilterbank_fft // change to filename
#define OCTFILEHELP "This function calls the C-library \n\
                     F = comp_ifilterbank_fft(c,G,a) \n\
                     Yeah."


#include "ltfat_oct_template_helper.h"
// octave_idx_type 32 or 64 bit signed integer


static inline void
fwd_ifilterbank_fft(const Complex *c[], const Complex *filt[],
                    const ltfatInt L, const ltfatInt W,
                    const ltfatInt a[], const ltfatInt M,
                    Complex *f)
{
    ifilterbank_fft_d(reinterpret_cast<const fftw_complex **>(c),
                      reinterpret_cast<const fftw_complex **>(filt),
                      L, W, a, M,
                      reinterpret_cast<fftw_complex *>(f));
}

static inline void
fwd_ifilterbank_fft(const FloatComplex *c[], const FloatComplex *filt[],
                    const ltfatInt L, const ltfatInt W,
                    const ltfatInt a[], const ltfatInt M,
                    FloatComplex *f)
{
    ifilterbank_fft_s(reinterpret_cast<const fftwf_complex **>(c),
                      reinterpret_cast<const fftwf_complex **>(filt),
                      L, W, a, M,
                      reinterpret_cast<fftwf_complex *>(f));
}

template <class LTFAT_TYPE, class LTFAT_REAL, class LTFAT_COMPLEX>
octave_value_list octFunction(const octave_value_list& args, int nargout)
{
    // Input data
    Cell c = args(0).cell_value();
    // Cell aray containing impulse responses
    Cell G = args(1).cell_value();
    // Subsampling factors
    Matrix aDouble = args(2).matrix_value();

    // Output length
    const octave_idx_type L  = G.elem(0).rows();
    // Number of channels
    const octave_idx_type W  = c.elem(0).columns();
    // Number of filters
    const octave_idx_type M = G.numel();

    // Output signal
    MArray<LTFAT_COMPLEX> F(dim_vector(L, W));


    // Allocating temporary arrays
    // Output subband lengths
    // Impulse responses pointers
    OCTAVE_LOCAL_BUFFER (const LTFAT_COMPLEX*, GPtrs, M);
    // Output subbands pointers
    OCTAVE_LOCAL_BUFFER (const LTFAT_COMPLEX*, cPtrs, M);
    // Output cell elements array,
    OCTAVE_LOCAL_BUFFER (MArray<LTFAT_COMPLEX>, cElems, M);
    //
    OCTAVE_LOCAL_BUFFER (MArray<LTFAT_COMPLEX>, GElems, M);
    OCTAVE_LOCAL_BUFFER (ltfatInt, a, M);

    for (octave_idx_type m = 0; m < M; m++)
    {
        a[m] = (ltfatInt) aDouble(m);
        GElems[m] = ltfatOctArray<LTFAT_COMPLEX>(G.elem(m));
        GPtrs[m] = GElems[m].data();
        cElems[m] = ltfatOctArray<LTFAT_COMPLEX>(c.elem(m));
        cPtrs[m] = cElems[m].data();
    }


    fwd_ifilterbank_fft(cPtrs, GPtrs, L, W, a, M, F.fortran_vec());

    return octave_value(F);
}
