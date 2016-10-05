/************************************************************************/
/*                                                                      */
/*               Copyright 2014-2015 by Ullrich Koethe                  */
/*                                                                      */
/*    This file is part of the VIGRA2 computer vision library.          */
/*    The VIGRA2 Website is                                             */
/*        http://ukoethe.github.io/vigra2                               */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

#ifndef VIGRA2_IMAGEIO_IMPEXALPHA_HXX
#define VIGRA2_IMAGEIO_IMPEXALPHA_HXX

#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <vigra2/shape.hxx>
#include <vigra2/sized_int.hxx>

#include "imageinfo.hxx"
#include "impex.hxx"
#include "impexbase.hxx"

namespace vigra
{
/** \addtogroup VigraImpex
 * @{
*/
namespace detail
{

template <class ValueType, int N, class T>
void
    read_image_bands_and_alpha(Decoder* decoder,
        PointerND<N, T> p, ArrayIndex channels,
        PointerND<N, T> alpha)
{
    vigra_assert(p.ndim() == 3, "Number of dimensions must be 3.");
    vigra_assert(alpha.ndim() == 3, "Number of dimensions must be 3.");

    const ArrayIndex width(decoder->getWidth());
    const ArrayIndex height(decoder->getHeight());
    const ArrayIndex bands(decoder->getNumBands() - decoder->getNumExtraBands());
    const ArrayIndex offset(decoder->getOffset());

    vigra_precondition(decoder->getNumExtraBands() == 1,
        "vigra::detail::read_image_bands_and_alpha: expecting exactly one alpha band");
    vigra_precondition(bands == channels || bands == 1,
        "vigra::detail::read_image_bands_and_alpha: number of channels does not match");

    // OPTIMIZATION: Specialization for the most common cases
    // (1 or 3 channels)
    if (channels == 1)
    {
        for (ArrayIndex y = 0; y != height; ++y, p.inc(0), alpha.inc(0))
        {
            decoder->nextScanline();

            auto scanline0 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(0));
            auto scanline1 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(1));

            for (ArrayIndex x = 0; x < width; ++x, p.inc(1), alpha.inc(1))
            {
                *p = *scanline0;
                *alpha = *scanline1;

                scanline0 += offset;
                scanline1 += offset;
            }

            p.move(1, -width);
            alpha.move(1, -width);
        }
    }
    else if (channels == 3)
    {
        const ValueType* scanline_0;
        const ValueType* scanline_1;
        const ValueType* scanline_2;
        const ValueType* scanline_3;

        for (ArrayIndex y = 0; y != height; ++y, p.inc(0), alpha.inc(0))
        {
            decoder->nextScanline();

            scanline_0 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(0));

            if (bands == 1)
            {
                scanline_1 = scanline_0;
                scanline_2 = scanline_0;
                scanline_3 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(1));
            }
            else
            {
                scanline_1 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(1));
                scanline_2 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(2));
                scanline_3 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(3));
            }

            for (ArrayIndex x = 0; x < width; ++x, p.inc(1), alpha.inc(1))
            {
                *p = *scanline_0;
                p.inc(2);
                *p = *scanline_1;
                p.inc(2);
                *p = *scanline_2;
                p.move(2, -2);

                *alpha = *scanline_3;

                scanline_0 += offset;
                scanline_1 += offset;
                scanline_2 += offset;
                scanline_3 += offset;
            }

            p.move(1, -width);
            alpha.move(1, -width);
        }
    }
    else
    {
        std::vector<const ValueType*> scanlines(channels + 1);

        for (ArrayIndex y = 0; y != height; ++y, p.inc(0), alpha.inc(0))
        {
            decoder->nextScanline();

            scanlines[0] = static_cast<const ValueType*>(decoder->currentScanlineOfBand(0));

            if (bands == 1)
            {
                for (ArrayIndex i = 1; i != channels; ++i)
                {
                    scanlines[i] = scanlines[0];
                }
                scanlines[channels] = static_cast<const ValueType*>(decoder->currentScanlineOfBand(1));
            }
            else
            {
                for (ArrayIndex i = 1; i != channels; ++i)
                {
                    scanlines[i] = static_cast<const ValueType*>(decoder->currentScanlineOfBand(i));
                }
                scanlines[channels] = static_cast<const ValueType*>(decoder->currentScanlineOfBand(channels));
            }

            for (ArrayIndex x = 0; x < width; ++x, p.inc(1), alpha.inc(1))
            {
                for (ArrayIndex i = 0; i != channels; ++i, p.inc(2))
                {
                    *p = *scanlines[i];
                    scanlines[i] += offset;
                }
                p.move(2, -channels);

                *alpha = *scanlines[channels];
                scanlines[channels] += offset;
            }

            p.move(1, -width);
            alpha.move(1, -width);
        }
    }
}

template <class T1, class T2, int N>
void
    importImageAlpha(const ImageImportInfo& import_info,
        ArrayViewND<N, T1> view, ArrayViewND<N, T2> alpha)
{
    std::unique_ptr<Decoder> decoder(vigra::decoder(import_info));

    vigra_precondition(decoder->getNumExtraBands() == 1,
        "importImageAlpha(): expecting exactly one alpha channel in file.");

    const ArrayIndex bands = decoder->getNumBands() - 1;

    // FIXME: ensure C order.
    auto ca = view.ensureChannelAxis(2);
    vigra_precondition(ca.ndim() == 3,
        "importImageAlpha(): wrong number of axes in target array.");
    vigra_precondition(ca.shape(0) == decoder->getHeight() && ca.shape(1) == decoder->getWidth(),
        "importImageAlpha(): shape mismatch between input and output.");
    vigra_precondition(ca.shape(2) == bands || bands == 1,
        "importImageAlpha(): Number of channels in input and destination image don't match.");
    auto p = ca.pointer_nd();
    ArrayIndex channels = ca.shape(2);

    auto aa = alpha.ensureChannelAxis(2);
    vigra_precondition(aa.ndim() == 3,
        "importImageAlpha(): wrong number of axes in alpha array.");
    vigra_precondition(aa.shape(0) == decoder->getHeight() && aa.shape(1) == decoder->getWidth(),
        "importImageAlpha(): shape mismatch between data and alpha array.");
    vigra_precondition(aa.shape(2) == 1,
        "importImageAlpha(): alpha array must have exactly one channel.");
    auto a = aa.pointer_nd();

    switch (pixel_t_of_string(decoder->getPixelType()))
    {
    case UNSIGNED_INT_8:
        read_image_bands_and_alpha<uint8_t>(decoder.get(), p, channels, a);
        break;
    case UNSIGNED_INT_16:
        read_image_bands_and_alpha<uint16_t>(decoder.get(), p, channels, a);
        break;
    case UNSIGNED_INT_32:
        read_image_bands_and_alpha<uint32_t>(decoder.get(), p, channels, a);
        break;
    case SIGNED_INT_16:
        read_image_bands_and_alpha<int16_t>(decoder.get(), p, channels, a);
        break;
    case SIGNED_INT_32:
        read_image_bands_and_alpha<int32_t>(decoder.get(), p, channels, a);
        break;
    case IEEE_FLOAT_32:
        read_image_bands_and_alpha<float>(decoder.get(), p, channels, a);
        break;
    case IEEE_FLOAT_64:
        read_image_bands_and_alpha<double>(decoder.get(), p, channels, a);
        break;
    default:
        vigra_fail("detail::importImageAlpha(): should never be reached.");
    }

    decoder->close();
}

template<class ValueType, class T, int N, class ImageScaler>
void
    write_image_bands_and_alpha(Encoder* encoder, PointerND<N, T> p, PointerND<N, T> alpha,
        const Shape<N> &shape,
        const ImageScaler& image_scaler, const ImageScaler& alpha_scaler)
{
    typedef RequiresExplicitCast<ValueType> explicit_cast;

    vigra_precondition(allLessEqual(1, shape),
        "vigra::detail::write_image_bands: invalid shape");

    encoder->setWidth(shape[1]);
    encoder->setHeight(shape[0]);
    encoder->setNumBands(shape[2] + 1);
    encoder->finalizeSettings();

    const unsigned offset(encoder->getOffset()); // correct offset only _after_ finalizeSettings()

                                                    // OPTIMIZATION: Specialization for the most common cases
                                                    // (1 or 3 channels)
    if (shape[2] == 1)
    {
        for (ArrayIndex y = 0; y != shape[0]; ++y, p.inc(0), alpha.inc(0))
        {
            auto scanline0 = static_cast<ValueType*>(encoder->currentScanlineOfBand(0));
            auto scanline1 = static_cast<ValueType*>(encoder->currentScanlineOfBand(1));

            for (ArrayIndex x = 0; x < shape[1]; ++x, p.inc(1), alpha.inc(1))
            {
                *scanline0 = explicit_cast::cast(image_scaler(*p));
                *scanline1 = explicit_cast::cast(alpha_scaler(*alpha));

                scanline0 += offset;
                scanline1 += offset;
            }

            p.move(1, -shape[1]);
            alpha.move(1, -shape[1]);

            encoder->nextScanline();
        }
    }
    else if (shape[2] == 3)
    {
        ValueType* scanline_0;
        ValueType* scanline_1;
        ValueType* scanline_2;
        ValueType* scanline_3;

        for (ArrayIndex y = 0; y != shape[0]; ++y, p.inc(0), alpha.inc(0))
        {
            scanline_0 = static_cast<ValueType*>(encoder->currentScanlineOfBand(0));
            scanline_1 = static_cast<ValueType*>(encoder->currentScanlineOfBand(1));
            scanline_2 = static_cast<ValueType*>(encoder->currentScanlineOfBand(2));
            scanline_3 = static_cast<ValueType*>(encoder->currentScanlineOfBand(3));

            for (ArrayIndex x = 0; x < shape[1]; ++x, p.inc(1))
            {
                *scanline_0 = explicit_cast::cast(image_scaler(*p));
                p.inc(2);
                *scanline_1 = explicit_cast::cast(image_scaler(*p));
                p.inc(2);
                *scanline_2 = explicit_cast::cast(image_scaler(*p));
                p.move(2, -2);

                *scanline_3 = explicit_cast::cast(alpha_scaler(*alpha));

                scanline_0 += offset;
                scanline_1 += offset;
                scanline_2 += offset;
                scanline_3 += offset;
            }

            p.move(1, -shape[1]);
            alpha.move(1, -shape[1]);

            encoder->nextScanline();
        }
    }
    else
    {
        std::vector<ValueType*> scanlines(shape[2] + 1);

        for (ArrayIndex y = 0; y != shape[0]; ++y, p.inc(0), alpha.inc(0))
        {
            for (ArrayIndex i = 0; i != shape[2]; ++i)
            {
                scanlines[i] = static_cast<ValueType*>(encoder->currentScanlineOfBand(i));
            }
            scanlines[shape[2]] = static_cast<ValueType*>(encoder->currentScanlineOfBand(shape[2]));

            for (ArrayIndex x = 0; x < shape[1]; ++x, p.inc(1), alpha.inc(1))
            {
                for (ArrayIndex i = 0; i != shape[2]; ++i, p.inc(2))
                {
                    *scanlines[i] = explicit_cast::cast(image_scaler(*p));
                    scanlines[i] += offset;
                }
                p.move(2, -shape[2]);

                *scanlines[shape[2]] = explicit_cast::cast(alpha_scaler(*alpha));
                scanlines[shape[2]] += offset;
            }

            p.move(1, -shape[1]);
            alpha.move(1, -shape[1]);

            encoder->nextScanline();
        }
    }
}

template <int N, class T1, class T2>
void
    exportImageAlpha(ArrayViewND<N, T1> view, ArrayViewND<N, T2> alpha,
        const ImageExportInfo& export_info)
{
    std::unique_ptr<Encoder> encoder(vigra::encoder(export_info));

    // FIXME: ensure C order.
    auto ca = view.ensureChannelAxis(2);
    vigra_precondition(ca.ndim() == 3,
        "exportImageAlpha(): wrong number of axes in source array.");
    auto p = ca.pointer_nd();
    ArrayIndex channels = ca.shape(2);
    vigra_precondition(isBandNumberSupported(encoder->getFileType(), channels + 1),
        "exportImageAlpha(): file format does not support requested number of bands (color channels).");

    auto aa = alpha.ensureChannelAxis(2);
    vigra_precondition(aa.ndim() == 3,
        "exportImageAlpha(): wrong number of axes in alpha arra.y");
    vigra_precondition(aa.shape(2) == 1,
        "exportImageAlpha(): alpha array must have exactly one channel.");
    auto a = aa.pointer_nd();

    using value_type = typename std::decay<typename decltype(ca)::value_type>::type;

    std::string pixel_type = export_info.getPixelType();
    const bool downcast = negotiatePixelType(encoder->getFileType(),
        TypeAsString<value_type>::result(), pixel_type);
    const pixel_t type = pixel_t_of_string(pixel_type);

    encoder->setPixelType(pixel_type);

    const range_t image_source_range = find_source_value_range(export_info, ca);
    const range_t alpha_source_range = find_source_value_range(export_info, aa);
    const range_t destination_range = find_destination_value_range(export_info, type);

    if ((downcast || export_info.hasForcedRangeMapping()) &&
        (image_source_range.first != destination_range.first ||
            image_source_range.second != destination_range.second ||
            alpha_source_range.first != destination_range.first ||
            alpha_source_range.second != destination_range.second))
    {
        const linear_transform image_rescaler(image_source_range, destination_range);
        const linear_transform alpha_rescaler(alpha_source_range, destination_range);

        switch (type)
        {
        case UNSIGNED_INT_8:
            write_image_bands_and_alpha<uint8_t>(encoder.get(), p, a, ca.shape(), image_rescaler, alpha_rescaler);
            break;
        case UNSIGNED_INT_16:
            write_image_bands_and_alpha<uint16_t>(encoder.get(), p, a, ca.shape(), image_rescaler, alpha_rescaler);
            break;
        case UNSIGNED_INT_32:
            write_image_bands_and_alpha<uint32_t>(encoder.get(), p, a, ca.shape(), image_rescaler, alpha_rescaler);
            break;
        case SIGNED_INT_16:
            write_image_bands_and_alpha<int16_t>(encoder.get(), p, a, ca.shape(), image_rescaler, alpha_rescaler);
            break;
        case SIGNED_INT_32:
            write_image_bands_and_alpha<int32_t>(encoder.get(), p, a, ca.shape(), image_rescaler, alpha_rescaler);
            break;
        case IEEE_FLOAT_32:
            write_image_bands_and_alpha<float>(encoder.get(), p, a, ca.shape(), image_rescaler, alpha_rescaler);
            break;
        case IEEE_FLOAT_64:
            write_image_bands_and_alpha<double>(encoder.get(), p, a, ca.shape(), image_rescaler, alpha_rescaler);
            break;
        default:
            vigra_fail("vigra::detail::exportImageAlpha(): should never be reached");
        }
    }
    else
    {
        switch (type)
        {
        case UNSIGNED_INT_8:
            write_image_bands_and_alpha<uint8_t>(encoder.get(), p, a, ca.shape(), identity(), identity());
            break;
        case UNSIGNED_INT_16:
            write_image_bands_and_alpha<uint16_t>(encoder.get(), p, a, ca.shape(), identity(), identity());
            break;
        case UNSIGNED_INT_32:
            write_image_bands_and_alpha<uint32_t>(encoder.get(), p, a, ca.shape(), identity(), identity());
            break;
        case SIGNED_INT_16:
            write_image_bands_and_alpha<int16_t>(encoder.get(), p, a, ca.shape(), identity(), identity());
            break;
        case SIGNED_INT_32:
            write_image_bands_and_alpha<int32_t>(encoder.get(), p, a, ca.shape(), identity(), identity());
            break;
        case IEEE_FLOAT_32:
            write_image_bands_and_alpha<float>(encoder.get(), p, a, ca.shape(), identity(), identity());
            break;
        case IEEE_FLOAT_64:
            write_image_bands_and_alpha<double>(encoder.get(), p, a, ca.shape(), identity(), identity());
            break;
        default:
            vigra_fail("vigra::detail::exportImageAlpha(): should never be reached");
        }
    }

    encoder->close();
}

} // end namespace detail


/** \brief Read the image specified by the given \ref
    vigra::ImageImportInfo object including its alpha channel.

    See \ref importImage() for more information.

    <B>Declarations</B>

    pass 2D array views:
    \code
    namespace vigra {
    template <class T1, class T2>
    void
    importImageAlpha(ImageImportInfo const & import_info,
    ArrayViewND<2,T1> image,
    ArrayViewND<2,T2> alpha);
    }
    \endcode

    \deprecatedAPI{importImageAlpha}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
    template <class ImageIterator, class ImageAccessor,
    class AlphaIterator, class AlphaAccessor>
    void
    importImageAlpha(const ImageImportInfo& importInfo,
    ImageIterator imageIterator, ImageAccessor imageAccessor,
    AlphaIterator alphaIterator, AlphaAccessor alphaAccessor)
    }
    \endcode
    Use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
    template <class ImageIterator, class ImageAccessor,
    class AlphaIterator, class AlphaAccessor>
    void
    importImageAlpha(const ImageImportInfo& importInfo,
    const std::pair<ImageIterator, ImageAccessor>& image,
    const std::pair<AlphaIterator, AlphaAccessor>& alpha)
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <B>\#include</B> \<vigra/impexalpha.hxx\><br/>
    Namespace: vigra

    \code
    typedef uint8_t value_t;
    ImageImportInfo info("zorro.tif");

    if (info.isGrayscale())
    {
    MultiArray<2, value_t> alpha(info.shape());
    MultiArray<2, value_t> image(info.shape());

    importImageAlpha(info, image, alpha);
    ...
    }
    else
    {
    MultiArray<2, value_t>            alpha(info.shape());
    MultiArray<2, RGBValue<value_t> > image(info.shape());

    importImageAlpha(info, image, alpha);
    ...
    }
    \endcode

    \deprecatedUsage{importImageAlpha}
    \code
    typedef uint8_t value_t;
    ImageImportInfo info("zorro.tif");

    if (info.isGrayscale())
    {
    BasicImage<value_t> alpha(info.size());
    BasicImage<value_t> image(info.size());

    importImageAlpha(info,
    image.upperLeft(), image.accessor(),
    alpha.upperLeft(), alpha.accessor());
    ...
    }
    else
    {
    BasicImage<value_t> alpha(info.size());
    BasicImage<vigra::RGBValue<value_t> > image(info.size());

    importImageAlpha(info,
    image.upperLeft(), image.accessor(),
    alpha.upperLeft(), alpha.accessor());
    ...
    }
    \endcode
    \deprecatedEnd

    <B>Preconditions</B>

    - The same preconditions hold as for importImage(), however the
    only image formats that support alpha channels are
    + TIFF and
    + PNG.
    In particular, JPEG does <B>not</B> support alpha channels.
    - The alpha channel always is scalar-valued, i.e. comprises a
    single band.
*/
    doxygen_overloaded_function(template <...> void importImageAlpha)

template <class T1, class T2>
inline void
importImageAlpha(ImageImportInfo const & import_info,
                 ArrayViewND<2, T1> image, ArrayViewND<2, T2> alpha)
{
    // TODO get the memory order from image also considering axistags,
    // if available.
    vigra_precondition(import_info.shape() == image.shape(),
        "importImageAlpha(): shape mismatch between input and output.");
    vigra_precondition(import_info.shape() == alpha.shape(),
        "importImageAlpha(): shape mismatch between input and alpha array.");
    detail::importImageAlpha(import_info, image, alpha);
}

template <class T1, class T2, class A1, class A2>
inline void
importImageAlpha(char const * name,
                 ArrayND<2, T1, A1> & image, ArrayND<2, T2, A2> & alpha)
{
    ImageImportInfo info(name);
    image.resize(info.shape());
    alpha.resize(info.shape());
    detail::importImageAlpha(info, image, alpha);
}

template <class T1, class T2, class A1, class A2>
inline void
importImageAlpha(std::string const & name,
                 ArrayND<2, T1, A1> & image, ArrayND<2, T2, A2> & alpha)
{
    importImageAlpha(name.c_str(), image, alpha);
}

// FIXME: implement importImageAlpha(ArrayViewND<3, T>, ...)

/** \brief Write the image and its alpha channel to a file.

    See \ref exportImage() for more information.

    <B>Declarations</B>

    pass 2D array views:
    \code
    namespace vigra {
    template <class T1, class T2>
    void
    exportImageAlpha(ArrayViewND<2,T1> const & image,
    ArrayViewND<2,T2> const & alpha,
    ImageExportInfo const & export_info);

    template <class T1, class T2>
    void
    exportImageAlpha(ArrayViewND<2,T1> const & image,
    ArrayViewND<2,T2> const & alpha,
    char const * filename)

    template <class T1, class T2>
    void
    exportImageAlpha(ArrayViewND<2,T1> const & image,
    ArrayViewND<2,T2> const & alpha,
    std::string const & filename)
    }
    \endcode

    \deprecatedAPI{exportImageAlpha}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
    template <class ImageIterator, class ImageAccessor,
    class AlphaIterator, class AlphaAccessor>
    void
    exportImageAlpha(ImageIterator imageUpperLeft, ImageIterator imageLowerRight, ImageAccessor imageAccessor,
    AlphaIterator alphaUpperLeft, AlphaAccessor alphaAccessor,
    const ImageExportInfo& exportInfo)
    }
    \endcode
    Use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
    template <class ImageIterator, class ImageAccessor,
    class AlphaIterator, class AlphaAccessor>
    void
    exportImageAlpha(const std::tuple<ImageIterator, ImageIterator, ImageAccessor>& image,
    const std::pair<AlphaIterator, AlphaAccessor>& alpha,
    const ImageExportInfo& exportInfo)
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <B>\#include</B> \<vigra/impexalpha.hxx\><br/>
    Namespace: vigra

    \code
    typedef uint8_t value_t;

    MultiArray<2, value_t>            alpha(width, height);
    MultiArray<2, RGBValue<value_t> > image(width, height);

    ... // do some image processing

    // specify the output filename
    exportImageAlpha(image, alpha, "zorro.tif");

    // use a ImageExportInfo if you need more control over the export
    exportImageAlpha(image, alpha, ImageExportInfo("zorro.tif").setPixelType("FLOAT"));
    \endcode

    \deprecatedUsage{exportImageAlpha}
    \code
    typedef uint8_t value_t;
    ImageExportInfo info("zorro.tif");

    if (info.isGrayscale())
    {
    BasicImage<value_t> alpha;
    BasicImage<value_t> image;

    ...

    exportImageAlpha(image.upperLeft(), image.lowerRight(), image.accessor(),
    alpha.upperLeft(), alpha.accessor(),
    info);
    }
    else
    {
    BasicImage<value_t> alpha;
    BasicImage<vigra::RGBValue<value_t> > image;

    ...

    exportImageAlpha(image.upperLeft(), image.lowerRight(), image.accessor(),
    alpha.upperLeft(), alpha.accessor(),
    info);
    }
    \endcode
    \deprecatedEnd

    <B>Preconditions</B>

    - The same preconditions hold as for exportImage(), however the
    only image formats that support alpha channels are
    + TIFF and
    + PNG.
    In particular, JPEG does <B>not</B> support alpha channels.
    - The alpha channel always is scalar-valued, i.e. comprises a
    single band.
*/
doxygen_overloaded_function(template <...> void exportImageAlpha)

template <class T1, class T2>
inline void
exportImageAlpha(ArrayViewND<2, T1> const & image,
                 ArrayViewND<2, T2> const & alpha,
                 ImageExportInfo const & export_info)
{
    // FIXME: ensure C-order, possibly using axistags
    try
    {
        detail::exportImageAlpha(image, alpha, export_info);
    }
    catch (Encoder::TIFFCompressionException&)
    {
        ImageExportInfo info(export_info);
        info.setCompression("");
        detail::exportImageAlpha(image, alpha, info);
    }
}

template <class T1, class T2>
inline void
exportImageAlpha(ArrayViewND<2, T1> const & image,
                 ArrayViewND<2, T2> const & alpha,
                 char const * name)
{
    ImageExportInfo export_info(name);
    exportImageAlpha(image, alpha, export_info);
}

template <class T1, class T2>
inline void
exportImageAlpha(ArrayViewND<2, T1> const & image,
                 ArrayViewND<2, T2> const & alpha,
                 std::string const & name)
{
    ImageExportInfo export_info(name.c_str());
    exportImageAlpha(image, alpha, export_info);
}

// FIXME: implement exportImageAlpha(ArrayViewND<3, T>, ...)

/** @} */

} // end namespace vigra

#endif // VIGRA2_IMAGEIO_IMPEXALPHA_HXX
