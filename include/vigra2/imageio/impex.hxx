/************************************************************************/
/*                                                                      */
/*               Copyright 2001-2002 by Gunnar Kedenburg                */
/*        Copyright 2012 Christoph Spiel and Ullrich Koethe             */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
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


/**
 * \file  impex.hxx
 * \brief image import and export functions
 *
 * This module implements functions importImage() and exportImage().
 * The matching implementation for any given datatype is selected by
 * template meta code.
 *
 */

#ifndef VIGRA2_IMAGEIO_IMPEX_HXX
#define VIGRA2_IMAGEIO_IMPEX_HXX

#include <vigra2/config.hxx>
#include <vigra2/shape.hxx>
#include <vigra2/sized_int.hxx>
#include <vigra2/array_nd.hxx>

#include <tuple>
#include <type_traits>
#include <utility>

#include "imageinfo.hxx"
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
        read_image_bands(Decoder* decoder, PointerND<N,T> p, ArrayIndex channels)
        {
            vigra_assert(p.ndim() == 3, "Number of dimensions must be 3.");
            const ArrayIndex width(decoder->getWidth());
            const ArrayIndex height(decoder->getHeight());
            const ArrayIndex bands(decoder->getNumBands());
            const ArrayIndex offset(decoder->getOffset());

            // OPTIMIZATION: Specialization for the most common case
            // of an RGB-image, i.e. 3 channels.
            if (channels == 1)
            {
                for (ArrayIndex y = 0; y != height; ++y, p.inc(0))
                {
                    decoder->nextScanline();

                    auto scanline = static_cast<const ValueType*>(decoder->currentScanlineOfBand(0));

                    for (ArrayIndex x = 0; x < width; ++x, p.inc(1), scanline += offset)
                    {
                        *p = *scanline;
                    }

                    p.move(1,-width);
                }
            }
            else if (channels == 3)
            {
                const ValueType* scanline_0;
                const ValueType* scanline_1;
                const ValueType* scanline_2;

                for (ArrayIndex y = 0; y != height; ++y, p.inc(0))
                {
                    decoder->nextScanline();

                    scanline_0 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(0));

                    if(bands == 1)
                    {
                        scanline_1 = scanline_0;
                        scanline_2 = scanline_0;
                    }
                    else
                    {
                        scanline_1 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(1));
                        scanline_2 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(2));
                    }

                    for (ArrayIndex x = 0; x < width; ++x, p.inc(1))
                    {
                        *p = *scanline_0;
                        p.inc(2);
                        *p = *scanline_1;
                        p.inc(2);
                        *p = *scanline_2;
                        p.move(2,-2);

                        scanline_0 += offset;
                        scanline_1 += offset;
                        scanline_2 += offset;
                    }

                    p.move(1,-width);
                }
            }
            else
            {
                std::vector<const ValueType*> scanlines(channels);

                for (ArrayIndex y = 0; y != height; ++y, p.inc(0))
                {
                    decoder->nextScanline();

                    scanlines[0] = static_cast<const ValueType*>(decoder->currentScanlineOfBand(0));

                    if(bands == 1)
                    {
                        for (ArrayIndex i = 1; i != channels; ++i)
                        {
                            scanlines[i] = scanlines[0];
                        }
                    }
                    else
                    {
                        for (ArrayIndex i = 1; i != channels; ++i)
                        {
                            scanlines[i] = static_cast<const ValueType*>(decoder->currentScanlineOfBand(i));
                        }
                    }

                    for (ArrayIndex x = 0; x < width; ++x, p.inc(1))
                    {
                        for (ArrayIndex i = 0; i != channels; ++i, p.inc(2))
                        {
                            *p = *scanlines[i];
                            scanlines[i] += offset;
                        }
                        p.move(2,-channels);
                    }

                    p.move(1,-width);
                }
            }
        }

        template <class ValueType, int N>
        void
        importImage(const ImageImportInfo& import_info, ArrayViewND<N,ValueType> view)
        {
            std::unique_ptr<Decoder> decoder(vigra::decoder(import_info));

            // FIXME: ensure C order.
            auto ca = view.ensureChannelAxis(2);
            vigra_precondition(ca.ndim() == 3, "importImage(): wrong number of axes in target array");
            vigra_precondition(ca.shape(0) == decoder->getHeight() && ca.shape(1) == decoder->getWidth(),
                "importImage(): shape mismatch between input and output.");
            vigra_precondition(ca.shape(2) == decoder->getNumBands() || decoder->getNumBands() == 1,
                "importImage(): Number of channels in input and destination image don't match.");
            auto p = ca.pointer_nd();
            ArrayIndex channels = ca.shape(2);

            switch (pixel_t_of_string(decoder->getPixelType()))
            {
            case UNSIGNED_INT_8:
                read_image_bands<uint8_t>(decoder.get(), p, channels);
                break;
            case UNSIGNED_INT_16:
                read_image_bands<uint16_t>(decoder.get(), p, channels);
                break;
            case UNSIGNED_INT_32:
                read_image_bands<uint32_t>(decoder.get(), p, channels);
                break;
            case SIGNED_INT_16:
                read_image_bands<int16_t>(decoder.get(), p, channels);
                break;
            case SIGNED_INT_32:
                read_image_bands<int32_t>(decoder.get(), p, channels);
                break;
            case IEEE_FLOAT_32:
                read_image_bands<float>(decoder.get(), p, channels);
                break;
            case IEEE_FLOAT_64:
                read_image_bands<double>(decoder.get(), p, channels);
                break;
            default:
                vigra_fail("detail::importImage<scalar>: not reached");
            }

            decoder->close();
        }

        template<class ValueType,
                 class ImageIterator, class ImageAccessor, class ImageScaler>
        void
        write_image_band(Encoder* encoder,
                         ImageIterator image_upper_left, ImageIterator image_lower_right, ImageAccessor image_accessor,
                         const ImageScaler& image_scaler)
        {
            typedef typename ImageIterator::row_iterator ImageRowIterator;

            typedef RequiresExplicitCast<ValueType> explicit_cast;

            vigra_precondition(image_lower_right.x >= image_upper_left.x,
                               "vigra::detail::write_image_band: negative width");
            vigra_precondition(image_lower_right.y >= image_upper_left.y,
                               "vigra::detail::write_image_band: negative height");

            const unsigned width(static_cast<unsigned>(image_lower_right.x - image_upper_left.x));
            const unsigned height(static_cast<unsigned>(image_lower_right.y - image_upper_left.y));

            encoder->setWidth(width);
            encoder->setHeight(height);
            encoder->setNumBands(1);
            encoder->finalizeSettings();

            const unsigned offset(encoder->getOffset()); // correct offset only _after_ finalizeSettings()

            // IMPLEMENTATION NOTE: We avoid calling the default
            // constructor to allow classes ImageIterator that do not
            // define one.
            ImageIterator image_iterator(image_upper_left);

            for (unsigned y = 0U; y != height; ++y)
            {
                ValueType* scanline = static_cast<ValueType*>(encoder->currentScanlineOfBand(0));

                ImageRowIterator is(image_iterator.rowIterator());
                const ImageRowIterator is_end(is + width);

                while (is != is_end)
                {
                    *scanline = explicit_cast::cast(image_scaler(image_accessor(is)));
                    scanline += offset;
                    ++is;
                }

                encoder->nextScanline();

                ++image_iterator.y;
            }
        }


        template<class ValueType,
                 class ImageIterator, class ImageAccessor, class ImageScaler>
        void
        write_image_bands(Encoder* encoder,
                          ImageIterator image_upper_left, ImageIterator image_lower_right, ImageAccessor image_accessor,
                          const ImageScaler& image_scaler)
        {
            typedef typename ImageIterator::row_iterator ImageRowIterator;
            typedef RequiresExplicitCast<ValueType> explicit_cast;

            vigra_precondition(image_lower_right.x >= image_upper_left.x,
                               "vigra::detail::write_image_bands: negative width");
            vigra_precondition(image_lower_right.y >= image_upper_left.y,
                               "vigra::detail::write_image_bands: negative height");

            const unsigned width(static_cast<unsigned>(image_lower_right.x - image_upper_left.x));
            const unsigned height(static_cast<unsigned>(image_lower_right.y - image_upper_left.y));
            const unsigned accessor_size(image_accessor.size(image_upper_left));

            encoder->setWidth(width);
            encoder->setHeight(height);
            encoder->setNumBands(accessor_size);
            encoder->finalizeSettings();

            const unsigned offset(encoder->getOffset()); // correct offset only _after_ finalizeSettings()

            // IMPLEMENTATION NOTE: We avoid calling the default
            // constructor to allow classes ImageIterator that do not
            // define one.
            ImageIterator image_iterator(image_upper_left);

            // OPTIMIZATION: Specialization for the most common case
            // of an RGB-image, i.e. 3 channels.
            if (accessor_size == 3U)
            {
                ValueType* scanline_0;
                ValueType* scanline_1;
                ValueType* scanline_2;

                for (unsigned y = 0U; y != height; ++y)
                {
                    scanline_0 = static_cast<ValueType*>(encoder->currentScanlineOfBand(0));
                    scanline_1 = static_cast<ValueType*>(encoder->currentScanlineOfBand(1));
                    scanline_2 = static_cast<ValueType*>(encoder->currentScanlineOfBand(2));

                    ImageRowIterator is(image_iterator.rowIterator());
                    const ImageRowIterator is_end(is + width);

                    while (is != is_end)
                    {
                        *scanline_0 = explicit_cast::cast(image_scaler(image_accessor.getComponent(is, 0)));
                        *scanline_1 = explicit_cast::cast(image_scaler(image_accessor.getComponent(is, 1)));
                        *scanline_2 = explicit_cast::cast(image_scaler(image_accessor.getComponent(is, 2)));

                        scanline_0 += offset;
                        scanline_1 += offset;
                        scanline_2 += offset;

                        ++is;
                    }

                    encoder->nextScanline();

                    ++image_iterator.y;
                }
            }
            else
            {
                std::vector<ValueType*> scanlines(accessor_size);

                for (unsigned y = 0U; y != height; ++y)
                {
                    for (unsigned i = 0U; i != accessor_size; ++i)
                    {
                        scanlines[i] = static_cast<ValueType*>(encoder->currentScanlineOfBand(i));
                    }

                    ImageRowIterator is(image_iterator.rowIterator());
                    const ImageRowIterator is_end(is + width);

                    while (is != is_end)
                    {
                        for (unsigned i = 0U; i != accessor_size; ++i)
                        {
                            *scanlines[i] = explicit_cast::cast(image_scaler(image_accessor.getComponent(is, static_cast<int>(i))));
                            scanlines[i] += offset;
                        }
                        ++is;
                    }

                    encoder->nextScanline();

                    ++image_iterator.y;
                }
            }
        }

        template<class ValueType, class T, int N, class ImageScaler>
        void
        write_image_bands(Encoder* encoder, PointerND<N,T> p, const Shape<N> &shape, const ImageScaler& image_scaler)
        {
            typedef RequiresExplicitCast<ValueType> explicit_cast;

            vigra_precondition(allLessEqual(1,shape),
                               "vigra::detail::write_image_bands: invalid shape");

            encoder->setWidth(shape[1]);
            encoder->setHeight(shape[0]);
            encoder->setNumBands(shape[2]);
            encoder->finalizeSettings();

            const unsigned offset(encoder->getOffset()); // correct offset only _after_ finalizeSettings()

            // OPTIMIZATION: Specialization for the most common case
            // of an RGB-image, i.e. 3 channels.
            if (shape[2] == 1)
            {
                for (ArrayIndex y = 0; y != shape[0]; ++y, p.inc(0))
                {
                    auto scanline = static_cast<ValueType*>(encoder->currentScanlineOfBand(0));

                    for (ArrayIndex x = 0; x < shape[1]; ++x, p.inc(1), scanline += offset)
                    {
                        *scanline = explicit_cast::cast(image_scaler(*p));
                    }

                    p.move(1,-shape[1]);

                    encoder->nextScanline();
                }
            }
            else if (shape[2] == 3)
            {
                ValueType* scanline_0;
                ValueType* scanline_1;
                ValueType* scanline_2;

                for (ArrayIndex y = 0; y != shape[0]; ++y, p.inc(0))
                {
                    scanline_0 = static_cast<ValueType*>(encoder->currentScanlineOfBand(0));
                    scanline_1 = static_cast<ValueType*>(encoder->currentScanlineOfBand(1));
                    scanline_2 = static_cast<ValueType*>(encoder->currentScanlineOfBand(2));

                    for (ArrayIndex x = 0; x < shape[1]; ++x, p.inc(1))
                    {
                        *scanline_0 = explicit_cast::cast(image_scaler(*p));
                        p.inc(2);
                        *scanline_1 = explicit_cast::cast(image_scaler(*p));
                        p.inc(2);
                        *scanline_2 = explicit_cast::cast(image_scaler(*p));
                        p.move(2,-2);

                        scanline_0 += offset;
                        scanline_1 += offset;
                        scanline_2 += offset;
                    }

                    p.move(1,-shape[1]);

                    encoder->nextScanline();
                }
            }
            else
            {
                std::vector<ValueType*> scanlines(shape[2]);

                for (ArrayIndex y = 0; y != shape[0]; ++y, p.inc(0))
                {
                    for (ArrayIndex i = 0; i != shape[2]; ++i)
                    {
                        scanlines[i] = static_cast<ValueType*>(encoder->currentScanlineOfBand(i));
                    }

                    for (ArrayIndex x = 0; x < shape[1]; ++x, p.inc(1))
                    {
                        for (ArrayIndex i = 0; i != shape[2]; ++i, p.inc(2))
                        {
                            *scanlines[i] = explicit_cast::cast(image_scaler(*p));
                            scanlines[i] += offset;
                        }

                        p.move(2,-shape[2]);
                    }

                    p.move(1,-shape[1]);

                    encoder->nextScanline();
                }
            }
        }

        template <int N, class T>
        void
        exportImage(ArrayViewND<N,T> view, const ImageExportInfo& export_info)
        {
            std::unique_ptr<Encoder> encoder(vigra::encoder(export_info));

            // FIXME: ensure C order.
            auto ca = view.ensureChannelAxis(2);
            vigra_precondition(ca.ndim() == 3, "exportImage(): wrong number of axes in source array");
            auto p = ca.pointer_nd();
            ArrayIndex channels = ca.shape(2);

            using value_type = typename std::decay<typename decltype(ca)::value_type>::type;

            std::string pixel_type = export_info.getPixelType();
            const bool downcast    = negotiatePixelType(encoder->getFileType(),
                                                        TypeAsString<value_type>::result(), pixel_type);
            const pixel_t type     = pixel_t_of_string(pixel_type);

            encoder->setPixelType(pixel_type);

            vigra_precondition(isBandNumberSupported(encoder->getFileType(), channels),
                               "exportImage(): file format does not support requested number of bands (color channels)");

            const range_t image_source_range = find_source_value_range(export_info, ca);
            const range_t destination_range  = find_destination_value_range(export_info, type);

            if ((downcast || export_info.hasForcedRangeMapping()) &&
                (image_source_range.first != destination_range.first || image_source_range.second != destination_range.second))
            {
                const linear_transform image_rescaler(image_source_range, destination_range);

                switch (type)
                {
                case UNSIGNED_INT_8:
                    write_image_bands<uint8_t>(encoder.get(), p, ca.shape(), image_rescaler);
                    break;
                case UNSIGNED_INT_16:
                    write_image_bands<uint16_t>(encoder.get(), p, ca.shape(), image_rescaler);
                    break;
                case UNSIGNED_INT_32:
                    write_image_bands<uint32_t>(encoder.get(), p, ca.shape(), image_rescaler);
                    break;
                case SIGNED_INT_16:
                    write_image_bands<int16_t>(encoder.get(), p, ca.shape(), image_rescaler);
                    break;
                case SIGNED_INT_32:
                    write_image_bands<int32_t>(encoder.get(), p, ca.shape(), image_rescaler);
                    break;
                case IEEE_FLOAT_32:
                    write_image_bands<float>(encoder.get(), p, ca.shape(), image_rescaler);
                    break;
                case IEEE_FLOAT_64:
                    write_image_bands<double>(encoder.get(), p, ca.shape(), image_rescaler);
                    break;
                default:
                    vigra_fail("vigra::detail::exportImage<scalar>: not reached");
                }
            }
            else
            {
                switch (type)
                {
                case UNSIGNED_INT_8:
                    write_image_bands<uint8_t>(encoder.get(), p, ca.shape(), identity());
                    break;
                case UNSIGNED_INT_16:
                    write_image_bands<uint16_t>(encoder.get(), p, ca.shape(), identity());
                    break;
                case UNSIGNED_INT_32:
                    write_image_bands<uint32_t>(encoder.get(), p, ca.shape(), identity());
                    break;
                case SIGNED_INT_16:
                    write_image_bands<int16_t>(encoder.get(), p, ca.shape(), identity());
                    break;
                case SIGNED_INT_32:
                    write_image_bands<int32_t>(encoder.get(), p, ca.shape(), identity());
                    break;
                case IEEE_FLOAT_32:
                    write_image_bands<float>(encoder.get(), p, ca.shape(), identity());
                    break;
                case IEEE_FLOAT_64:
                    write_image_bands<double>(encoder.get(), p, ca.shape(), identity());
                    break;
                default:
                    vigra_fail("vigra::detail::exportImage<scalar>: not reached");
                }
            }

            encoder->close();
        }


        template <class ImageIterator, class ImageAccessor>
        void
        exportImage(ImageIterator image_upper_left, ImageIterator image_lower_right, ImageAccessor image_accessor,
                    const ImageExportInfo& export_info,
                    /* isScalar? */ std::false_type)
        {
            typedef typename ImageAccessor::value_type ImageBaseType;
            typedef typename ImageBaseType::value_type ImageValueType;

            VIGRA_UNIQUE_PTR<Encoder> encoder(vigra::encoder(export_info));

            std::string pixel_type(export_info.getPixelType());
            const bool downcast(negotiatePixelType(encoder->getFileType(), TypeAsString<ImageValueType>::result(), pixel_type));
            const pixel_t type(pixel_t_of_string(pixel_type));

            encoder->setPixelType(pixel_type);

            vigra_precondition(isBandNumberSupported(encoder->getFileType(), image_accessor.size(image_upper_left)),
                               "exportImage(): file format does not support requested number of bands (color channels)");

            const range_t image_source_range(find_source_value_range(export_info,
                                                                     image_upper_left, image_lower_right, image_accessor));
            const range_t destination_range(find_destination_value_range(export_info, pixel_t_of_string(pixel_type)));

            if ((downcast || export_info.hasForcedRangeMapping()) &&
                (image_source_range.first != destination_range.first || image_source_range.second != destination_range.second))
            {
                const linear_transform image_rescaler(image_source_range, destination_range);

                switch (type)
                {
                case UNSIGNED_INT_8:
                    write_image_bands<uint8_t>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case UNSIGNED_INT_16:
                    write_image_bands<uint16_t>(encoder.get(),
                                              image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case UNSIGNED_INT_32:
                    write_image_bands<uint32_t>(encoder.get(),
                                              image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case SIGNED_INT_16:
                    write_image_bands<int16_t>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case SIGNED_INT_32:
                    write_image_bands<int32_t>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case IEEE_FLOAT_32:
                    write_image_bands<float>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case IEEE_FLOAT_64:
                    write_image_bands<double>(encoder.get(),
                                              image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                default:
                    vigra_fail("vigra::detail::exportImage<non-scalar>: not reached");
                }
            }
            else
            {
                switch (type)
                {
                case UNSIGNED_INT_8:
                    write_image_bands<uint8_t>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case UNSIGNED_INT_16:
                    write_image_bands<uint16_t>(encoder.get(),
                                              image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case UNSIGNED_INT_32:
                    write_image_bands<uint32_t>(encoder.get(),
                                              image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case SIGNED_INT_16:
                    write_image_bands<int16_t>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case SIGNED_INT_32:
                    write_image_bands<int32_t>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case IEEE_FLOAT_32:
                    write_image_bands<float>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case IEEE_FLOAT_64:
                    write_image_bands<double>(encoder.get(),
                                              image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                default:
                    vigra_fail("vigra::detail::exportImage<non-scalar>: not reached");
                }
            }

            encoder->close();
        }
    }  // end namespace detail

    /**
    \brief Read an image from a file.

    If the first parameter is \ref vigra::ImageImportInfo, this function assumes that the destination
    image has already the appropriate shape. If the first parameter is a string, the destination
    must be a \ref vigra::ArrayND reference, which will be reshaped automatically.

    If the input image has only a single band, but the destination has multiple bands (e.g. is an RGB
    image), all bands will receive the same data. When a multi-band file is read into a single-band
    destination array, only the first band is read. Any other mismatch between the number of bands in
    input and output is an error and will throw a precondition exception.

    <B>Declarations</B>

    pass 2D array views:
    \code
    namespace vigra {
        // read the data into an array view of appropriate size
        template <class T, class S>
        void
        importImage(ImageImportInfo const & import_info,
                    ArrayViewND<2, T, S> image);

        // resize the given array and then read the data
        template <class T, class A>
        void
        importImage(char const * filename,
                    ArrayND<2, T, A> & image);

        template <class T, class A>
        void
        importImage(std::string const & filename,
                    ArrayND<2, T, A> & image);
    }
    \endcode

    \deprecatedAPI{importImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        importImage(const ImageImportInfo& importInfo,
                    ImageIterator imageIterator, Accessor imageAccessor)
    }
    \endcode
    Use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class ImageIterator, class Accessor>
        void
        importImage(const ImageImportInfo& importInfo,
                    const pair<ImageIterator, Accessor>& image)
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <B>\#include</B> \<vigra/impex.hxx\><br/>
    Namespace: vigra

    \code
    ImageImportInfo info("myimage.gif");

    if (info.isGrayscale())
    {
        // create byte image of appropriate size
        ArrayND<2, unsigned char> image(info.width(), info.height());

        importImage(info, image);
        ...
    }
    else
    {
        // create byte RGB image of appropriate size
        ArrayND<2, RGBValue<unsigned char> > image(info.width(), info.height());

        importImage(info, image);
        ...
    }
    \endcode
    When the type of input image is already known, this can be shortened:
    \code
    // create empty float image
    ArrayND<2, float> image;

    // resize image and read the data
    importImage("myimage.png", image);
    \endcode

    \deprecatedUsage{importImage}
    \code
        ImageImportInfo info("myimage.gif");

        if (info.isGrayscale())
        {
            // create byte image of appropriate size
            BImage image(info.width(), info.height());

            importImage(info, destImage(image));
            ...
        }
        else
        {
            // create byte RGB image of appropriate size
            BRGBImage image(info.width(), info.height());

            importImage(info, destImage(image));
            ...
        }
    \endcode
    \deprecatedEnd

    <B>Preconditions</B>

    - The image file must be readable.
    - The required support library must be installed (if the table doesn't specify an external library,
      VIGRA supports the format natively).
    - The file type must be one of the following:

    <table cellspacing="10">
    <tr align="left">
    <th>Type</th><th> Extension </th><th> Name                           </th><th> Support Library </th>
    </tr><tr>
       <td> BMP  </td><td> bmp       </td><td> Microsoft Windows bitmap image file                        </td><td> </td>
       </tr><tr>
       <td> EXR  </td><td> exr       </td><td> OpenEXR high dynamic range image format                    </td><td> libopenexr </td>
       </tr><tr>
       <td> GIF  </td><td> gif       </td><td> CompuServe graphics interchange format, 8-bit color        </td><td> </td>
       </tr><tr>
       <td> HDR  </td><td> hdr       </td><td> Radiance RGBE high dynamic range image format              </td><td> </td>
       </tr><tr>
       <td> JPEG </td><td> jpg, jpeg </td><td> Joint Photographic Experts Group JFIF format, 24-bit color </td><td> libjpeg </td>
       </tr><tr>
       <td> PBM  </td><td> pbm       </td><td> Portable bitmap format (black and white)                   </td><td> </td>
       </tr><tr>
       <td> PGM  </td><td> pgm       </td><td> Portable graymap format (gray scale)                       </td><td> </td>
       </tr><tr>
       <td> PNG  </td><td> png       </td><td> Portable Network Graphic                                   </td><td> libpng </td>
       </tr><tr>
       <td> PNM  </td><td> pnm       </td><td> Portable anymap                                            </td><td> </td>
       </tr><tr>
       <td> PPM  </td><td> ppm       </td><td> Portable pixmap format (color)                             </td><td> </td>
       </tr><tr>
       <td> SUN  </td><td> ras       </td><td> SUN Rasterfile                                             </td><td> </td>
       </tr><tr>
       <td> TIFF </td><td> tif, tiff </td><td> Tagged Image File Format                                   </td><td> libtiff </td>
       </tr><tr>
       <td> VIFF </td><td> xv        </td><td> Khoros Visualization image file                            </td><td> </td>
       </table>
*/
    doxygen_overloaded_function(template <...> void importImage)


    template <class ImageIterator, class ImageAccessor>
    inline void
    importImage(const ImageImportInfo& import_info,
                ImageIterator image_iterator, ImageAccessor image_accessor)
    {
        typedef typename ImageAccessor::value_type ImageValueType;
        typedef typename NumericTraits<ImageValueType>::isScalar is_scalar;

        detail::importImage(import_info,
                    image_iterator, image_accessor,
                    is_scalar());
    }


    template <class ImageIterator, class ImageAccessor>
    inline void
    importImage(ImageImportInfo const & import_info,
                std::pair<ImageIterator, ImageAccessor> image)
    {
        importImage(import_info,
                    image.first, image.second);
    }

    template <class T>
    inline void
    importImage(ImageImportInfo const & import_info,
                ArrayViewND<2, T> image)
    {
        bool nontrivial_axistags = nontrivialAxisTags(image.axistags());
        vigra_precondition(nontrivial_axistags || image.axistags() == tags::axis_unknown,"importImage(): inconsistent axistags.");
        const auto p = nontrivial_axistags
            ? detail::permutationToOrder(image.axistags(),C_ORDER)
            : detail::permutationToOrder(image.strides(),C_ORDER);
        auto view = image.transpose(p);
        vigra_precondition(import_info.shape() == view.shape(),
            "importImage(): shape mismatch between input and output.");
        detail::importImage(import_info, view);
    }

    template <class T>
    inline void
    importImage(ImageImportInfo const & import_info,
                ArrayViewND<3, T> image)
    {
        vigra_precondition(image.hasChannelAxis(), "importImage(): channel axis must be marked.");
        bool nontrivial_axistags = nontrivialAxisTags(image.axistags());
        vigra_precondition(nontrivial_axistags || channelOnlyAxisTags(image.axistags()),"importImage(): inconsistent axistags.");
        const auto p = nontrivial_axistags
            ? detail::permutationToOrder(image.axistags(),C_ORDER)
            : detail::permutationToOrder(image.strides(),C_ORDER);
        auto view = image.transpose(p);
        vigra_precondition(import_info.shape() == view.shape().erase(view.channelAxis()),
            "importImage(): shape mismatch between input and output.");
        detail::importImage(import_info, view);
    }

    template <class T, class A>
    inline void
    importImage(ImageImportInfo const & info,
                ArrayND<2, T, A> & image,
                MemoryOrder order = C_ORDER)
    {
        image.resize(info.shape(order),defaultAxistags(2, false, order));
        importImage(info, image.view());
    }

    template <class T, class A>
    inline void
    importImage(std::string const & name,
                ArrayND<2, T, A> & image,
                MemoryOrder order = C_ORDER)
    {
        importImage(ImageImportInfo(name.c_str()), image, order);
    }

    template <class T, class A>
    inline void
    importImage(ImageImportInfo const & info,
                ArrayND<3, T, A> & image,
                MemoryOrder order = C_ORDER)
    {
        image.resize(info.shape(order).insert(order==C_ORDER ? 2 : 0, info.numBands()),defaultAxistags(3, true, order));
        importImage(info, image.view());
    }

    template <class T, class A>
    inline void
    importImage(std::string const & name,
                ArrayND<3, T, A> & image,
                MemoryOrder order = C_ORDER)
    {
        importImage(ImageImportInfo(name.c_str()), image, order);
    }

    /** \brief Write an image to a file.

    The file can be specified either by a file name or by a \ref vigra::ImageExportInfo object.
    In the latter case, you have much more control about how the file is written. By default,
    the file format to be created is guessed from the filename extension. This can be
    overridden by an explicit file type in the ImageExportInfo object. If the file format
    supports compression (e.g. JPEG and TIFF), default compression parameters are used
    which can be overridden by the ImageExportInfo object.

    If the file format to be created supports the pixel type of the source image, this
    pixel type will be kept in the file (e.g. <tt>float</tt> can be stored by TIFF without
    conversion) unless the ImageExportInfo object
    explicitly requests a different storage type. If the array's pixel type is not supported by
    the file format, the pixel values are transformed to the range 0..255 and
    converted to <tt>unsigned char</tt>, unless another mapping is explicitly requested by
    the ImageExportInfo object.

    Currently, the following file formats are supported.  The pixel types given in brackets
    are those that are written without conversion:
        - BMP: Microsoft Windows bitmap image file (pixel types: UINT8 as gray and RGB);
        - GIF: CompuServe graphics interchange format, 8-bit color (pixel types: UINT8 as gray and RGB);
        - JPEG: Joint Photographic Experts Group JFIF format, compressed 24-bit color
                (pixel types: UINT8 as gray and RGB), only available if libjpeg is installed;
        - PNG: Portable Network Graphic (pixel types: UINT8 and UINT16 with up to 4 channels),
                only available if libpng is installed;
        - PBM: Portable bitmap format (black and white);
        - PGM: Portable graymap format (pixel types: UINT8, INT16, INT32 as gray scale);
        - PNM: Portable anymap (pixel types: UINT8, INT16, INT32 as gray and RGB);
        - PPM: Portable pixmap format (pixel types: UINT8, INT16, INT32 as RGB);
        - SUN: SUN Rasterfile (pixel types: UINT8 as gray and RGB);
        - TIFF: Tagged Image File Format
              (pixel types: UINT8, INT16, INT32, FLOAT, DOUBLE with up to 4 channels),
              only available if libtiff is installed;
        - VIFF: Khoros Visualization image file
              (pixel types: UINT8, INT16, INT32, FLOAT, DOUBLE with arbitrary many channels);

    <B>Declarations</B>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T>
        void
        exportImage(ArrayViewND<2, T> const & image,
                    ImageExportInfo const & export_info);

        template <class T>
        void
        exportImage(ArrayViewND<2, T> const & image,
                    char const * filename);

        template <class T>
        void
        exportImage(ArrayViewND<2, T> const & image,
                    std::string const & filename);
    }
    \endcode

    \deprecatedAPI{exportImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class ImageIterator, class ImageAccessor>
        void
        exportImage(ImageIterator imageUpperLeft, ImageIterator imageLowerRight, ImageAccessor imageAccessor,
                    const ImageExportInfo& exportInfo)
    }
    \endcode
    Use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
        namespace vigra {
            template <class ImageIterator, class ImageAccessor>
            void exportImage(ImageIterator imageUpperLeft, ImageIterator imageLowerRight, ImageAccessor imageAccessor,
                             const ImageExportInfo& exportInfo)
        }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <B>\#include</B> \<vigra/impex.hxx\><br/>
    Namespace: vigra

    \code
    ArrayND<2, RGBValue<unsigned char> > image(width, height);
    ...

    // write as JPEG image, using compression quality 80
    exportImage(image,
                ImageExportInfo("my-image.jpg").setCompression("80"));

    // Force it to a particular pixel type.  The pixel type must be supported by the
    // desired image file format, otherwise an \ref vigra::PreconditionViolation
    // exception will be thrown.
    exportImage(image,
                ImageExportInfo("my-INT16-image.tif").setPixelType("INT16"));
    \endcode

    \deprecatedUsage{exportImage}
    \code
    BRGBImage image(width, height);
    ...

    // write as JPEG image, using compression quality 80
    exportImage(srcImageRange(image),
                ImageExportInfo("my-image.jpg").setCompression("80"));

    // Force it to a particular pixel type.  The pixel type must be supported by the
    // desired image file format, otherwise an \ref vigra::PreconditionViolation
    // exception will be thrown.
    exportImage(srcImageRange(image),
                ImageExportInfo("my-INT16-image.tif").setPixelType("INT16"));
    \endcode
    \deprecatedEnd

    <B>Preconditions</B>

    - The image file must be writable and
    - the file type must be one of the supported file types.
*/
    doxygen_overloaded_function(template <...> void exportImage)


    template <class ImageIterator, class ImageAccessor>
    inline void
    exportImage(ImageIterator image_upper_left, ImageIterator image_lower_right, ImageAccessor image_accessor,
                const ImageExportInfo& export_info)
    {
        typedef typename ImageAccessor::value_type ImageValueType;
        typedef typename NumericTraits<ImageValueType>::isScalar is_scalar;

        try
        {
            detail::exportImage(image_upper_left, image_lower_right, image_accessor,
                        export_info,
                        is_scalar());
        }
        catch (Encoder::TIFFCompressionException&)
        {
            ImageExportInfo info(export_info);

            info.setCompression("");
            detail::exportImage(image_upper_left, image_lower_right, image_accessor,
                                   info,
                                   is_scalar());
        }
    }

    template <class ImageIterator, class ImageAccessor>
    inline void
    exportImage(std::tuple<ImageIterator, ImageIterator, ImageAccessor> image,
                ImageExportInfo const & export_info)
    {
        exportImage(image.first, image.second, image.third,
                    export_info);
    }

    template <class T>
    inline void
    exportImage(ArrayViewND<2, T> const & image,
                ImageExportInfo const & export_info)
    {
        detail::exportImage(image, export_info);
    }

    template <class T>
    inline void
    exportImage(ArrayViewND<2, T> const & image,
                char const * name)
    {
        ImageExportInfo export_info(name);
        detail::exportImage(image, export_info);
    }

    template <class T>
    inline void
    exportImage(ArrayViewND<2, T> const & image,
                std::string const & name)
    {
        ImageExportInfo export_info(name.c_str());
        detail::exportImage(image, export_info);
    }

/** @} */

} // end namespace vigra

#endif // VIGRA_IMPEX_HXX
