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
/* Modifications by Pablo d'Angelo
 * updated to vigra 1.4 by Douglas Wilkins
 * as of 18 Febuary 2006:
 *  - Added support for obtaining extra bands beyond RGB.
 *  - Added support for a position field that indicates the start of this
 *    image relative to some global origin.
 *  - Added support for x and y resolution fields.
 * Modifications by Andrew Mihal, as of 16 October 2004
 *  - Added include for vigra/diff2d.hxx
 *  - Added support for ICC profiles
 */

#ifndef VIGRA2_IMAGEIO_IMPEX_TIFF_HXX
#define VIGRA2_IMAGEIO_IMPEX_TIFF_HXX

#include <vector>
#include <vigra2/shape.hxx>
#include <vigra2/imageio/codec.hxx>

namespace vigra {

    struct TIFFCodecFactory : public CodecFactory
    {
        CodecDesc getCodecDesc() const;
        VIGRA_UNIQUE_PTR<Decoder> getDecoder() const;
        VIGRA_UNIQUE_PTR<Encoder> getEncoder() const;
    };

    class TIFFDecoderImpl;
    class TIFFEncoderImpl;

    class TIFFDecoder : public Decoder
    {
        TIFFDecoderImpl * pimpl;

    public:

        TIFFDecoder() : pimpl(0) {}

        ~TIFFDecoder();

        std::string getFileType() const;
        unsigned int getWidth() const;
        unsigned int getHeight() const;
        unsigned int getNumBands() const;

        unsigned int getNumExtraBands() const;

        unsigned int getNumImages() const;
        void setImageIndex(unsigned int);
        unsigned int getImageIndex() const;

        Shape<2> getPosition() const;
        Shape<2> getCanvasSize() const;
        float getXResolution() const;
        float getYResolution() const;

        const void * currentScanlineOfBand( unsigned int ) const;
        void nextScanline();

        std::string getPixelType() const;
        unsigned int getOffset() const;

        void init( const std::string &, unsigned int );
        void init( const std::string & fileName)
        {
            init(fileName, 0);
        }

        void close();
        void abort();
    };

    class TIFFEncoder : public Encoder
    {
        TIFFEncoderImpl * pimpl;

    public:

        TIFFEncoder() : pimpl(0) {}

        ~TIFFEncoder();

        std::string getFileType() const;
        void setWidth( unsigned int );
        void setHeight( unsigned int );
        void setNumBands( unsigned int );

        void setCompressionType( const std::string &, int = -1 );
        void setPixelType( const std::string & );

        void setPosition( const Shape<2> & pos );
        void setCanvasSize( const Shape<2> & pos );
        void setXResolution( float xres );
        void setYResolution( float yres );

        unsigned int getOffset() const;

        void finalizeSettings();

        void * currentScanlineOfBand( unsigned int );
        void nextScanline();

        void setICCProfile(const ICCProfile & data);

        void init( const std::string &, const std::string & );
        void init( const std::string & fileName)
        {
            init(fileName, "w");
        }

        void close();
        void abort();
    };
}

#endif // VIGRA_IMPEX_TIFF_HXX
