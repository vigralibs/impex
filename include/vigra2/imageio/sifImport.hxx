/************************************************************************/
/*                                                                      */
/*       Copyright 2010 by Joachim Schleicher and Ullrich Koethe        */
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


/*
 *  Opens an Andor .sif file as MultiImageView.
 *  The width, height and number of images are extracted
 *  from the ASCII encoded variable length header.
 *
 *  Based on the Java-Code from
 *  http://rsb.info.nih.gov/ij/plugins/open-sif.html
 *  written by
 *  L. Stirling Churchman (stirling at stanford.edu)
 *  Philippe Carl (pcarl at uni-muenster.de)
 *  Yoshiyuki Arai (arai at phys1.med.osaka-u.ac.jp)
 *
 *  Currently tested SIF versions: 4.16.12005.0
 *                                 4.16.30001.0
 *                                 4. 6.    3.0
*/

#ifndef VIGRA_SIFIMPORT_HXX
#define VIGRA_SIFIMPORT_HXX

#include <fstream>
#include <cstring>
#include <cstddef>
#include <vector>

#include <vigra2/array_nd.hxx>

namespace vigra {


 /** \addtogroup VigraSIFImport Import of Images from Andor Cameras

    Read an Andor SIF file into a ArrayViewND.
*/
//@{

/********************************************************/
/*                                                      */
/*                   SIFImportInfo                      */
/*                                                      */
/********************************************************/
/** \brief Extracts image properties from an Andor SIF file header.

See \ref readSIF() for a usage example. This object must be
used to read the image header of an Andor SIF file
and enquire its properties.

<b>\#include</b> \<vigra/sifImport.hxx\><br>
Namespace: vigra
*/
class SIFImportInfo
{
    public:
        /** Construct SIFImportInfo object.

            The header of the Andor SIF file \a filename is accessed to
            read the image properties.

            \code
            SIFImportInfo info(filename);
            \endcode
         */
        VIGRA_EXPORT SIFImportInfo(const char* filename);

        /** Get the width in pixels.
         */
        VIGRA_EXPORT int width() const;

        /** Get the height in pixels.
         */
        VIGRA_EXPORT int height() const;

        /** Get the stacksize, that is the number of
            images contained in the dataset.
         */
        VIGRA_EXPORT int stacksize() const;

        /** Get the number of dimensions of the dataset represented by this info object.
         */
        VIGRA_EXPORT ArrayIndex numDimensions() const;

        /** Get the shape of the dataset represented by this info object.
         */
        VIGRA_EXPORT std::vector<size_t> const & shape() const;

        /** Get the shape (length) of the dataset along dimension \a dim.
         */
        VIGRA_EXPORT ArrayIndex shapeOfDimension(const int dim) const;

        /** Get the offset to the beginning of the actual data.
            Everything before this point belongs to the
            variable length header.
         */
        VIGRA_EXPORT std::ptrdiff_t getOffset() const;

        /** Get the filename of this SIF object.
         */
        VIGRA_EXPORT const char * getFileName() const;

        /** Output all information such as shutter, Temperature etc.
           as human readable output.

        <b> Usage:</b>

        <b>\#include</b> \<vigra/sifImport.hxx\><br>
        Namespace: vigra

        \code
        SIFImportInfo info(filename);
        std::cout << info << std::endl; // print infos to the console

        \endcode
         */
        VIGRA_EXPORT friend std::ostream& operator<<(std::ostream& os, const SIFImportInfo& info);

    private:
        const char* m_filename;
        std::vector<size_t> m_dims;
        std::ptrdiff_t m_offset;
        int mod;
        int left, right, bottom, top;
        int xbin, ybin, xres, yres;
        int headerlen;
        double readout;
        double temperature1, temperature2;
        long long d;
        std::string cycleTime, temperature, exposureTime, EMGain,
        verticalShiftSpeed, version, model, originalFilename, preAmpGain;
        size_t filesize;

};




    /** \brief Read the image data specified by the given \ref vigra::SIFImportInfo object
                and write them into the given 'array'.

    The array must have the correct number of dimensions and shape for the dataset
    represented by 'info'.

    <b> Declaration:</b>

    \code
    namespace vigra {
        void
        readSIF(const SIFImportInfo &info, ArrayViewND<3, float> array);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/sifImport.hxx\><br>
    Namespace: vigra

    \code
    SIFImportInfo info(filename);

    // create a 3D array of appropriate size
    typedef MultiArray<3, float>::difference_type Shape;
    MultiArray<3, float> in(Shape(info.width(), info.height(), info.stacksize()));

    readSIF(info, in);
    \endcode
*/
VIGRA_EXPORT void readSIF(const SIFImportInfo &info, ArrayViewND<3, float> array);

template <unsigned int N, class T>
void readSIF(const SIFImportInfo &, ArrayViewND<N, T>)
{
    vigra_precondition(false, "readSIF(): Destination array must be ArrayViewND<3, float>.");
}

inline void readSIF(const SIFImportInfo &info, ArrayViewND<3, float> array)
{
    readSIF(info, ArrayViewND<3, float>(array));
}

/**
    \brief Read parts of the image data from an Andor SIF file specified with an SIFImportInfo object
    and write them into the MultiArray array.

    \code
    SIFImportInfo info(filename);

    // create a 3D array of appropriate size
    MultiArray<3, float> in(Shape<3>(info.width(), info.height(), 1));

    readBlock(info, Shape<3>(0,0,0), Shape<3>(w,h,1), im); // read the first frame only

    \endcode
*/
VIGRA_EXPORT void readSIFBlock(const SIFImportInfo &info, Shape<3> offset, Shape<3> shape, ArrayViewND<3, float> array);

template <unsigned int N, class T>
void readSIFBlock(const SIFImportInfo &, Shape<3>, Shape<3>, ArrayViewND<N, T>)
{
    vigra_precondition(false, "readSIFBlock(): Destination array must be ArrayViewND<3, float>.");
}

inline void readSIFBlock(const SIFImportInfo &info, Shape<3> offset, Shape<3> shape, ArrayViewND<3, float> array)
{
    readSIFBlock(info, offset, shape, ArrayViewND<3, float>(array));
}

VIGRA_EXPORT std::ostream& operator<<(std::ostream& os, const SIFImportInfo& info);

//@}

} // namespace vigra

#endif // VIGRA_SIFIMPORT_HXX
