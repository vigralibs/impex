/************************************************************************/
/*                                                                      */
/*           Copyright 2004-2012 by Ullrich Koethe                      */
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

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <vigra2/imageio/impex.hxx>
#include <vigra2/imageio/impexalpha.hxx>
#include <vigra2/unittest.hxx>
#include <vigra2/array_nd.hxx>

// FIXME: is tiff.hxx still necessary?
#if 0
#if HasTIFF
# include <vigra2/imageio/tiff.hxx>
#endif
#endif

using namespace vigra;

template <class Image>
void failCodec(Image const & img, ImageExportInfo const & info)
{
    try {
        exportImage (img, info);
        failTest( "Failed to throw exception." );
    }
    catch( vigra::ContractViolation & e )
    {
        std::string expected = "\nPrecondition violation!\n";
        expected += "did not find a matching codec for the given file extension";
        const bool rc = std::strncmp( expected.c_str(), e.what(), expected.length() ) == 0;
        should(rc);
    }
}

class ByteImageExportImportTest
{
    typedef vigra::ArrayND<2, uint8_t> Image;
    typedef vigra::ArrayViewND<2, uint8_t> View;

public:

    ByteImageExportImportTest ()
    {
        vigra::ImageImportInfo info ("lenna.xv");

        const int w = info.width ();
        const int h = info.height ();

        img.resize({h,w});

        importImage (info, img);
    }

    void testListFormatsExtensions()
    {
        const std::string formats = impexListFormats();
        const std::string extensions = impexListExtensions();

        const char * refFormats = "BMP "
#if defined(HasEXR)
        "EXR "
#endif
        "GIF HDR "
#if defined(HasJPEG)
        "JPEG "
#endif
#if defined(HasPNG)
        "PNG "
#endif
        "PNM SUN "
#if defined(HasTIFF)
        "TIFF "
#endif
        "VIFF";
        shouldEqual(formats, refFormats);

        const char * refExtensions = "bmp "
#if defined(HasEXR)
        "exr "
#endif
        "gif hdr "
#if defined(HasJPEG)
        "jpeg jpg "
#endif
        "pbm pgm "
#if defined(HasPNG)
        "png "
#endif
        "pnm ppm ras "
#if defined(HasTIFF)
        "tif tiff "
#endif
        "xv";
        shouldEqual(extensions, refExtensions);
    }

    void testIsImage()
    {
        should(isImage("lenna.xv"));
        should(!isImage("no-image.txt"));
        should(!isImage("filename-does-not-exist.gif"));
    }

    void testFile (const char *filename);

    void testGIF ()
    {
        exportImage (View(img), "res.gif");

        vigra::ImageImportInfo info ("res.gif");

        should (info.width () == img.shape(1));
        should (info.height () == img.shape(0));
        should (info.isGrayscale ());
        should (info.pixelType () == vigra::ImageImportInfo::UINT8);

        Image res (info.height (), info.width ());

        importImage (info, View (res));

        auto i = img.begin ();
        auto i1 = res.begin ();

        float sum = 0;
        for (; i != img.end (); ++i, ++i1)
            sum += abs (*i - *i1);

        should (sum / (info.width () * info.height ()) < 0.1);

        ArrayND<2, uint8_t> res1;
        importImage("res.gif", res1);
        should(res1 == View(res));
    }

    void testGrayToRGB()
    {
        ArrayND<2, TinyArray<uint8_t,3> > rgb;

        importImage("lenna.xv", rgb);

        should (rgb.shape(0) == img.shape(0));
        should (rgb.shape(1) == img.shape(1));
        shouldEqualSequence(img.begin(), img.end(), rgb.bindChannel(0).begin());
        shouldEqualSequence(img.begin(), img.end(), rgb.bindChannel(1).begin());
        shouldEqualSequence(img.begin(), img.end(), rgb.bindChannel(2).begin());
    }

    void testJPEG ()
    {
        vigra::ImageExportInfo exportinfo ("res.jpg");
#if !defined(HasJPEG)
        failCodec(img, exportinfo);
#else
        exportinfo.setCompression ("JPEG QUALITY=100");
        exportImage (img, exportinfo);

        vigra::ImageImportInfo info ("res.jpg");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isGrayscale ());
        should (info.pixelType () == vigra::ImageImportInfo::UINT8);

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        float sum = 0;
        for (; i != img.end (); ++i, ++i1)
            sum += abs (*i - *i1);

        should (sum / (info.width () * info.height ()) < 0.1);
#endif
    }

    void testTIFF ()
    {
        vigra::ImageExportInfo exportinfo ("res.tif");
#if !defined(HasTIFF)
        failCodec(img, exportinfo);
#else
        exportinfo.setCompression ("LZW");
        exportImage (img, exportinfo);

        vigra::ImageImportInfo info ("res.tif");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.height (), info.width ());

        importImage (info, res);
        shouldEqual(res.axistags(0),tags::axis_y);
        shouldEqual(res.axistags(1),tags::axis_x);

        auto i = img.begin ();
        auto i1 = res.begin ();

        for (; i != img.end (); ++i, ++i1)
            should (*i == *i1);
#if 0
        TiffImage * tiff = TIFFOpen("res2.tif", "w");
        createTiffImage(View(img), tiff);
        TIFFClose(tiff);

        uint32 w, h;
        tiff = TIFFOpen("res2.tif", "r");
        TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
        TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);
        shouldEqual(w, img.shape (1));
        shouldEqual(h, img.shape (0));

        MultiArray<2, unsigned char> res2(w,h);
        importTiffImage(tiff, res2);
        TIFFClose(tiff);

        shouldEqualSequence(res2.begin(), res2.end(), img.data());

        // test bilevel
        MultiArray<2, unsigned char> bilevel;
        importImage("bilevel.tiff", bilevel);
        shouldEqual(bilevel.shape(), Shape2(1577, 1083));
        UInt8 m, M;
        bilevel.minmax(&m, &M);
        shouldEqual(m, 0);
        shouldEqual(M, 1);
        shouldEqual(bilevel.sum<int>(), 1653050); // 96% white pixels
#endif

#endif
    }

    void testTIFFSequence()
    {
#if defined(HasTIFF)
        for (int i=0; i < 3; ++i)
        {
            std::string fileName = std::string("lenna_") + vigra::asString(i) + ".tif";
            vigra::ImageImportInfo ininfo (fileName.c_str());
            Image inimg(ininfo.width(), ininfo.height());
            importImage(ininfo, inimg);
            vigra::ImageExportInfo outinfo ("resseq.tif", i==0?"w":"a");
            exportImage(inimg, outinfo);
        }

        for (int j=0; j < 3; ++j)
        {
            vigra::ImageImportInfo ininfo ("resseq.tif", j);
            Image inimg(ininfo.width(), ininfo.height());
            std::string fileName = std::string("lenna_") + vigra::asString(j) + ".tif";
            importImage(ininfo, inimg);
            vigra::ImageImportInfo originfo (fileName.c_str());
            Image origimg(originfo.width(), originfo.height());
            importImage(originfo, origimg);

            auto it = inimg.begin ();
            auto it1 = origimg.begin ();
            for (; it != inimg.end (); ++it, ++it1)
                should (*it == *it1);
        }
#endif
    }

    void testBMP ()
    {
        testFile ("res.bmp");
    }

    void testPGM ()
    {
        testFile ("res.pgm");
    }

    void testPNM ()
    {
        testFile ("res.pnm");
    }

    void testPNM2 ()
    {
        vigra::ImageExportInfo exportinfo ("res.pgm");
        exportinfo.setCompression ("ASCII");
        exportImage (img, exportinfo);

        vigra::ImageImportInfo info ("res.pgm");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));
        should (info.getFileType () == std::string ("PNM"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        for (; i != img.end (); ++i, ++i1)
            should (*i == *i1);
    }

    void testPNG ()
    {
#if !defined(HasPNG)
        failCodec(img, vigra::ImageExportInfo("res.png"));
#else
        testFile ("res.png");
#endif
    }

    void testSUN ()
    {
        testFile ("res.ras");
    }

    void testVIFF1 ()
    {
        testFile ("res.xv");
    }

    void testVIFF2 ()
    {
        vigra::ImageExportInfo exportinfo ("res.foo");
        exportinfo.setFileType ("VIFF");
        exportImage (img, exportinfo);

        vigra::ImageImportInfo info ("res.foo");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));
        should (info.getFileType () == std::string ("VIFF"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        for (; i != img.end (); ++i, ++i1)
            should (*i == *i1);
    }

    Image img;
};

void
ByteImageExportImportTest::testFile (const char *filename)
{
    exportImage (img, vigra::ImageExportInfo (filename));

    vigra::ImageImportInfo info (filename);

    should (info.width () == img.shape (1));
    should (info.height () == img.shape (0));
    should (info.isGrayscale ());
    should (info.getPixelType () == std::string ("UINT8"));

    Image res (info.height (), info.width ());

    importImage (info, res);

    auto i = img.begin ();
    auto i1 = res.begin ();

    for (; i.isValid(); ++i, ++i1)
        should (*i == *i1);
}


class ByteRGBImageExportImportTest
{
    typedef ArrayND<2, TinyArray<uint8_t,3> > Image;
public:
    ByteRGBImageExportImportTest ()
    {
        vigra::ImageImportInfo info ("lennargb.xv");

        int w = info.width ();
        int h = info.height ();

        img.resize ({h, w});

        importImage (info, img);
    }

    void testFile (const char *fileName);

    void testGIF ()
    {
        exportImage (img, vigra::ImageExportInfo ("resrgb.gif"));

        vigra::ImageImportInfo info ("resrgb.gif");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isColor ());
        should (info.pixelType () == vigra::ImageImportInfo::UINT8);

        Image res (info.height (), info.width ());
        Image ref (info.height (), info.width ());

        importImage (info, res);
        importImage (vigra::ImageImportInfo("lenna_gifref.xv"), ref);
        shouldEqual(res.axistags(), AxisTags<2>(tags::axis_y, tags::axis_x));

        auto i = ref.begin ();
        auto i1 = res.begin ();

        double sum = 0.0;
        for (; i != ref.end (); ++i, ++i1)
                sum += norm(*i - *i1);

        should (sum / (info.width () * info.height ()) < 4.0);  // use rather large tolerance to make the
                                                                // test portable
        ref.expandElements(2) = 0;
        importImage(info, ref.expandElements(2));
        should(res == ref);

        ArrayND<3, uint8_t> banded;
        importImage(info, banded);
        shouldEqual(banded.shape(), info.shape().insert(2, 3));
        shouldEqual(banded.axistags(), AxisTags<3>(tags::axis_y, tags::axis_x, tags::axis_c));
        should(banded == ref.expandElements(2));

        importImage(info, banded, F_ORDER);
        shouldEqual(banded.shape(), reversed(info.shape().insert(2, 3)));
        shouldEqual(banded.axistags(), AxisTags<3>(tags::axis_c, tags::axis_x, tags::axis_y));
        should(transpose(banded) == ref.expandElements(2));

        Image noncontiguous1(reversed(info.shape()), AxisTags<2>(tags::axis_x, tags::axis_y));
        importImage(info, noncontiguous1.view());
        should(transpose(noncontiguous1) == res);

        auto view1 = noncontiguous1.expandElements(0);
        view1 = 0;
        importImage(info, view1);
        should(transpose(view1) == res.expandElements(2));

        Image noncontiguous2(reversed(info.shape()), F_ORDER);
        importImage(info, noncontiguous2.view());
        should(transpose(noncontiguous2) == res);

        auto view2 = noncontiguous2.expandElements(0);
        view2 = 0;
        importImage(info, view2);
        should(transpose(view2) == res.expandElements(2));
    }

    void testJPEG ()
    {
#if !defined(HasJPEG)
        failCodec(img, vigra::ImageExportInfo ("res.jpg").setCompression ("JPEG QUALITY=100"));
#else
        exportImage (img,
                     vigra::ImageExportInfo ("res.jpg").setCompression ("JPEG QUALITY=100"));

        vigra::ImageImportInfo info ("res.jpg");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isColor ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        double sum = 0.0;
        for (; i != img.end (); ++i, ++i1)
            {
                sum += norm(*i - *i1);
            }
        should (sum / (info.width () * info.height ()) < 2.0);
#endif
    }

    void testTIFF ()
    {
#if !defined(HasTIFF)
        failCodec(img, vigra::ImageExportInfo ("res.tif").setCompression ("LZW"));
#else
        exportImage (img,
                     vigra::ImageExportInfo ("res.tif").
                     setCompression ("LZW"));

        vigra::ImageImportInfo info ("res.tif");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isColor ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        for (; i != img.end (); ++i, ++i1)
            {
                should (*i == *i1);
            }
#if 0
        TiffImage * tiff = TIFFOpen("res2.tif", "w");
        createTiffImage(MultiArrayView<2, RGBValue<unsigned char> >(img), tiff);
        TIFFClose(tiff);

        uint32 w, h;
        tiff = TIFFOpen("res2.tif", "r");
        TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &w);
        TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &h);
        shouldEqual(w, img.shape (1));
        shouldEqual(h, img.shape (0));

        MultiArray<2, RGBValue<unsigned char> > res2(w,h);
        importTiffImage(tiff, res2);
        TIFFClose(tiff);

        shouldEqualSequence(res2.begin(), res2.end(), img.data());
#endif
#endif
    }

    void testBMP ()
    {
        testFile ("res.bmp");
    }

    void testPPM ()
    {
        testFile ("res.ppm");
    }

    void testPNM ()
    {
        testFile ("res.pnm");
    }


    void testPNM2 ()
    {
        exportImage (img,
                     vigra::ImageExportInfo ("res.ppm").setCompression ("ASCII"));

        vigra::ImageImportInfo info ("res.ppm");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isColor ());
        should (info.getPixelType () == std::string ("UINT8"));
        should (info.getFileType () == std::string ("PNM"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        for (; i != img.end (); ++i, ++i1)
            should (*i == *i1);
    }

    void testPNG ()
    {
#if !defined(HasPNG)
        failCodec(img, vigra::ImageExportInfo("res.png"));
#else
        testFile ("res.png");
#endif
    }

    void testSUN ()
    {
        testFile ("res.ras");
    }

    void testVIFF1 ()
    {
        testFile ("res.xv");
    }

    void testVIFF2 ()
    {
        exportImage (img,
                     vigra::ImageExportInfo ("res.foo").setFileType ("VIFF"));

        vigra::ImageImportInfo info ("res.foo");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isColor ());
        should (info.getPixelType () == std::string ("UINT8"));
        should (info.getFileType () == std::string ("VIFF"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        for (; i != img.end (); ++i, ++i1)
            should (*i == *i1);
    }

    Image img;
};

void
ByteRGBImageExportImportTest::testFile (const char *fileName)
{
    exportImage (img, vigra::ImageExportInfo (fileName));

    vigra::ImageImportInfo info (fileName);

    should (info.width () == img.shape (1));
    should (info.height () == img.shape (0));
    should (info.isColor ());
    should (info.getPixelType () == std::string ("UINT8"));

    Image res (info.height (), info.width ());

    importImage (info, res);

    auto i = img.begin ();
    auto i1 = res.begin ();

    for (; i != img.end (); ++i, ++i1)
        should (*i == *i1);
}

class CanvasSizeTest
{
  public:
    void testTIFFCanvasSize ()
    {
        vigra::ImageExportInfo exportinfo ("res.tif");
        ArrayND<2,TinyArray<float,3> > img(Shape<2>{1, 1});
#if !defined(HasTIFF)
        failCodec(img, exportinfo);
#else
        img(0,0) = 1;
        exportinfo.setCompression ("LZW");
        // TODO: fix axis order.
        Shape<2> canvasSize(3, 8);
        exportinfo.setCanvasSize (canvasSize);
        exportImage (img, exportinfo);

        vigra::ImageImportInfo info ("res.tif");

        should (info.getCanvasSize () == canvasSize);
#endif
    }
};

class PositionTest
{
  public:
    void testFile(const char* filename)
    {
        ImageExportInfo exportinfo(filename);
        ArrayND<2,TinyArray<float,4> > img(Shape<2>{1, 1});
        img(0, 0) = 1;

        const Shape<2> position(0, 100);
        exportinfo.setPosition(position);
        exportinfo.setXResolution(1.0);
        exportinfo.setYResolution(1.0);

        exportImage(img, exportinfo);

        ImageImportInfo info(filename);

        should (info.getPosition() == position);
    }

    void testEXRPosition()
    {
#if !defined(HasEXR)
        ArrayND<2,TinyArray<float,3> > img(Shape<2>{1, 1});
        failCodec(img, vigra::ImageExportInfo("res.exr"));
#else
        testFile("res.exr");
#endif
    }

    void testTIFFPosition()
    {
#if !defined(HasTIFF)
        ArrayND<2,TinyArray<float,3> > img(Shape<2>{1, 1});
        failCodec(img, vigra::ImageExportInfo("res.tif"));
#else
        testFile("res.tif");
#endif
    }

    void testPNGPosition()
    {
#if !defined(HasPNG)
        ArrayND<2,TinyArray<float,3> > img(Shape<2>{1, 1});
        failCodec(img, vigra::ImageExportInfo("res.png"));
#else
        testFile("res.png");
#endif
    }
};

class PNGInt16Test
{
  public:
    void testByteOrder()
    {
        ArrayND<2,uint16_t> i(1,1);
        i(0,0) = 1;
        exportImage(i, ImageExportInfo("res.png"));
        ImageImportInfo info("res.png");
        shouldEqual(info.width(), 1);
        shouldEqual(info.height(), 1);
        shouldEqual(info.numBands(), 1);
        shouldEqual(info.isGrayscale(), true);
        shouldEqual(std::string(info.getPixelType()), std::string("UINT16"));
        i(0,0) = 0;
        importImage(info, i);
        shouldEqual(i(0,0), 1);

        // DGSW: Note that this produces a PNG standard conformant image
        //       but both Imagemagick 'identify' and photoshop CS2 see
        //       the data incorrectly
        ArrayND<2, TinyArray<uint16_t,3> > rgb(1,1);
        // Using unsigned values 0xff01, 0xfff1, 0xfffd
        rgb(0,0) = {65281,65521,65533};
        // Using unsigned values 0x7f01, 0x7ff1, 0x7ffd
        // rgb(0,0) = RGBValue<unsigned short>(32513,32753,32765);
        exportImage(rgb, ImageExportInfo("res.png"));
        ImageImportInfo rgbinfo("res.png");
        shouldEqual(rgbinfo.width(), 1);
        shouldEqual(rgbinfo.height(), 1);
        shouldEqual(rgbinfo.numBands(), 3);
        shouldEqual(std::string(rgbinfo.getPixelType()), std::string("UINT16"));
        rgb(0,0) = {0,0,0};
        importImage(rgbinfo, rgb);
        shouldEqual(rgb(0,0), (TinyArray<uint16_t,3>(65281u,65521u,65533u)));
//        shouldEqual(rgb(0,0), RGBValue<unsigned short>(32513,32753,32765));
    }
};

class FloatImageExportImportTest
{
    typedef ArrayND<2, double> Image;
    std::string rereadType;

public:

    FloatImageExportImportTest ()
    : rereadType("DOUBLE")
    {
        vigra::ImageImportInfo info ("lenna.xv");

        int w = info.width ();
        int h = info.height ();

        img.resize ({h, w});

        importImage (info, img);

        vigra::ImageImportInfo rinfo ("lennafloat.xv");

        reread.resize ({h,w});

        importImage (rinfo, reread);
    }

    void testGIF()
    {
        exportImage (img, vigra::ImageExportInfo ("res.gif"));

        vigra::ImageImportInfo info ("res.gif");

        should (info.width () == reread.shape (1));
        should (info.height () == reread.shape (0));
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = reread.begin ();
        auto i1 = res.begin ();

        double sum = 0.0;
        for (; i != reread.end (); ++i, ++i1)
            sum += abs (*i - *i1);
        should (sum / (info.width () * info.height ()) < 0.1);
    }

    void testJPEG ()
    {
#if !defined(HasJPEG)
        failCodec(img, vigra::ImageExportInfo ("res.jpg").setCompression ("JPEG QUALITY=100"));
#else
        exportImage (img,
                     vigra::ImageExportInfo ("res.jpg").setCompression ("JPEG QUALITY=100"));

        vigra::ImageImportInfo info ("res.jpg");

        should (info.width () == reread.shape (1));
        should (info.height () == reread.shape (0));
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = reread.begin ();
        auto i1 = res.begin ();

        double sum = 0.0;
        for (; i != reread.end (); ++i, ++i1)
            sum += abs (*i - *i1);
        should (sum / (info.width () * info.height ()) < 0.1);
#endif
    }

    void testPNG ()
    {
#if !defined(HasPNG)
        failCodec(img, vigra::ImageExportInfo ("res.png"));
#else
        exportImage (img,
                     vigra::ImageExportInfo ("res_.png"));

        vigra::ImageImportInfo info ("res_.png");

        should (info.width () == reread.shape (1));
        should (info.height () == reread.shape (0));
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = reread.begin ();
        auto i1 = res.begin ();

        double sum = 0.0;
        for (; i != reread.end (); ++i, ++i1)
            sum += abs (*i - *i1);
        should (sum / (info.width () * info.height ()) < 0.1);
#endif
    }

    void testTIFF ()
    {
#if !defined(HasTIFF)
        failCodec(img, vigra::ImageExportInfo ("res.tif").setCompression ("LZW"));
#else
        exportImage (img,
                     vigra::ImageExportInfo ("res.tif").setCompression ("LZW"));

        vigra::ImageImportInfo info ("res.tif");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isGrayscale ());
        should (info.getPixelType () == rereadType);

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        shouldEqualSequence(i, img.end(), i1);
#endif
    }

    void testTIFFForcedRange ()
    {
#if !defined(HasTIFF)
        failCodec(img, vigra::ImageExportInfo ("res.tif").setForcedRangeMapping(0, 255, 1, 2));
#else
        exportImage (img,
                     vigra::ImageExportInfo ("res.tif").setForcedRangeMapping(0, 255, 1, 2));

        vigra::ImageImportInfo info ("res.tif");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isGrayscale ());
        should (info.getPixelType () == rereadType);

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        for (; i != img.end (); ++i, ++i1)
            shouldEqualTolerance(*i / 255.0, *i1 - 1.0, 1e-12);
#endif
    }

    void testBMP ()
    {
        exportImage (img, vigra::ImageExportInfo ("res.bmp"));

        vigra::ImageImportInfo info ("res.bmp");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = reread.begin ();
        auto i1 = res.begin ();

        for (; i != reread.end (); ++i, ++i1)
            should (*i == *i1);
    }

    void testSUN ()
    {
        exportImage (img, vigra::ImageExportInfo ("res.ras"));

        vigra::ImageImportInfo info ("res.ras");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isGrayscale ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = reread.begin ();
        auto i1 = res.begin ();

        for (; i != reread.end (); ++i, ++i1)
            {
                should (*i == *i1);
            }
    }

    void testVIFF ()
    {
        exportImage (img, vigra::ImageExportInfo ("res.xv"));

        vigra::ImageImportInfo info ("res.xv");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isGrayscale ());
        should (info.getPixelType () == rereadType);

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        for (; i != img.end (); ++i, ++i1)
            {
                should (*i == *i1);
            }
    }

    Image img, reread;
};

class FloatRGBImageExportImportTest
{
    typedef ArrayND<2,TinyArray<float,3> > Image;
    Image img, reread;

public:

    FloatRGBImageExportImportTest ()
    {
        vigra::ImageImportInfo info ("lennargb.xv");

        int w = info.width ();
        int h = info.height ();

        img.resize ({h,w});

        importImage (info, img);

        vigra::ImageImportInfo rinfo ("lennafloatrgb.xv");

        reread.resize ({h,w});

        importImage (rinfo, reread);
    }

    void testJPEG ()
    {
#if !defined(HasJPEG)
        failCodec(img, vigra::ImageExportInfo ("res.jpg").setCompression ("JPEG QUALITY=100"));
#else
        exportImage (img,
                     vigra::ImageExportInfo ("res.jpg").setCompression ("JPEG QUALITY=100"));

        vigra::ImageImportInfo info ("res.jpg");

        should (info.width () == reread.shape (1));
        should (info.height () == reread.shape (0));
        should (info.isColor ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = reread.begin ();
        auto i1 = res.begin ();

        float sum = 0.0f;
        for (; i != reread.end (); ++i, ++i1)
            sum += norm(*i - *i1);
        should (sum / (info.width () * info.height ()) < 2.0f);
#endif
    }

    void testTIFF ()
    {
#if !defined(HasTIFF)
        failCodec(img, vigra::ImageExportInfo ("res.tif"));
#else
        exportImage (img,
                     vigra::ImageExportInfo ("res.tif"));

        vigra::ImageImportInfo info ("res.tif");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isColor ());
        should (info.getPixelType () == std::string ("FLOAT"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        for (; i != img.end (); ++i, ++i1)
            shouldEqual (*i, *i1);
#endif
    }

    void testTIFFForcedRange ()
    {
#if !defined(HasTIFF)
        failCodec(img, vigra::ImageExportInfo ("res.tif").setForcedRangeMapping(0, 255, 1, 2));
#else
        exportImage (img,
                     vigra::ImageExportInfo ("res.tif").setForcedRangeMapping(0,255,1,2));

        vigra::ImageImportInfo info ("res.tif");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isColor ());
        should (info.getPixelType () == std::string ("FLOAT"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        for (; i != img.end (); ++i, ++i1)
        {
            should(closeAtTolerance(*i/255.0f, *i1-1.0f, 1e-4));
            // shouldEqualTolerance(acc.red(i)/255.0f, (*i1)[0]-1.0f, 1e-4);
            // shouldEqualTolerance(acc.green(i)/255.0f, (*i1)[1]-1.0f, 1e-4);
            // shouldEqualTolerance(acc.blue(i)/255.0f, (*i1)[2]-1.0f, 1e-4);
        }
#endif
    }

    void testBMP ()
    {
        exportImage (img, vigra::ImageExportInfo ("res.bmp"));

        vigra::ImageImportInfo info ("res.bmp");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isColor ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = reread.begin ();
        auto i1 = res.begin ();

        float sum = 0.0f;
        for (; i != reread.end (); ++i, ++i1)
            sum += norm(*i - *i1);
        should (sum / (info.width () * info.height ()) < 2.0f);
    }

    void testSUN ()
    {
        exportImage (img, vigra::ImageExportInfo ("res.ras"));

        vigra::ImageImportInfo info ("res.ras");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isColor ());
        should (info.getPixelType () == std::string ("UINT8"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = reread.begin ();
        auto i1 = res.begin ();

        float sum = 0.0f;
        for (; i != reread.end (); ++i, ++i1)
            sum += norm(*i - *i1);
        should (sum / (info.width () * info.height ()) < 2.0f);
    }

    void testVIFF ()
    {
        exportImage (img, vigra::ImageExportInfo ("res.xv"));

        vigra::ImageImportInfo info ("res.xv");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isColor ());
        should (info.getPixelType () == std::string ("FLOAT"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        for (; i != img.end (); ++i, ++i1)
            should (*i == *i1);
    }

    void testHDR ()
    {
        vigra::ImageExportInfo exi("res.hdr");

        exportImage (img, exi) ;

        vigra::ImageImportInfo info ("res.hdr");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isColor ());
        should (info.getPixelType () == std::string ("FLOAT"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        for (; i != img.end (); ++i, ++i1)
            shouldEqual (*i, *i1);
    }

};

class Vector4ExportImportTest
{
  public:

    typedef ArrayND<2, TinyArray<float,4>> Image;
    typedef ArrayND<2, TinyArray<uint8_t,4>> BImage;
    Image img, reread;
    BImage breread, breference;

public:

    Vector4ExportImportTest ()
    : img(2,3),
      reread(2,3),
      breread(2,3),
      breference(2,3)
    {
        double scale = 255.0 / 11.0;
        double offset = 5.5;
        for(int y = 0; y<3; ++y)
        {
            for(int x=0; x<2; ++x)
            {
                img(x,y)[0] = 2*y+x + 0.5f;
                img(x,y)[1] = -img(x,y)[0];
                img(x,y)[2] = 0.0;
                img(x,y)[3] = 0.5;
                for(int b=0; b<4; ++b)
                {
                    breference(x,y)[b] =
                        NumericTraits<uint8_t>::fromRealPromote(scale*(img(x,y)[b]+offset));
                }
            }
        }
    }

    void failingTest (char const * filename,
                      char const * message = "exportImage(): file format does not support requested number of bands (color channels)")
    {
        try
        {
            exportImage( img, vigra::ImageExportInfo( filename ) );
            failTest( "Failed to throw exception." );
        }
        catch( vigra::ContractViolation & e )
        {
            std::string expected = "\nPrecondition violation!\n";
            expected += message;
            const bool rc = std::strncmp( expected.c_str(), e.what(), expected.length() ) == 0;
            should(rc);
        }
    }

    void testJPEG ()
    {
#if !defined(HasJPEG)
        failingTest("res.jpg", "did not find a matching codec for the given file extension");
#else
        failingTest("res.jpg");
#endif
    }

    void testGIF ()
    {
        failingTest("res.gif");
    }

    void testBMP ()
    {
        failingTest("res.bmp");
    }

    void testPNM ()
    {
        failingTest("res.pnm");
    }

    void testSUN ()
    {
        failingTest("res.ras");
    }

    void testVIFF ()
    {
        exportImage (img,
                     vigra::ImageExportInfo ("res.xv"));

        vigra::ImageImportInfo info ("res.xv");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.numBands () == 4);
        should (info.getPixelType () == std::string ("FLOAT"));

        importImage (info, reread);

        auto i = img.begin ();
        auto i1 = reread.begin ();

        for (; i != img.end (); ++i, ++i1)
            shouldEqual (*i, *i1);
    }

    void testTIFF ()
    {
#if !defined(HasTIFF)
        failCodec(img, vigra::ImageExportInfo ("res.tif"));
#else
        exportImage (img,
                     vigra::ImageExportInfo ("res.tif"));

        vigra::ImageImportInfo info ("res.tif");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.numBands () == 4);
        should (info.getPixelType () == std::string ("FLOAT"));

        importImage (info, reread);

        auto i = img.begin ();
        auto i1 = reread.begin ();

        for (; i != img.end (); ++i, ++i1)
            shouldEqual (*i, *i1);
#endif
    }

    void testEXR ()
    {
#if !defined(HasEXR)
        failCodec(img, vigra::ImageExportInfo ("res.exr"));
#else
        exportImage (img,
                     vigra::ImageExportInfo ("res.exr"));

        vigra::ImageImportInfo info ("res.exr");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.isColor ());
        should (info.getPixelType () == std::string ("FLOAT"));

        Image res (info.height (), info.width ());

        importImage (info, res);

        auto i = img.begin ();
        auto i1 = res.begin ();

        for (; i != img.end (); ++i, ++i1)
            shouldEqual (*i, *i1);
#endif
    }

    void testPNG ()
    {
#if !defined(HasPNG)
        failCodec(img, vigra::ImageExportInfo("res.png"));
#else
        exportImage (img,
                     vigra::ImageExportInfo ("res.png"));

        vigra::ImageImportInfo info ("res.png");

        should (info.width () == img.shape (1));
        should (info.height () == img.shape (0));
        should (info.numBands () == 4);
        should (info.getPixelType () == std::string ("UINT8"));

        importImage (info, breread);

        auto i = breference.begin ();
        auto i1 = breread.begin ();

        for (; i != breference.end (); ++i, ++i1)
        {
            should (norm(*i- *i1) <= 1.0);
        }
#endif
    }
};

class ImageExportImportFailureTest
{
    ArrayND<2, uint8_t> img;

public:

    ImageExportImportFailureTest()
        : img( 3, 3 )
    {}

    // gif

    void testGIFExport()
    {
        testExport("gif");
    }

    void testGIFImport()
    {
        testImport("gif");
    }

    // jpeg

    void testJPEGExport()
    {
#if !defined(HasJPEG)
        testExport("jpg", "did not find a matching codec for the given file extension");
#else
        testExport("jpg");
#endif
    }

    void testJPEGImport()
    {
        testImport("jpg");
    }

    // tiff

    void testTIFFExport()
    {
#if !defined(HasTIFF)
        testExport("tiff", "did not find a matching codec for the given file extension");
#else
        testExport("tiff");
#endif
    }

    void testTIFFImport()
    {
        testImport("tiff");
    }

    // exr

    void testEXRExport()
    {
#if !defined(HasEXR)
        testExport("exr", "did not find a matching codec for the given file extension");
#else
        testExport("exr");
#endif
    }

    // viff

    void testVIFFExport()
    {
        testExport("xv");
    }

    void testVIFFImport()
    {
        testImport("xv");
    }

    // sun

    void testSUNExport()
    {
        testExport("ras");
    }

    void testSUNImport()
    {
        testImport("ras");
    }

    // pnm

    void testPNMExport()
    {
        testExport("pnm");
    }

    void testPNMImport()
    {
        testImport("pnm");
    }

    // png

    void testPNGExport()
    {
#if !defined(HasPNG)
        testExport("png", "did not find a matching codec for the given file extension");
#else
        testExport("png");
#endif
    }

    void testPNGImport()
    {
        testImport("png");
    }

    // bmp

    void testBMPExport()
    {
        testExport("bmp");
    }

    void testBMPImport()
    {
        testImport("bmp");
    }

    // test implementation

    void testImport( const char * fext )
    {
        std::string fname = "foo.";
        fname += fext;
        try {
            vigra::ImageImportInfo info( fname.c_str() );
            failTest( "Failed to throw exception." );
        }
        catch( vigra::ContractViolation & e ) {
            std::string expected = "\nPrecondition violation!\n";
            expected += "Unable to open file '";
            expected += fname;
            expected += "'.";
            const bool rc = std::strncmp( expected.c_str(), e.what(), expected.length() ) == 0;
            should(rc);
        }
    }

    void testExport( const char * fext ,
                     const char * message = 0)
    {
        std::string fname = "intentionalFailure/foo.";
        fname += fext;
        try {
            exportImage( img, vigra::ImageExportInfo( fname.c_str() ) );
            failTest( "Failed to throw exception." );
        }
        catch( vigra::ContractViolation & e ) {
            std::string expected = "\nPrecondition violation!\n";
            if(message)
            {
                expected += message;
            }
            else
            {
                expected += "Unable to open file '";
                expected += fname;
                expected += "'.";
            }
            const bool rc = std::strncmp( expected.c_str(), e.what(), expected.length() ) == 0;
            should(rc);
        }
    }

    void testShapeMismatch ()
    {
        ArrayND<2, TinyArray<uint8_t,3> > rgb(1,1);

        try {
            importImage(ImageImportInfo("lennargb.xv"), rgb.view());
            failTest( "Failed to throw exception." );
        }
        catch( vigra::ContractViolation & e ) {
            std::string expected = "\nPrecondition violation!\nimportImage(): shape mismatch between input and output.";
            const bool rc = std::strncmp( expected.c_str(), e.what(), expected.length() ) == 0;
            should(rc);
        }

        ArrayND<2, TinyArray<uint8_t,4> > vec4;

        try {
            importImage("lennargb.xv", vec4);
            failTest( "Failed to throw exception." );
        }
        catch( vigra::ContractViolation & e ) {
            std::string expected = "\nPrecondition violation!\nimportImage(): Number of channels in input and destination image don't match.";
            const bool rc = std::strncmp( expected.c_str(), e.what(), expected.length() ) == 0;
            should(rc);
        }
    }
};
#if 0
class GrayscaleImportExportAlphaTest
{
public:
    GrayscaleImportExportAlphaTest()
    {
#if defined(HasTIFF)
        vigra::ImageImportInfo info("lenna_masked_gray.tif");

        image_.resize(info.size());
        alpha_.resize(info.size());
        importImageAlpha(info, image_, alpha_);
#else
        image_.resize(Shape<2>(20,10));
        alpha_.resize(image_.shape());

        image_.init(10);
        alpha_.init(255);
#endif
    }

    void testFile(const char* filename);
    void testFileMultiArray(const char* filename);

    void testTIFF()
    {
        const char filename[] = "res.tif";

#if defined(HasTIFF)
        testFile(filename);
        testFileMultiArray(filename);
#else
        failCodec(image_, vigra::ImageExportInfo(filename));
#endif
    }

    void testPNG()
    {
        const char filename[] = "res.png";

#if defined(HasPNG)
        testFile(filename);
        testFileMultiArray(filename);
#else
        failCodec(image_, vigra::ImageExportInfo(filename));
#endif
    }

private:
    ArrayND<2, uint8_t> image_;
    ArrayND<2, uint8_t> alpha_;
};

void
GrayscaleImportExportAlphaTest::testFile(const char* filename)
{
    exportImageAlpha(image_, alpha_, vigra::ImageExportInfo(filename));

    vigra::ImageImportInfo info(filename);

    should(info.width() == image_.shape(1));
    should(info.height() == image_.shape(0));
    should(!info.isColor());
    should(!strcmp(info.getPixelType(), "UINT8"));
    should(info.numBands() == 2);

    should(info.width() == alpha_.shape(1));
    should(info.height() == alpha_.shape(0));
    should(info.numExtraBands() == 1);

    ArrayND<2, uint8_t> image(info.size());
    ArrayND<2, uint8_t> alpha(info.size());

    importImageAlpha(info, image, alpha);

    for (auto x = alpha_.begin(), xx = alpha_.begin(); x != alpha_.end(); ++x, ++xx)
    {
        should(*x == 255);
        should(*x == *xx);
    }

    for (auto x = image_.begin(), xx = image.begin(); x != image_.end(); ++x, ++xx)
    {
        should(*x == *xx);
    }
}

void
GrayscaleImportExportAlphaTest::testFileMultiArray(const char* filename)
{
    typedef ArrayViewND<2, uint8_t> View;

    exportImageAlpha(View(image_), View(alpha_), vigra::ImageExportInfo(filename));

    vigra::ImageImportInfo info(filename);
    ArrayND<2, uint8_t> image(info.shape()),
                                 alpha(info.shape());

    importImageAlpha(info, image, alpha);

    auto xx = alpha.begin();
    for (auto x = alpha_.begin(); x != alpha_.end(); ++x, ++xx)
    {
        should(*x == 255);
        should(*x == *xx);
    }

    xx = image.begin();
    for (auto x = image_.begin(); x != image_.end(); ++x, ++xx)
    {
        should(*x == *xx);
    }
}

class RGBImportExportAlphaTest
{
public:
    RGBImportExportAlphaTest()
    {
#if defined(HasTIFF)
        vigra::ImageImportInfo info("lenna_masked_color.tif");

        image_.resize(info.size());
        alpha_.resize(info.size());
        importImageAlpha(info, image_, alpha_);
#else
        image_.resize(Shape<2>(20,10));
        alpha_.resize(image_.shape());

        image_.init(TinyArray<uint8_t, 3>(10,20,30));
        alpha_.init(255);
#endif
    }

    void testFile(const char* filename);

    void testTIFF()
    {
        const char filename[] = "res.tif";

#if defined(HasTIFF)
        testFile(filename);
#else
        failCodec(image_, vigra::ImageExportInfo(filename));
#endif
    }

    void testPNG()
    {
        const char filename[] = "res.png";

#if defined(HasPNG)
        testFile(filename);
#else
        failCodec(image_, vigra::ImageExportInfo(filename));
#endif
    }

private:
    ArrayND<2, TinyArray<uint8_t,3>> image_;
    ArrayND<2, uint8_t> alpha_;
};

void
RGBImportExportAlphaTest::testFile(const char* filename)
{
    exportImageAlpha(image_, alpha_, vigra::ImageExportInfo(filename));

    vigra::ImageImportInfo info(filename);

    should(info.width() == image_.shape(1));
    should(info.height() == image_.shape(0));
    should(info.isColor());
    should(!strcmp(info.getPixelType(), "UINT8"));
    should(info.numBands() == 4);

    should(info.width() == alpha_.shape(1));
    should(info.height() == alpha_.shape(0));
    should(info.numExtraBands() == 1);

    ArrayND<2, TinyArray<uint8_t,3>> image(info.size());
    ArrayND<2, uint8_t> alpha(info.size());

    importImageAlpha(info, image, alpha);

    for (auto x = alpha_.begin(), xx = alpha_.begin(); x != alpha_.end(); ++x, ++xx)
    {
        should(*x == 255);
        should(*x == *xx);
    }

    for (auto x = image_.begin(), xx = image.begin(); x != image_.end(); ++x, ++xx)
    {
        should(*x == *xx);
    }
}
#endif
struct ImageImportExportTestSuite : public vigra::test_suite
{
    ImageImportExportTestSuite()
        : vigra::test_suite("ImageImportExportTestSuite")
    {
        // general tests
        add(testCase(&ByteImageExportImportTest::testListFormatsExtensions));
        add(testCase(&ByteImageExportImportTest::testIsImage));

        // grayscale byte images
        add(testCase(&ByteImageExportImportTest::testGIF));
        add(testCase(&ByteImageExportImportTest::testJPEG));
        add(testCase(&ByteImageExportImportTest::testTIFF));
        add(testCase(&ByteImageExportImportTest::testTIFFSequence));
        add(testCase(&ByteImageExportImportTest::testBMP));
        add(testCase(&ByteImageExportImportTest::testPGM));
        add(testCase(&ByteImageExportImportTest::testPNM));
        add(testCase(&ByteImageExportImportTest::testPNM2));
        add(testCase(&ByteImageExportImportTest::testPNG));
        add(testCase(&ByteImageExportImportTest::testSUN));
        add(testCase(&ByteImageExportImportTest::testVIFF1));
        add(testCase(&ByteImageExportImportTest::testVIFF2));
        add(testCase(&ByteImageExportImportTest::testGrayToRGB));

        // rgb byte images
        add(testCase(&ByteRGBImageExportImportTest::testGIF));
        add(testCase(&ByteRGBImageExportImportTest::testJPEG));
        add(testCase(&ByteRGBImageExportImportTest::testTIFF));
        add(testCase(&ByteRGBImageExportImportTest::testBMP));
        add(testCase(&ByteRGBImageExportImportTest::testPPM));
        add(testCase(&ByteRGBImageExportImportTest::testPNM));
        add(testCase(&ByteRGBImageExportImportTest::testPNM2));
        add(testCase(&ByteRGBImageExportImportTest::testPNG));
        add(testCase(&ByteRGBImageExportImportTest::testSUN));
        add(testCase(&ByteRGBImageExportImportTest::testVIFF1));
        add(testCase(&ByteRGBImageExportImportTest::testVIFF2));

#if defined(HasPNG)
        // 16-bit PNG
        add(testCase(&PNGInt16Test::testByteOrder));
#endif

        add(testCase(&CanvasSizeTest::testTIFFCanvasSize));

        add(testCase(&PositionTest::testEXRPosition));
        add(testCase(&PositionTest::testTIFFPosition));
        add(testCase(&PositionTest::testPNGPosition));

        // grayscale float images
        add(testCase(&FloatImageExportImportTest::testGIF));
        add(testCase(&FloatImageExportImportTest::testJPEG));
        add(testCase(&FloatImageExportImportTest::testPNG));
        add(testCase(&FloatImageExportImportTest::testTIFF));
        add(testCase(&FloatImageExportImportTest::testTIFFForcedRange));
        add(testCase(&FloatImageExportImportTest::testBMP));
        add(testCase(&FloatImageExportImportTest::testSUN));
        add(testCase(&FloatImageExportImportTest::testVIFF));

        // 4-band images
        add(testCase(&Vector4ExportImportTest::testJPEG));
        add(testCase(&Vector4ExportImportTest::testGIF));
        add(testCase(&Vector4ExportImportTest::testBMP));
        add(testCase(&Vector4ExportImportTest::testPNM));
        add(testCase(&Vector4ExportImportTest::testSUN));
        add(testCase(&Vector4ExportImportTest::testVIFF));
        add(testCase(&Vector4ExportImportTest::testTIFF));
        add(testCase(&Vector4ExportImportTest::testEXR));
        add(testCase(&Vector4ExportImportTest::testPNG));

        // rgb float images
        add(testCase(&FloatRGBImageExportImportTest::testJPEG));
        add(testCase(&FloatRGBImageExportImportTest::testTIFF));
        add(testCase(&FloatRGBImageExportImportTest::testTIFFForcedRange));
        add(testCase(&FloatRGBImageExportImportTest::testBMP));
        add(testCase(&FloatRGBImageExportImportTest::testSUN));
        add(testCase(&FloatRGBImageExportImportTest::testVIFF));
        add(testCase(&FloatRGBImageExportImportTest::testHDR));

        // failure tests
        add(testCase(&ImageExportImportFailureTest::testGIFExport));
        add(testCase(&ImageExportImportFailureTest::testGIFImport));
        add(testCase(&ImageExportImportFailureTest::testJPEGExport));
        add(testCase(&ImageExportImportFailureTest::testJPEGImport));
        add(testCase(&ImageExportImportFailureTest::testTIFFExport));
        add(testCase(&ImageExportImportFailureTest::testTIFFImport));
        add(testCase(&ImageExportImportFailureTest::testBMPExport));
        add(testCase(&ImageExportImportFailureTest::testBMPImport));
        add(testCase(&ImageExportImportFailureTest::testPNMExport));
        add(testCase(&ImageExportImportFailureTest::testPNMImport));
        add(testCase(&ImageExportImportFailureTest::testPNGExport));
        add(testCase(&ImageExportImportFailureTest::testPNGImport));
        add(testCase(&ImageExportImportFailureTest::testSUNExport));
        add(testCase(&ImageExportImportFailureTest::testSUNImport));
        add(testCase(&ImageExportImportFailureTest::testVIFFExport));
        add(testCase(&ImageExportImportFailureTest::testVIFFImport));
        add(testCase(&ImageExportImportFailureTest::testShapeMismatch));

        // alpha-channel tests
        // add(testCase(&GrayscaleImportExportAlphaTest::testTIFF));
        // add(testCase(&GrayscaleImportExportAlphaTest::testPNG));
        // add(testCase(&RGBImportExportAlphaTest::testTIFF));
        // add(testCase(&RGBImportExportAlphaTest::testPNG));
    }
};


int main (int argc, char ** argv)
{
    ImageImportExportTestSuite test;
    const int failed = test.run(vigra::testsToBeExecuted(argc, argv));
    std::cout << test.report() << std::endl;

    return failed != 0;
}
