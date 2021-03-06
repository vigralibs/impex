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

#ifndef VIGRA2_IMAGEIO_IMPEXBASE_HXX
#define VIGRA2_IMAGEIO_IMPEXBASE_HXX

#include <limits>
#include <string>
#include <type_traits>
#include <utility>

#include <vigra2/sized_int.hxx>
#include <vigra2/utilities.hxx>


namespace vigra
{
    typedef enum
    {
        UNSIGNED_INT_8,
        UNSIGNED_INT_16,
        UNSIGNED_INT_32,
        SIGNED_INT_16,
        SIGNED_INT_32,
        IEEE_FLOAT_32,
        IEEE_FLOAT_64
    } pixel_t;


    namespace detail
    {

        template <class VALUETYPE>
        class FindMinMax
        {
           public:

                /** the functor's argument type
                */
            typedef VALUETYPE argument_type;

                /** the functor's result type
                */
            typedef VALUETYPE result_type;

                /** \deprecated use argument_type
                */
            typedef VALUETYPE value_type;

                /** init min and max
                */
            FindMinMax()
            : min( std::numeric_limits<value_type>::max() ),
              max( std::numeric_limits<value_type>::min() ),
              count(0)
            {}

                /** (re-)init functor (clear min, max)
                */
            void reset()
            {
                count = 0;
            }

                /** update min and max
                */
            void operator()(argument_type const & v)
            {
                if(count)
                {
                    if(v < min) min = v;
                    if(max < v) max = v;
                }
                else
                {
                    min = v;
                    max = v;
                }
                ++count;
            }

#if 0
                /** update min and max with components of RGBValue<VALUETYPE>
                */
            void operator()(RGBValue<VALUETYPE> const & v)
            {
                operator()(v.red());
                operator()(v.green());
                operator()(v.blue());
            }

                /** merge two statistics
                */
            void operator()(FindMinMax const & v)
            {
                if(v.count)
                {
                    if(count)
                    {
                        if(v.min < min) min = v.min;
                        if((this->max) < v.max) max = v.max;
                    }
                    else
                    {
                        min = v.min;
                        max = v.max;
                    }
                }
                count += v.count;
            }
#endif

                /** the current min
                */
            VALUETYPE min;

                /** the current max
                */
            VALUETYPE max;

                /** the number of values processed so far
                */
            unsigned int count;

        };


        inline static pixel_t
        pixel_t_of_string(const std::string& pixel_type)
        {
            if (pixel_type == "BILEVEL")
            {
                return UNSIGNED_INT_8;
            }
            else if (pixel_type == "UINT8")
            {
                return UNSIGNED_INT_8;
            }
            else if (pixel_type == "UINT16")
            {
                return UNSIGNED_INT_16;
            }
            else if (pixel_type == "UINT32")
            {
                return UNSIGNED_INT_32;
            }
            else if (pixel_type == "INT16")
            {
                return SIGNED_INT_16;
            }
            else if (pixel_type == "INT32")
            {
                return SIGNED_INT_32;
            }
            else if (pixel_type == "FLOAT")
            {
                return IEEE_FLOAT_32;
            }
            else if (pixel_type == "DOUBLE")
            {
                return IEEE_FLOAT_64;
            }
            else
            {
                vigra_fail("vigra_ext::detail::pixel_t_of_string: unknown pixel type");
                return UNSIGNED_INT_8; // NOT REACHED
            }
        }


        struct identity
        {
            template <typename T>
            T operator()(T x) const
            {
                return x;
            }
        };


        typedef std::pair<double, double> range_t;


        class linear_transform
        {
        public:
            linear_transform(const range_t& source, const range_t& destination) :
                scale_((destination.second - destination.first) / (source.second - source.first)),
                offset_(destination.first / scale_ - source.first)
            {}

            template <typename T>
            double operator()(T x) const
            {
                return scale_ * (static_cast<double>(x) + offset_);
            }

        private:
            const double scale_;
            const double offset_;
        };


        template <class Iterator, class Accessor>
        inline static range_t
        find_value_range(Iterator upper_left, Iterator lower_right, Accessor accessor,
                         /* is_scalar? */ std::true_type)
        {
            typedef typename Accessor::value_type value_type;

            FindMinMax<value_type> extrema;

            inspectImage(upper_left, lower_right, accessor, extrema);

            return range_t(static_cast<double>(extrema.min), static_cast<double>(extrema.max));
        }


        template <class Iterator, class Accessor>
        inline static range_t
        find_value_range(Iterator upper_left, Iterator lower_right, Accessor accessor,
                         /* is_scalar? */ std::false_type)
        {
            typedef typename Accessor::ElementAccessor element_accessor;
            typedef typename element_accessor::value_type value_type;

            const int number_of_bands(static_cast<int>(accessor.size(upper_left)));
            FindMinMax<value_type> extrema;

            for (int i = 0; i != number_of_bands; ++i)
            {
                element_accessor band(i, accessor);

                inspectImage(upper_left, lower_right, band, extrema);
            }

            return range_t(static_cast<double>(extrema.min), static_cast<double>(extrema.max));
        }


        template <class SourceIterator, class SourceAccessor>
        inline static range_t
        find_source_value_range(const ImageExportInfo& export_info,
                                SourceIterator upper_left, SourceIterator lower_right, SourceAccessor accessor)
        {
            if (export_info.getFromMin() < export_info.getFromMax())
            {
                return range_t(export_info.getFromMin(), export_info.getFromMax());
            }
            else
            {
                typedef typename SourceAccessor::value_type SourceValueType;
                typedef typename NumericTraits<SourceValueType>::isScalar is_scalar;

                const range_t range(find_value_range(upper_left, lower_right, accessor, is_scalar()));

                if (range.first < range.second)
                {
                    return range_t(range.first, range.second);
                }
                else
                {
                    return range_t(range.first, range.first + 1.0);
                }
            }
        }

        template <class T, int N>
        inline static range_t
        find_source_value_range(const ImageExportInfo& export_info, ArrayViewND<N,T> view)
        {
            if (export_info.getFromMin() < export_info.getFromMax())
            {
                return range_t(export_info.getFromMin(), export_info.getFromMax());
            }
            else
            {
                auto minmax = view.minmax();
                const range_t range(minmax[0],minmax[1]);

                if (range.first < range.second)
                {
                    return range_t(range.first, range.second);
                }
                else
                {
                    return range_t(range.first, range.first + 1.0);
                }
            }
        }


        template <typename T>
        inline static range_t
        find_destination_value_range(const ImageExportInfo& export_info)
        {
            if (export_info.getToMin() < export_info.getToMax())
            {
                return range_t(export_info.getToMin(), export_info.getToMax());
            }
            else
            {
                return range_t(static_cast<double>(NumericTraits<T>::min()),
                               static_cast<double>(NumericTraits<T>::max()));
            }
        }


        inline static range_t
        find_destination_value_range(const ImageExportInfo& export_info, pixel_t pixel_type)
        {
            switch (pixel_type)
            {
            case UNSIGNED_INT_8: return find_destination_value_range<uint8_t>(export_info);
            case UNSIGNED_INT_16: return find_destination_value_range<uint16_t>(export_info);
            case UNSIGNED_INT_32: return find_destination_value_range<uint32_t>(export_info);
            case SIGNED_INT_16: return find_destination_value_range<int16_t>(export_info);
            case SIGNED_INT_32: return find_destination_value_range<int32_t>(export_info);
            case IEEE_FLOAT_32: return find_destination_value_range<float>(export_info);
            case IEEE_FLOAT_64: return find_destination_value_range<double>(export_info);
            default:
                vigra_fail("vigra_ext::detail::find_destination_value_range: not reached");
                return range_t(0.0, 0.0); // NOT REACHED
            }
        }
    } // end namespace detail
} // end namespace vigra


#endif // VIGRA_IMPEXBASE_HXX
