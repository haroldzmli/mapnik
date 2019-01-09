/*****************************************************************************
 *
 * This file is part of Mapnik (c++ mapping toolkit)
 *
 * Copyright (C) 2017 Artem Pavlenko
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *****************************************************************************/

#include "gdal_datasource.hpp"
#include "gdal_featureset.hpp"

// mapnik
#include <mapnik/debug.hpp>
#include <mapnik/boolean.hpp>
#include <mapnik/geom_util.hpp>
#include <mapnik/timer.hpp>
#include <mapnik/value/types.hpp>

#include <gdal_version.h>

// boost
#include <boost/format.hpp>
#include <boost/make_shared.hpp>
#include <boost/tokenizer.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_adapt_adt_attributes.hpp>
#include <mutex>

using mapnik::datasource;
using mapnik::parameters;

DATASOURCE_PLUGIN(gdal_datasource)

using mapnik::box2d;
using mapnik::coord2d;
using mapnik::query;
using mapnik::featureset_ptr;
using mapnik::layer_descriptor;
using mapnik::datasource_exception;

static std::once_flag once_flag;

extern "C" MAPNIK_EXP void on_plugin_load()
{
    // initialize gdal formats
    std::call_once(once_flag,[](){
        GDALAllRegister();
    });
}

gdal_datasource::gdal_datasource(parameters const& params)
    : datasource(params),
      dataset_(nullptr, &GDALClose),
      desc_(gdal_datasource::name(), "utf-8"),
      nodata_value_(params.get<double>("nodata")),
      nodata_tolerance_(*params.get<double>("nodata_tolerance",1e-12))
{
    MAPNIK_LOG_DEBUG(gdal) << "gdal_datasource: Initializing...";

#ifdef MAPNIK_STATS
    mapnik::progress_timer __stats__(std::clog, "gdal_datasource::init");
#endif

    boost::optional<std::string> file = params.get<std::string>("file");
    if (!file) throw datasource_exception("missing <file> parameter");

    boost::optional<std::string> base = params.get<std::string>("base");
    if (base)
    {
        dataset_name_ = *base + "/" + *file;
    }
    else
    {
        dataset_name_ = *file;
    }

    shared_dataset_ = *params.get<mapnik::boolean_type>("shared", false);
    band_ = *params.get<mapnik::value_integer>("band", -1);
    
    // Maximum memory limitation for image will be simply based on the maximum
    // area we allow for an image. The true memory footprint therefore will vary based
    // on the type of imagery that exists. This is not the maximum size of an image
    // on disk but rather the maximum size we will load into mapnik from GDAL.
    // max_im_area based on 50 mb limit for RGBA
    max_image_area_ = *params.get<mapnik::value_integer>("max_image_area", (50*1024*1024) / 4);

#if GDAL_VERSION_NUM >= 1600
    if (shared_dataset_)
    {
        auto ds = GDALOpenShared(dataset_name_.c_str(), GA_ReadOnly);
        dataset_.reset(static_cast<GDALDataset*>(ds));
    }
    else
#endif
    {
        auto ds = GDALOpen(dataset_name_.c_str(), GA_ReadOnly);
        dataset_.reset(static_cast<GDALDataset*>(ds));
    }

    if (! dataset_)
    {
        throw datasource_exception(CPLGetLastErrorMsg());
    }

    MAPNIK_LOG_DEBUG(gdal) << "gdal_featureset: opened Dataset=" << dataset_.get();

    nbands_ = dataset_->GetRasterCount();
    width_ = dataset_->GetRasterXSize();
    height_ = dataset_->GetRasterYSize();
    desc_.add_descriptor(mapnik::attribute_descriptor("nodata", mapnik::Double));

    double tr[6];
    bool bbox_override = false;
    boost::optional<std::string> bbox_s = params.get<std::string>("extent");
    if (bbox_s)
    {
        MAPNIK_LOG_DEBUG(gdal) << "gdal_datasource: BBox Parameter=" << *bbox_s;

        bbox_override = extent_.from_string(*bbox_s);
        if (! bbox_override)
        {
            throw datasource_exception("GDAL Plugin: bbox parameter '" + *bbox_s + "' invalid");
        }
    }

    if (bbox_override)
    {
        tr[0] = extent_.minx();
        tr[1] = extent_.width() / (double)width_;
        tr[2] = 0;
        tr[3] = extent_.maxy();
        tr[4] = 0;
        tr[5] = -extent_.height() / (double)height_;
        MAPNIK_LOG_DEBUG(gdal) << "gdal_datasource extent override gives Geotransform="
                               << tr[0] << "," << tr[1] << ","
                               << tr[2] << "," << tr[3] << ","
                               << tr[4] << "," << tr[5];
    }
    else
    {
        if (dataset_->GetGeoTransform(tr) != CPLE_None)
        {
            MAPNIK_LOG_DEBUG(gdal) << "gdal_datasource GetGeotransform failure gives="
                                   << tr[0] << "," << tr[1] << ","
                                   << tr[2] << "," << tr[3] << ","
                                   << tr[4] << "," << tr[5];
        }
        else
        {
            MAPNIK_LOG_DEBUG(gdal) << "gdal_datasource Geotransform="
                                   << tr[0] << "," << tr[1] << ","
                                   << tr[2] << "," << tr[3] << ","
                                   << tr[4] << "," << tr[5];
        }
    }

    // TODO - We should throw for true non-north up images, but the check
    // below is clearly too restrictive.
    // https://github.com/mapnik/mapnik/issues/970
    /*
      if (tr[2] != 0 || tr[4] != 0)
      {
      throw datasource_exception("GDAL Plugin: only 'north up' images are supported");
      }
    */

    dx_ = tr[1];
    dy_ = tr[5];

    if (! bbox_override)
    {
        double x0 = tr[0];
        double y0 = tr[3];
        double x1 = tr[0] + width_ * dx_ + height_ *tr[2];
        double y1 = tr[3] + width_ *tr[4] + height_ * dy_;

        /*
          double x0 = tr[0] + (height_) * tr[2]; // minx
          double y0 = tr[3] + (height_) * tr[5]; // miny

          double x1 = tr[0] + (width_) * tr[1]; // maxx
          double y1 = tr[3] + (width_) * tr[4]; // maxy
        */

        extent_.init(x0, y0, x1, y1);
    }
	GDALDataType band_type = dataset_->GetRasterBand(1)->GetRasterDataType();
    boost::optional<std::string> rgbaStr = params.get<std::string>("rgba");
    int count = 0;
    if (rgbaStr)
    {
		boost::char_separator<char> sep(", ");
		double d;
		boost::tokenizer<boost::char_separator<char> > tok(*rgbaStr, sep);

		for (boost::tokenizer<boost::char_separator<char> >::iterator beg = tok.begin();
				beg != tok.end(); ++beg)
		{
			std::string item = mapnik::util::trim_copy(*beg);
			// note: we intentionally do not use mapnik::util::conversions::string2double
			// here to ensure that shapeindex can statically compile mapnik::box2d without
			// needing to link to libmapnik
			boost::spirit::qi::double_type double_;
			boost::spirit::ascii::space_type space;
			std::string::const_iterator str_beg = item.begin();
			std::string::const_iterator str_end = item.end();
			bool r = boost::spirit::qi::phrase_parse(str_beg,
					str_end,
					double_,
					space,
					d);
			if (!(r && (str_beg == str_end)))
			{
				break;
			}
			int bandIndex = static_cast<int>(d);
			bandInfo(band_type,dataset_->GetRasterBand(bandIndex),count);
			count++;
		}
   }
   else
   {
	   int nband = 0;
	   if(nbands_ > 3)
		   nband = 3;
	   else
		   nband = nbands_;
	   for (int i = 0; i < nbands_; ++i) {
		   bandInfo(band_type,dataset_->GetRasterBand(i+1),count);
		   count++;
		}
	}

    MAPNIK_LOG_DEBUG(gdal) << "gdal_datasource: Raster Size=" << width_ << "," << height_;
    MAPNIK_LOG_DEBUG(gdal) << "gdal_datasource: Raster Extent=" << extent_;

}

gdal_datasource::~gdal_datasource()
{
    MAPNIK_LOG_DEBUG(gdal) << "gdal_featureset: Closing Dataset=" << dataset_.get();
}

datasource::datasource_t gdal_datasource::type() const
{
    return datasource::Raster;
}

const char * gdal_datasource::name()
{
    return "gdal";
}

box2d<double> gdal_datasource::envelope() const
{
    return extent_;
}

boost::optional<mapnik::datasource_geometry_t> gdal_datasource::get_geometry_type() const
{
    return boost::optional<mapnik::datasource_geometry_t>();
}

layer_descriptor gdal_datasource::get_descriptor() const
{
    return desc_;
}

featureset_ptr gdal_datasource::features(query const& q) const
{
#ifdef MAPNIK_STATS
    mapnik::progress_timer __stats__(std::clog, "gdal_datasource::features");
#endif

    return std::make_shared<gdal_featureset>(*dataset_,
                                              band_,
                                              gdal_query(q),
                                              extent_,
                                              width_,
                                              height_,
                                              nbands_,
                                              dx_,
                                              dy_,
                                              nodata_value_,
                                              nodata_tolerance_,
											  bandInfo_,
                                              max_image_area_);
}

featureset_ptr gdal_datasource::features_at_point(coord2d const& pt, double tol) const
{
#ifdef MAPNIK_STATS
    mapnik::progress_timer __stats__(std::clog, "gdal_datasource::features_at_point");
#endif

    return std::make_shared<gdal_featureset>(*dataset_,
                                              band_,
                                              gdal_query(pt),
                                              extent_,
                                              width_,
                                              height_,
                                              nbands_,
                                              dx_,
                                              dy_,
                                              nodata_value_,
                                              nodata_tolerance_,
											  bandInfo_,
                                              max_image_area_);
}
void gdal_datasource::bandInfo(GDALDataType band_type, GDALRasterBand * band,int count)
{
	int bHaveNoData;
	float fNoData = (float) band->GetNoDataValue(&bHaveNoData);
	std::pair<double,double> maxmin;
	int bGotMin, bGotMax;
	double adfMinMax[2];
	adfMinMax[0] = band->GetMinimum(&bGotMin);
	adfMinMax[1] = band->GetMaximum(&bGotMax);
	if (!(bGotMin && bGotMax))
		GDALComputeRasterMinMax(band, TRUE, adfMinMax);
	switch (band_type) {
	case GDT_UInt16: {
		unsigned short *maxMinBuff = (unsigned short *) CPLCalloc(
				sizeof(unsigned short), 1024 * 1024);
		band->RasterIO(GF_Read, 0, 0, width_, height_, maxMinBuff, 1024, 1024,
				GDT_UInt16, 0, 0, NULL);
		maxmin = histogramAccumlateMinMax(band_type, maxMinBuff, 1024, 1024, fNoData,adfMinMax);
		CPLFree(maxMinBuff);
		break;
	}
	case GDT_Int16: {
			short *maxMinBuff = (short *) CPLCalloc(
					sizeof(short), 1024 * 1024);
			band->RasterIO(GF_Read, 0, 0, width_, height_, maxMinBuff, 1024, 1024,
					GDT_Int16, 0, 0, NULL);
			maxmin = histogramAccumlateMinMax(band_type, maxMinBuff, 1024, 1024, fNoData,adfMinMax);
			CPLFree(maxMinBuff);
			break;
		}
	case GDT_Float32: {
		float *maxMinBuff = (float *) CPLCalloc(
				sizeof(float), 1024 * 1024);
		band->RasterIO(GF_Read, 0, 0, width_, height_, maxMinBuff, 1024, 1024,
				GDT_Float32, 0, 0, NULL);
		maxmin = histogramAccumlateMinMax(band_type, maxMinBuff, 1024, 1024, fNoData,adfMinMax);
		CPLFree(maxMinBuff);
		break;
	}
	}
	std::map<std::string, double> temp;
	temp.insert(std::pair<std::string, double>("max", maxmin.first));
	temp.insert(std::pair<std::string, double>("min", maxmin.second));
	temp.insert(std::pair<std::string, double>("band", band->GetBand()));
	temp.insert(std::pair<std::string, double>("nodata", fNoData));
	if (count == 0) {
		bandInfo_.insert(std::make_pair("red", temp));
	}
	if (count == 1) {
		bandInfo_.insert(std::make_pair("green", temp));
	}
	if (count == 2) {
		bandInfo_.insert(std::make_pair("blue", temp));
	}
	if (count == 3) {
		bandInfo_.insert(std::make_pair("alpha", temp));
	}
}
std::pair<double,double> gdal_datasource::histogramAccumlateMinMax(GDALDataType band_type,void *pBuf,int width,int height,double fNoData,double adfMinMax[])
{
	int length = static_cast<int>(adfMinMax[1] - adfMinMax[0]) + 1;
	if(length < 500)
		length = 500;
	int num[length];
	memset(num, 0, sizeof(num));
	//todo 1/myBinSize may nan .this is not normal
	double myBinSize = (adfMinMax[1] - adfMinMax[0]) / length;
	int nonNullCount = 0;
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			double tmpData = 0;
			switch (band_type) {
			case GDT_UInt16: {
				unsigned short v = ((unsigned short *) pBuf)[x * height + y];
				tmpData = static_cast<double>(v);
				break;
			}
			case GDT_Int16: {
				short v = ((short *) pBuf)[x * height + y];
				tmpData = static_cast<double>(v);
				break;
			}
			case GDT_Float32: {
				float v = ((float *) pBuf)[x * height + y];
				tmpData = static_cast<double>(v);
				break;
			}
			}
			if (tmpData == fNoData)
				continue;
			int myBinIndex = static_cast<int>(std::floor(tmpData - adfMinMax[0])
					/ myBinSize);
			if (myBinIndex < 0)
				myBinIndex = 0;
			if (myBinIndex > (length - 1))
				myBinIndex = length -1;
			num[myBinIndex] += 1;
			nonNullCount++;
		}
	}
	int myMincount = static_cast<int>(std::round(0.02 * nonNullCount));
	int myMaxcount = static_cast<int>(std::round(0.98 * nonNullCount));
	bool myLowerFound = false;
	int myCount = 0;
	double lowerValue = std::numeric_limits<double>::quiet_NaN();
	double upperValue = std::numeric_limits<double>::quiet_NaN();
	for (int myBin = 0; myBin < length; myBin++) {
		int myBinValue = num[myBin];
		myCount += myBinValue;
		if (!myLowerFound && myCount > myMincount) {
			lowerValue = adfMinMax[0] + myBin * myBinSize;
			myLowerFound = true;
		}
		if (myCount >= myMaxcount) {
			upperValue = adfMinMax[0] + myBin * myBinSize;
			break;
		}
	}
	//std::cout << "max" << upperValue << "min" << lowerValue << std::endl;
	return std::pair<double,double>(std::ceil(upperValue),std::floor(lowerValue));
	//return std::pair<double,double>(upperValue,lowerValue);
}

