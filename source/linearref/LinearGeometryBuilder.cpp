/**********************************************************************
 * $Id: ExtractLineByLocation.cpp 1938 2006-12-07 10:45:16Z strk $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2005-2006 Refractions Research Inc.
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation.
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: linearref/ExtractLineByLocation.java rev. 1.35
 *
 **********************************************************************/

#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/CoordinateArraySequence.h>
#include <geos/geom/LineString.h>
#include <geos/linearref/ExtractLineByLocation.h>
#include <geos/linearref/LengthIndexedLine.h>
#include <geos/linearref/LinearLocation.h>
#include <geos/linearref/LengthLocationMap.h>
#include <geos/linearref/LengthIndexOfPoint.h>
#include <geos/linearref/LinearGeometryBuilder.h>
#include <geos/util/IllegalArgumentException.h>

using namespace std;

using namespace geos::geom;

namespace geos
{

namespace linearref   // geos.linearref
{

/* public */
LinearGeometryBuilder::LinearGeometryBuilder(const GeometryFactory* geomFact) :
		geomFact(geomFact),
		ignoreInvalidLines(false),
		fixInvalidLines(false),
		coordList(0) {}

/* public */
void
LinearGeometryBuilder::setIgnoreInvalidLines(bool ignoreInvalidLines)
{
	this->ignoreInvalidLines = ignoreInvalidLines;
}

/* public */
void
LinearGeometryBuilder::setFixInvalidLines(bool fixInvalidLines)
{
	this->fixInvalidLines = fixInvalidLines;
}

/* public */
void
LinearGeometryBuilder::add(const Coordinate& pt)
{
	add(pt, true);
}

/* public */
void
LinearGeometryBuilder::add(const Coordinate& pt, bool allowRepeatedPoints)
{
	if (!coordList)
		coordList = new CoordinateArraySequence();
	coordList->add(pt, allowRepeatedPoints);
	lastPt = pt;
}

/* public */
Coordinate
LinearGeometryBuilder::getLastCoordinate() const
{
	return lastPt;
}

/* public */
void
LinearGeometryBuilder::endLine()
{
	if (!coordList)
	{
		return;
	}
	if (coordList->size() < 2)
	{
		if (ignoreInvalidLines)
		{
			if (coordList)
			{
				delete coordList;
				coordList = 0;
			}
			return;
		}
		else if (fixInvalidLines)
		{
			add((*coordList)[0]);
		}
	}

	LineString* line = 0;
	try
	{
		line = geomFact->createLineString(coordList);
	}
	catch (util::IllegalArgumentException ex)
	{
		// exception is due to too few points in line.
		// only propagate if not ignoring short lines
		if (! ignoreInvalidLines)
			throw ex;
	}

	if (line) lines.push_back(line);
	coordList = 0;
}

/* public */
Geometry*
LinearGeometryBuilder::getGeometry()
{
	// end last line in case it was not done by user
	endLine();

	// NOTE: lines elements are cloned
	return geomFact->buildGeometry(lines);
}

} // namespace geos.linearref
} // namespace geos
