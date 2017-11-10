/*
 * map.cc
 * Copyright (C) 2017 Roman C. Podolski <maito:roman.podolski@tum.de>
 *
 * Distributed under terms of the MIT license.
 */

#include "planner/map.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <glog/logging.h>
#include <iostream>

DEFINE_string(map, "data/map.wyp", "name of map");

namespace planner {
using boost::lexical_cast;
using boost::bad_lexical_cast;
using boost::iends_with;
using boost::erase_all;

Map::Map(std::string path) { load(path); }
// loads the whole map onto memory
void Map::parse_from_stream(std::istream *is) {
  std::string s = "";
  std::vector<std::string> line;
  bool datum_set = false;

  for (int i = 0; std::getline(*is, s); i++) {
    erase_all(s, " ");
    boost::algorithm::split(line, s, boost::algorithm::is_any_of(", "));

    try {
      if (s[0] == 'D') {
        LOG_IF(FATAL, i != 0) << "First line is expected to be a Datum entry";
        LOG_IF(FATAL, datum_set) << "File contains several Datum entry";
        _datum = std::pair<double, double>(lexical_cast<double>(line[1]),
                                           lexical_cast<double>(line[2]));
        datum_set = true;
      } else if (s[0] == 'P') {
        // do nothing
      } else {
        std::pair<double, double> m(lexical_cast<double>(line[1]),
                                    lexical_cast<double>(line[2]));
        if (s[0] == 'F') {
          add_fencepost(m);
        } else {
          lexical_cast<int>(line[0]); // sanity check
          add_waypoint(m);
        }
      }
    } catch (const bad_lexical_cast &) {
      LOG(FATAL) << "Fail to parse map at line: " << s;
    }
  }
}
void Map::load(std::string path) {
  // use the csv
  CHECK(iends_with(path, ".wyp") || iends_with(path, ".csv"))
      << "only .wyp or .csv files are accepted";
  std::filebuf fb;
  fb.open(path, std::ios::in);
  std::istream is(&fb);

  CHECK(is) << "Could not open map " << path;
  parse_from_stream(&is);
  fb.close();
}

void Map::save(std::string path) { LOG(FATAL) << "uninplemented"; }
} // namespace planner
