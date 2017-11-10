#ifndef TEX_H
#define TEX_H

#include <boost/format.hpp>
#include <iostream>
#include "planner/baseframe.h"
#include "planner/maneuver.h"
#include "planner/planner.h"

namespace planner {
namespace io {

class beamer_mapper {
 private:
  std::ostream &_os;
  std::string _axis_props;

  void write_header() {
    _os << "% Autor: Roman C. Podolski <roman.podolski@tum.de>\n";
    _os << "% this is a generated file, please do not alter\n\n";
    _os << "\\documentclass{beamer}\n";
    _os << "\\setbeamertemplate{navigation symbols}{}%\n";
    _os << "\\usepackage{tikz}\n";
    _os << "\\usetikzlibrary{decorations.markings,calc}\n";
    _os << "\\usepackage{pgfplots}\n";
    _os << "\\usepackage{pgfplotstable}\n";
    _os << "\\usepackage[svgnames]{xcolor}\n";
    _os << "% Settings for pgfplots\n";
    _os << "\\pgfplotsset{compat=1.9}\n";
    _os << "\\pgfplotsset{cycle list={black}}\n";
    _os << "\\begin{document}\n";
    _os << "\\begin{tikzpicture}[scale=1.0]\n";
    _os << boost::format("\\begin{axis}[%s]\n") % _axis_props;
  }

  void write_footer() {
    _os << "\\end{axis}\n";
    _os << "\\end{tikzpicture}\n";
    _os << "\\end{frame}\n";
    _os << "\\end{document}\n";
  }

 public:
  beamer_mapper(std::ostream &os, std::string axis_props = "axis equal")
      : _os(os), _axis_props{axis_props} {
    write_header();
  }
  virtual ~beamer_mapper() { write_footer(); }

  void map(Baseframe b, std::string props = "color=blue",
           const double ds = 1.0) {
    _os << boost::format("\\addplot[%s] coordinates {\n") % props;
    for (double s = 0.1; s < b.length(); s += ds) {
      auto p = b.P(s);
      _os << boost::format("(%d,%d)\n") % p.x() % p.y();
    }
    _os << "};\n";
  }

  void map(Maneuver m, std::string props = "color=black") {
    _os << boost::format("\\addplot[%s] coordinates {\n") % props;
    for (auto point : m.path()) {
      _os << boost::format("(%d,%d)\n") % point.x() % point.y();
    }
    _os << "};\n";
  }

  void map(std::vector<std::pair<double, double>> wp,
           std::string props = "color=blue,mark=x") {
    _os << boost::format("\\addplot[only marks, %s] coordinates {\n") % props;
    for (const auto &p : wp) {
      _os << boost::format("(%d,%d)\n") % p.first % p.second;
    }
    _os << "};\n";
  }

  void map(polygon p, std::string props = "") {
    _os << boost::format("\\addplot[%s] coordinates {\n") % props;
    for (auto point : p.outer()) {
      _os << boost::format("(%d,%d)\n") % point.x() % point.y();
    }
    _os << "} --cycle;\n";
  }

  void map(const pose &pose, std::string props = "") {
    const double l = 3, w = 2;
    polygon car, foo, baa;
    car.outer().push_back(point(-l / 2, -w / 2));
    car.outer().push_back(point(+l / 2, -w / 2));
    car.outer().push_back(point(+l / 2, +w / 2));
    car.outer().push_back(point(-l / 2, +w / 2));
    car.outer().push_back(point(-l / 2, -w / 2));  // close the polygon

    bg::strategy::transform::rotate_transformer<planner::bg::radian, double, 2,
                                                2>
        rotate(-pose.theta);
    bg::strategy::transform::translate_transformer<double, 2, 2> translate(
        pose.x(), pose.y());
    bg::transform(car, foo, rotate);     // rotate first...
    bg::transform(foo, baa, translate);  // ...then translate
    map(baa, props);
  }

  void map(box b, std::string props = "") {
    polygon tmp;
    bg::convert(b, tmp);
    map(tmp, props);
  }

  void next_frame() {
    _os << "\\end{axis}\n\\end{tikzpicture}\n\\end{frame}\n\\begin{frame}\n";
    _os << "\\begin{tikzpicture}[scale=1.0]\n";
    _os << boost::format("\\begin{axis}[%s]\n") % _axis_props;
  }
};
}  // namespace io
}  // namespace planner

#endif /* TEX_H */
