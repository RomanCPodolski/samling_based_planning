#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <utility>
#include <vector>
#include "planner/planner.h"

using std::milli;
using boost::format;
using boost::math::constants::sixth_pi;
using boost::math::constants::two_pi;

planner::polygon get_cone(planner::point pos, double cone_radius = 0.3) {
  planner::polygon p;
  for (double i = 0; i < two_pi<double>(); i += sixth_pi<double>()) {
    double x = pos.x() + cone_radius * cos(i);
    double y = pos.y() + cone_radius * sin(i);
    p.outer().push_back(planner::point(x, y));
  }
  return p;
}

double bench_gain(int samples, int length) {
  FLAGS_min_manouver_length = length;

  double L = 1000;
  std::vector<std::pair<double, double>> wp = {};

  for (double i = 0; i < L; i++) {
    wp.push_back(std::make_pair(i, 0));
  }

  double s = 0;

  planner::Planner p(wp);
  planner::pose pos{p.baseframe().P(s), p.baseframe().theta(s), s, 0.0};

  std::vector<planner::polygon> obstacles;

  for (int i = 0; i < p.baseframe().length(); i += 5) {
    obstacles.push_back(get_cone(planner::point(i, 3.5)));
    obstacles.push_back(get_cone(planner::point(i, -3.5)));
  }

  double time_milli = 0.0;
  int idx;

  for (idx = 0; idx <= samples && s <= p.baseframe().length() - 1; idx++) {
    auto start = std::chrono::steady_clock::now();
    auto m_opt = p.update(pos.position, pos.theta, 5, 1, obstacles);
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    time_milli += std::chrono::duration<double, milli>(diff).count();
    pos = m_opt.path()[1];
    s = m_opt.path()[1].s;
  }
  double time_milli_mean = time_milli / idx;
  return time_milli_mean;
}

int main(int argc, char *argv[]) {
  std::ofstream os("bench_gain.tex");

  os << "% Autor: Roman C. Podolski <roman.podolski@tum.de>\n";
  os << "% this is a generated file, please do not alter\n\n";
  os << "\\documentclass{standalone}[preview]\n";
  os << "\\usepackage{tikz}\n";
  os << "\\usetikzlibrary{decorations.markings,calc}\n";
  os << "\\usepackage{pgfplots}\n";
  os << "\\usepackage{pgfplotstable}\n";
  os << "\\usepackage[svgnames]{xcolor}\n";
  os << "% Settings for pgfplots\n";
  os << "\\pgfplotsset{compat=1.9}\n";
  os << "\\pgfplotsset{cycle list={black}}\n";
  os << "\\begin{document}\n";
  os << "\\begin{tikzpicture}[scale=1.0]\n";
  os << boost::format("\\begin{axis}[%s]\n") % "";
  os << boost::format("\\addplot[%s] coordinates {\n") % "";
  for (double i = 0.0; i < 20; i += 0.1) {
    double time = bench_gain(100, i);
    std::cout << i << ":" << time << '\n';
    os << boost::format("(%d,%d)\n") % i % time;
  }
  os << "};";
  os << "\\end{axis}\n";
  os << "\\end{tikzpicture}\n";
  os << "\\end{document}\n";
}
