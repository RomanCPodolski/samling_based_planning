#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Roman C. Podolski <mailto:roman.podolski@tum.de>

APPNAME='sampling-based-planning'
VERSION='0.0.1'

top = '.'
out = 'build'

def options(opt):
    opt.add_option("-d","--doc", action='store_true', dest="doc", help='generate documentation')
    opt.load('compiler_c compiler_cxx boost cpplint')

    opt.recurse('src test')


def configure(conf):
    conf.load('compiler_c compiler_cxx boost cpplint clang_compilation_database doxygen')
    if not conf.env.DOXYGEN:
        conf.fatal('doxygen is required, please install it')

    conf.env.CPPLINT_BREAK = 6  # build will not break if linter errors exist

    conf.env.CPPLINT_ROOT='sampling-based-planning'

    conf.env.CPPLINT_FILTERS = ','.join((
        '-readability/streams', # glog
    ))

    conf.recurse('src test')

def build(bld):
    if bld.options.doc:
        bld(
            features='doxygen',
            doxyfile='Doxyfile',
            install_path="doc"
        )

    bld.recurse('src test')

def dist(ctx):
    ctx.excl ='**/.waf* **/*~ **/*.pyc **/*.swp **/.lock-w* build/* log RunFiles/* .ropeproject/*'
