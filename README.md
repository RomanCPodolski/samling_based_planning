SAEAuto
=======

SAE-Autonomous

Before you start to work on this project, please read the [contribution guidelines](CONTRIBUT.md).
## Dependencies:
```
libboost-system-dev

libncurses5-dev

libqgpsmm-dev

libgps-dev

libboost-thread-dev

libopencv-dev

libcv-dev

libcv2.4

gpsd


lighttpd - enable CGI, config to run as user who is a sudoer (e.g. pi).
Pages served from FrontEnd/Maps.

crond - use to run M2Mconnect.py on boot for car internet connection.

cpplint

protobuf (libprotobuf-dev, libprotobuf-compiler, libprotobuf-c-compiler)

protoc (libprotoc-dev)

glog

gtest*

Perl modules: JSON (libjson-perl), Template (libtemplate-perl)
```
## Compile
To compile Control.cpp:
```
python3 waf configure build or python waf configure build
```
*Set up gtest:
```
sudo apt-get install libgtest-dev cmake
cd /usr/src/gtest
sudo cmake CMakeLists.txt
sudo make
sudo cp *.a /usr/lib

also:
see https://www.eriksmistad.no/getting-started-with-google-test-on-ubuntu/

```
