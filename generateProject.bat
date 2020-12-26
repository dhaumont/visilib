set EMBREE_INSTALL_DIR=c:\Program Files\Intel\Embree3 x64
set CGAL_DIR=C:\dev\CGAL-4.13.1
set BOOST_INCLUDEDIR=C:\code\boost\boost_1_70_0\
set BOOST_ROOT=C:\code\boost\boost_1_70_0\
set LEDA_DIR=C:\Algorithmic Solutions\LEDA-6.3-free-win64-msc10-std

mkdir build
cd build
Cmake ..
cd ..

xcopy /y demo\external\freeglut\x64\freeglut.dll build\demo\Debug\
xcopy /y demo\external\freeglut\x64\freeglut.dll build\demo\Release\

if exist "%EMBREE_INSTALL_DIR%\bin\embree3.dll"  xcopy /y "%EMBREE_INSTALL_DIR%\bin\*.dll" build\demo\Debug\
if exist "%EMBREE_INSTALL_DIR%\bin\embree3.dll"  xcopy /y "%EMBREE_INSTALL_DIR%\bin\*.dll" build\demo\Release\

if exist "%LEDA_DIR%\leda_md.dll"  xcopy /y "%LEDA_DIR%\*.dll" build\demo\Debug\
if exist "%LEDA_DIR%\leda_md.dll"  xcopy /y "%LEDA_DIR%\*.dll" build\demo\Release\