cd c:\
md TxSSA-VS2005-32-shared-build
cd TxSSA-VS2005-32-shared-build
call C:\txssa\cmake\windows-batch\TxSSA-VS2005-32-shared.bat
"C:\Program Files (x86)\Microsoft Visual Studio 8\Common7\IDE\devenv.exe" TxSSA.sln /build Release /project PACKAGE.vcproj
move TxSSA_1.0_Windows.zip TxSSA_1.0_Windows-VS2005-32bit.zip

cd c:\
md TxSSA-VS2005-64-shared-build
cd TxSSA-VS2005-64-shared-build
call C:\txssa\cmake\windows-batch\TxSSA-VS2005-64-shared.bat
"C:\Program Files (x86)\Microsoft Visual Studio 8\Common7\IDE\devenv.exe" TxSSA.sln /build Release /project PACKAGE.vcproj
move TxSSA_1.0_Windows.zip TxSSA_1.0_Windows-VS2005-64bit.zip

cd c:\
md TxSSA-VS2012-64-shared-build
cd TxSSA-VS2012-64-shared-build
call C:\txssa\cmake\windows-batch\TxSSA-VS2012-64-shared.bat
"C:\Program Files (x86)\Microsoft Visual Studio 11.0\Common7\IDE\devenv.exe" TxSSA.sln /build Release /project PACKAGE.vcxproj
move TxSSA_1.0_Windows.zip TxSSA_1.0_Windows-VS2012-64bit.zip
