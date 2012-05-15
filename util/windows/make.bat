@echo off

rem This file will run several abinit built-in tests in a DOS/Windows box
rem in the same manner than "make tests" under unix.
rem The "run-basic-tests.pl" perl script is called; this one will start abinit
rem and periodically display the status file till the end of the process.
rem
rem Several WIN32 binary versions of Perl for Windows exist, among them:
rem     Microsoft Windows NT Resource Kit (C)
rem	binaries downloadable from www.cpan.org/ports/
rem WARNING ! perl.exe MUST be accessible through the DOS PATH

rem Copyright (C) 1999-2012 ABINIT group (LSi)
rem This file is distributed under the terms of the
rem GNU General Public License, see ~ABINIT/COPYING
rem or http://www.gnu.org/copyleft/gpl.txt .
rem For the initials of contributors, see ~ABINIT/Infos/contributors .


if not .%1 == .help goto :T1?
echo usage is: make [ tests ! test1 ! test2 ! test3 ! test4 ! test5 ]
echo           make [ test12 ! test123 ! test124 ! test125 ]   
goto :EOF
:T1?
rem change to internal tests directory
echo cd tests
cd tests
if not .%1 == .test1 goto :T2?
perl Scripts\run-basic-tests.pl Windows 1
cd ..
goto :EOF
:T2?
if not .%1 == .test2 goto :T3?
perl Scripts\run-basic-tests.pl Windows 2
cd ..
goto :EOF
:T3?
if not .%1 == .test3 goto :T4?
perl Scripts\run-basic-tests.pl Windows 3
cd ..
goto :EOF
:T4?
if not .%1 == .test4 goto :T5?
perl Scripts\run-basic-tests.pl Windows 4
cd ..
goto :EOF
:T5?
if not .%1 == .test5 goto :T12?
perl Scripts\run-basic-tests.pl Windows 5
cd ..
goto :EOF
:T12?
if not .%1 == .test12 goto :T123?
perl Scripts\run-basic-tests.pl Windows 1
perl Scripts\run-basic-tests.pl Windows 2
cd ..
goto :EOF
:T123?
if not .%1 == .test123 goto :T124?
perl Scripts\run-basic-tests.pl Windows 1
perl Scripts\run-basic-tests.pl Windows 2
perl Scripts\run-basic-tests.pl Windows 3
cd ..
goto :EOF
:T124?
if not .%1 == .test124 goto :T125?
perl Scripts\run-basic-tests.pl Windows 1
perl Scripts\run-basic-tests.pl Windows 2
perl Scripts\run-basic-tests.pl Windows 4
cd ..
goto :EOF
:T125?
if not .%1 == .test125 goto :Tous?
perl Scripts\run-basic-tests.pl Windows 1
perl Scripts\run-basic-tests.pl Windows 2
perl Scripts\run-basic-tests.pl Windows 5
cd ..
goto :EOF
:Tous?
if not .%1 == .tests goto :Error
perl Scripts\run-basic-tests.pl Windows 1
perl Scripts\run-basic-tests.pl Windows 2
perl Scripts\run-basic-tests.pl Windows 3
perl Scripts\run-basic-tests.pl Windows 4
perl Scripts\run-basic-tests.pl Windows 5
cd ..
goto :EOF
:Error
echo invalid option: %1
set ERRORLEVEL=8
cd ..
rem
:EOF
